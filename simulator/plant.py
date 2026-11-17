"""
The LBNH complex, end to end
============================

LNGP -> BOGC -> NRHR -> HLS, with the rich-LNG return loop closed.

Rigorous: thermodynamics, every flash, the MCHE composite curves and
pinch, the propane and mixed-refrigerant compressors, the nitrogen
rejection column, the helium separators, and the helium liquefier.

Estimated (and flagged in open_items.md): auxiliary refrigeration inside
NRHR - the column condenser and the cold-box deficit - costed at a fixed
fraction of Carnot instead of simulating a nitrogen loop.
"""
import json
import time

import numpy as np
import CoolProp.CoolProp as CP

import design_basis as DB
from props import PKG, C, vec, feed_kmolh, KMOLH, subpackage
from flash import flash_TP, flash_PH, bubble_T, dew_T, mixture_H
from lngp import liquefaction, compress, h_molar
from nru import nru_with_condenser
from hls import collins

K0 = 273.15
T_AMB = 313.15
EXERGY_EFF = 0.30           # cryogenic refrigeration, fraction of Carnot
NRU = ["He", "CH4", "N2"]
IN = [C.index(c) for c in NRU]


def refrig_power(Q, T):
    """Shaft power to remove Q watts at T kelvin, at EXERGY_EFF of Carnot."""
    return abs(Q) * (T_AMB / T - 1.0) / EXERGY_EFF


def to9(n3):
    out = np.zeros(len(C))
    out[IN] = n3
    return out


def run(mr=None, P_lo=None, T_cond_C=-156.0, boilup=0.7424, verbose=True,
        max_recycle=6, nrhr="cycle", bog_cold=False, accel="direct",
        rtol=1e-4, mr_file="mr_best.json", hls_file="hls_best.json"):
    """nrhr: "cycle" simulates the N2 + He refrigeration of the cold box
    (coldbox.py, parameters from optimize_nrhr.py); "estimate" is the
    original 30 %-of-Carnot shortcut that went into the sealed
    predictions. bog_cold: recover the tank boil-off gas's cold against
    the compressed BOG (a post-comparison design improvement)."""
    t0 = time.time()
    say = print if verbose else (lambda *a, **k: None)
    best = json.load(open(mr_file))
    mr = mr or best["mr"]
    P_lo = P_lo or best["Plo"] * 1e5
    hls_p = json.load(open(hls_file))

    R = dict(inputs=dict(mr=mr, P_lo_bar=P_lo / 1e5, T_cond_C=T_cond_C,
                         boilup=boilup))

    # ---------------- LNGP ------------------------------------------
    ng = vec(feed_kmolh()) * KMOLH
    L = liquefaction(ng, vec(mr), P_lo)
    say(f"LNGP  {L['W']/1e6:.1f} MW  dTmin {L['dTmin']:.2f} K  [{time.time()-t0:.0f}s]")
    lng_hp = L["lng"]
    R["lngp"] = dict(W_MW=L["W"] / 1e6, W_mr_MW=L["W_mr"] / 1e6,
                     W_c3_MW=L["W_c3"] / 1e6, dTmin=L["dTmin"],
                     F_mr_kmolh=L["F_mr"] / KMOLH,
                     hhc_kmolh=L["hhc"].sum() / KMOLH,
                     hhc_comp={c: float(L["hhc"][i] / KMOLH) for i, c in enumerate(C)},
                     mche_MW=L["Q"] / 1e6)

    # ---------------- end flash ------------------------------------
    P_ng = (DB.FEED_P_BAR - 3.0 - 2.0) * 1e5
    z = lng_hp / lng_hp.sum()
    h_in = h_molar(DB.MCHE_OUTLET_T_C + K0, P_ng, z)
    P_t = DB.END_FLASH_P_BAR * 1e5
    T_ef, b, x, y = flash_PH(PKG, P_t, z, h_in, lo=80, hi=200)
    bog_flash = b * lng_hp.sum() * y
    lng_liq = (1 - b) * lng_hp.sum() * x
    say(f"end flash {T_ef-K0:.1f} C  BOG {bog_flash.sum()/KMOLH:.0f} kmol/h")

    # ---------------- recycle: tank BOG + rich LNG return ------------
    rich_ret = np.zeros(len(C))
    rich_flash_v = np.zeros(len(C))
    st = None
    pk3 = subpackage(NRU)
    for it in range(max_recycle):
        tank_liq = lng_liq + rich_ret
        # tank boil-off from heat leak: bubble-point vapour of the tank liquid
        Tb, xb, yb = bubble_T(PKG, P_t, tank_liq / tank_liq.sum(), lo=90, hi=200)
        tank_bog = yb * tank_liq.sum() * DB.TANK_BOG_PER_HOUR
        bog = bog_flash + tank_bog + rich_flash_v
        bog3 = bog[IN]
        lost = bog.sum() - bog3.sum()             # C2+ traces lumped to CH4
        bog3[1] += lost
        # compress to NRU pressure
        P_nru = DB.BOG_COMPRESSOR_P_BAR * 1e5
        z3 = bog3 / bog3.sum()
        W_bog, Qc_bog, _ = compress(to9(z3), bog3.sum(), 303.15, P_t, P_nru, 3)
        try:
            Td, *_ = dew_T(pk3, P_nru, z3)
        except ValueError:
            # fallback for unusual feeds: bisect on the vapour fraction
            lo_, hi_ = 90.0, 220.0
            for _ in range(40):
                mid_ = 0.5 * (lo_ + hi_)
                if flash_TP(pk3, mid_, P_nru, z3)[0] >= 0.9995:
                    hi_ = mid_
                else:
                    lo_ = mid_
            Td = hi_
        T_feed = Td - 8.0
        bf, *_ = flash_TP(pk3, T_feed, P_nru, z3)
        res = nru_with_condenser(bog3, P_nru, T_feed, bf, T_cond_C + K0,
                                 boilup, state=st)
        if res is None:
            raise RuntimeError("NRU did not converge")
        st = res
        col = res["col"].results()
        bottoms = col["l"][-1]                    # mol/s, 3 components
        # subcool the rich LNG in the cold box, then let it down to the tank
        T_sub = -158.0 + K0
        zb = bottoms / bottoms.sum()
        hsub = h_molar(T_sub, P_nru, to9(zb))
        Tr, br, xr, yr = flash_PH(PKG, P_t, to9(zb), hsub, lo=80, hi=200)
        new_ret = to9(xr[IN] * (1 - br) * bottoms.sum())
        new_rfv = to9(yr[IN] * br * bottoms.sum())
        err = (np.abs(new_ret - rich_ret).sum() + np.abs(new_rfv - rich_flash_v).sum()) / bog.sum()
        # tear-stream update: direct substitution, Wegstein or Broyden
        xk = np.concatenate([rich_ret, rich_flash_v])
        gk = np.concatenate([new_ret, new_rfv])
        if accel == "direct" or it == 0:
            xn = gk
            if accel == "broyden":
                Hinv = -np.eye(xk.size)          # first step = direct substitution
        elif accel == "wegstein":
            dx = xk - x_prev
            s_ = np.where(np.abs(dx) > 1e-12, (gk - g_prev) / np.where(np.abs(dx) > 1e-12, dx, 1), 0.0)
            q = np.clip(s_ / (s_ - 1.0 + 1e-300), -5.0, 0.0)
            xn = q * xk + (1 - q) * gk
        else:                                     # Broyden, "good" inverse update
            rk, rp = gk - xk, g_prev - x_prev
            dxv, dr = xk - x_prev, rk - rp
            den = dxv @ Hinv @ dr
            if abs(den) > 1e-300:
                Hinv = Hinv + np.outer(dxv - Hinv @ dr, dxv @ Hinv) / den
            xn = xk - Hinv @ rk
        xn = np.maximum(xn, 0.0)
        x_prev, g_prev = xk, gk
        rich_ret, rich_flash_v = xn[:len(C)], xn[len(C):]
        say(f"  recycle {it+1}: BOG {bog.sum()/KMOLH:.0f} kmol/h  rich LNG {rich_ret.sum()/KMOLH:.0f}"
            f"  err {err:.2e}  [{time.time()-t0:.0f}s]")
        R.setdefault("recycle_trace", []).append(float(err))
        if err < rtol:
            break

    R["bogc"] = dict(end_flash_T_C=T_ef - K0, end_flash_vf=float(b),
                     bog_flash_kmolh=bog_flash.sum() / KMOLH,
                     tank_bog_kmolh=tank_bog.sum() / KMOLH,
                     rich_flash_kmolh=rich_flash_v.sum() / KMOLH,
                     bog_total_kmolh=bog.sum() / KMOLH,
                     bog_comp={c: float(bog3[i] / bog3.sum()) for i, c in enumerate(NRU)},
                     W_bog_MW=W_bog / 1e6)
    vap = res["vapour"]
    R["nru"] = dict(T_cond_C=T_cond_C, boilup=boilup,
                    bottoms_kmolh=bottoms.sum() / KMOLH,
                    bottoms_N2=float(zb[2]), CH4_recovery=float(bottoms[1] / bog3[1]),
                    overhead_kmolh=vap.sum() / KMOLH,
                    overhead_He=float(vap[0] / vap.sum()),
                    Q_cond_MW=abs(res["Q_cond"]) / 1e6, Q_reb_MW=res["reb"] / 1e6,
                    T_top_C=col["T"][0] - K0, T_bot_C=col["T"][-1] - K0,
                    profile_T=[float(t - K0) for t in col["T"]],
                    profile_xN2=[float(v[2]) for v in col["x"]])

    # ---------------- helium separation ------------------------------
    stream = vap.copy()
    stages = []
    n2_liq = np.zeros(3)
    Q_hesep = 0.0
    T_prev = T_cond_C + K0
    for T in (90.0, 80.0, 70.0):
        z3 = stream / stream.sum()
        b3, x3, y3, _ = flash_TP(pk3, T, P_nru, z3)
        h_in = mixture_H(pk3, T_prev, P_nru, *flash_TP(pk3, T_prev, P_nru, z3)[:3])
        h_out = mixture_H(pk3, T, P_nru, b3, x3, y3)
        Q_hesep += stream.sum() * (h_in - h_out)
        liq = (1 - b3) * stream.sum() * x3
        n2_liq += liq
        stream = b3 * stream.sum() * y3
        stages.append(dict(T_K=T, vap_kmolh=stream.sum() / KMOLH,
                           He=float(stream[0] / stream.sum()),
                           N2=float(stream[2] / stream.sum())))
        T_prev = T
    crude = stream
    # everything the NRHR cold box has to balance, for coldbox.py
    top_v = col["v"][0]
    R["coldbox"] = dict(
        P=P_nru,
        bog=dict(F=float(bog3.sum()), z=(bog3 / bog3.sum()).tolist(),
                 T_in=DB.COOLER_T_C + K0, T_out=float(T_feed)),
        cond=dict(F=float(top_v.sum()), z=(top_v / top_v.sum()).tolist(),
                  T_in=float(col["T"][0]), T_out=T_cond_C + K0,
                  Q=float(abs(res["Q_cond"]))),
        reb=dict(T=float(col["T"][-1]), T_above=float(col["T"][-2]),
                 Q=float(res["reb"])),
        rich=dict(F=float(bottoms.sum()), z=zb.tolist(),
                  T_in=float(col["T"][-1]), T_out=T_sub),
        hesep=[], n2liq=[])
    s_in, Tp = vap.copy(), T_cond_C + K0
    for T in (90.0, 80.0, 70.0):
        zz = s_in / s_in.sum()
        b3, x3, y3, _ = flash_TP(pk3, T, P_nru, zz)
        R["coldbox"]["hesep"].append(dict(F=float(s_in.sum()), z=zz.tolist(),
                                          T_in=float(Tp), T_out=T))
        R["coldbox"]["n2liq"].append(dict(F=float((1 - b3) * s_in.sum()), T=T))
        s_in, Tp = b3 * s_in.sum() * y3, T
    he_prod = crude[0]                              # adsorber: He only
    he_rec = he_prod / bog3[0]
    # refrigeration: condensed N2 flashed to 1.3 bar boils at ~79 K and
    # carries the 90/80 K stages; the 70 K stage needs the helium loop.
    Tn = CP.PropsSI("T", "P", 1.3e5, "Q", 0, "Nitrogen")
    Q_n2_boil = n2_liq.sum() * (CP.PropsSI("Hmolar", "P", 1.3e5, "Q", 1, "Nitrogen")
                                - CP.PropsSI("Hmolar", "P", 1.3e5, "Q", 0, "Nitrogen"))
    Q_70 = stages and Q_hesep * 0.0
    R["hesep"] = dict(stages=stages, He_recovery=float(he_rec),
                      crude_He=float(crude[0] / crude.sum()),
                      he_product_kmolh=he_prod / KMOLH,
                      he_product_kgh=he_prod * 4.0026e-3 * 3600,
                      n2_recovered_kmolh=n2_liq.sum() / KMOLH,
                      Q_cold_MW=Q_hesep / 1e6, Q_n2_boil_MW=Q_n2_boil / 1e6)

    # ---------------- HLS ---------------------------------------------
    H = collins(he_prod, x1=hls_p["x1"], x2=hls_p["x2"], T_E1=hls_p["T_E1"],
                T_E2=hls_p["T_E2"], T_D=hls_p["T_D"], T_JT=hls_p["T_JT"],
                # the nitrogen has already boiled in the helium separator,
                # so it reaches the liquefier as vapour just above saturation
                n2_cold=(n2_liq.sum(), Tn + 0.5))
    R["hls"] = dict(ok=bool(H["ok"]), W_MW=H["W"] / 1e6, dTmin=H["dTmin"],
                    kWh_per_L=H["specific_kWh_per_L"], y_jt=H["y_jt"],
                    cycle_kmolh=H["m_c"] / KMOLH,
                    pinches={k: float(v) for k, v in H["pinches"].items()},
                    x1=hls_p["x1"], x2=hls_p["x2"])

    # ---------------- NRHR auxiliary refrigeration (estimated) --------
    W_cond = refrig_power(res["Q_cond"], T_cond_C + K0)
    deficit = max(Q_hesep - Q_n2_boil, 0.0)
    W_hesep = refrig_power(deficit, 80.0)
    R["nrhr_refrig"] = dict(W_cond_MW=W_cond / 1e6, W_hesep_MW=W_hesep / 1e6,
                            method=f"{EXERGY_EFF:.0%} of Carnot",
                            estimate_MW=(W_cond + W_hesep) / 1e6)
    if nrhr == "cycle":
        from coldbox import process_curves, evaluate as cb_eval, cascade, composites
        case = "bogcold" if bog_cold else "sealed"
        proc = process_curves(R["coldbox"], bog_cold=(T_ef, P_t) if bog_cold else None)
        par = json.load(open(f"nrhr_best_{case}.json"))["params"]
        CB = cb_eval(proc, par)
        if CB is None:
            # The helium cycle was sized exactly for the design feed. For a
            # different feed, re-tune it: a slightly colder expander inlet and
            # a larger helium flow, taking the first setting that closes.
            base = dict(par)
            for dT in (0.75, 1.5, 3.0, 5.0):
                for fac in (1.0, 1.1, 1.25, 1.5, 2.0, 3.0, 4.0):
                    par = dict(base, He_T_e=base["He_T_e"] - dT, m_he=base["m_he"] * fac)
                    CB = cb_eval(proc, par)
                    if CB is not None:
                        say(f"NRHR helium cycle re-tuned for this feed: expander inlet {dT:g} K colder, flow x{fac:g}")
                        break
                if CB is not None:
                    break
        if CB is None:
            raise RuntimeError("NRHR cold box cannot be closed")
        W_cond, W_hesep = CB["W_n2"], CB["W_he"]
        from coldbox import n2_cycle, he_cycle
        n2c = n2_cycle(par["P_hi"] * 1e5, par["P_lo"] * 1e5, par["T_e1"], par["T_e2"], par["f1"])
        hec = he_cycle(par["He_P_hi"] * 1e5, par["He_P_lo"] * 1e5, par["He_T_e"])
        allc = (proc + [c.scaled(CB["m_n2"]) for c in n2c["curves"]]
                + [c.scaled(par["m_he"]) for c in hec["curves"]])
        R["nrhr_refrig"] = dict(method="simulated N2 + He reverse-Brayton cycles",
                                case=case, W_n2_MW=W_cond / 1e6, W_he_MW=W_hesep / 1e6,
                                W_exp_MW=CB["W_exp"] / 1e6,
                                m_n2_kmolh=CB["m_n2"] / KMOLH,
                                m_he_kmolh=par["m_he"] / KMOLH, params=par,
                                n2=CB["n2"], he=CB["he"], pinch_T=CB["pinch_T"],
                                duties_MW={k: v / 1e6 for k, v in CB["duties"].items()},
                                deficit_MW=cascade(proc)[1] / 1e6,
                                composites=composites(allc),
                                estimate_MW=R["nrhr_refrig"]["estimate_MW"])
        say(f"NRHR cold box: N2 cycle {W_cond/1e6:.1f} MW  He cycle {W_hesep/1e6:.1f} MW")

    # ---------------- products & totals --------------------------------
    lng_prod = lng_liq + rich_ret - tank_bog * 0 # tank BOG leaves the tank
    lng_prod = lng_liq + rich_ret - tank_bog
    mwv = PKG.MW
    lng_kgh = float(lng_prod @ mwv) * 3.6
    W_tot = (L["W"] + W_bog + H["W"] + W_cond + W_hesep)
    R["products"] = dict(
        lng_kmolh=lng_prod.sum() / KMOLH, lng_tph=lng_kgh / 1000,
        lng_comp={c: float(lng_prod[i] / lng_prod.sum()) for i, c in enumerate(C)},
        lng_N2=float(lng_prod[C.index("N2")] / lng_prod.sum()),
        he_kgh=R["hesep"]["he_product_kgh"],
        he_Lh=R["hesep"]["he_product_kgh"] / 0.125,
        n2_kmolh=R["hesep"]["n2_recovered_kmolh"],
        hhc_kmolh=R["lngp"]["hhc_kmolh"])
    R["power_MW"] = dict(c3mr=L["W"] / 1e6, bog_compressor=W_bog / 1e6,
                         nru_condenser=W_cond / 1e6, he_separation=W_hesep / 1e6,
                         helium_liquefier=H["W"] / 1e6, total=W_tot / 1e6)
    R["kpi"] = dict(kWh_per_t_LNG=W_tot / 1e3 / (lng_kgh / 1000),
                    He_overall_recovery=float(he_prod / ng[0]),
                    CH4_to_LNG=float(lng_prod[1] / ng[1]))
    R["wall_s"] = time.time() - t0
    say(f"done in {R['wall_s']:.0f}s  total power {W_tot/1e6:.1f} MW")
    return R


if __name__ == "__main__":
    import sys
    mode = sys.argv[1] if len(sys.argv) > 1 else "estimate"
    if mode == "estimate":        # the sealed, blind flowsheet
        R = run(nrhr="estimate"); out = "plant_result.json"
    elif mode == "cycle":         # same flowsheet, refrigeration simulated
        R = run(nrhr="cycle"); out = "plant_result_cycle.json"
    elif mode == "bogcold":       # plus BOG cold recovery
        R = run(nrhr="cycle", bog_cold=True); out = "plant_result_bogcold.json"
    else:                         # simulated cycles + designs from the strong optimisers
        import os
        mf = "mr_best_strong.json" if os.path.exists("mr_best_strong.json") else "mr_best.json"
        R = run(nrhr="cycle", mr_file=mf, hls_file="hls_best_strong.json")
        R["inputs"]["design_files"] = [mf, "hls_best_strong.json"]
        out = "plant_result_strong.json"
    json.dump(R, open(out, "w"), indent=1, default=float)
