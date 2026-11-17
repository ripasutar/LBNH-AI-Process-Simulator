"""
Detailed flowsheet and stream table, generated from the AI simulator.

Solves the whole plant (LNG section on my Aspen settings + tank / boil-off recycle +
N2 rejection column + helium separation + adsorber + helium liquefier) and records
every main stream: temperature, pressure, flow, vapour fraction and composition.

    python stream_report.py            -> my design (Aspen settings)
Outputs (in this folder): streams_<case>.json, used by make_pfd.py and make_excel.py
"""
import json, os, sys, copy, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); SIM = os.path.dirname(HERE); OPT = os.path.join(SIM, "section_models")
for p in (SIM, OPT):
    if p not in sys.path: sys.path.insert(0, p)
os.chdir(SIM)
import CoolProp.CoolProp as CP
from props import PKG, C, KMOLH, subpackage
from flash import flash_TP, flash_PH, bubble_T
import design_basis as DB
import econ_opt as E
import he_opt as H
from lngp import precool_loads, h_molar

K0 = 273.15; MW = PKG.MW; NRU = ["He", "CH4", "N2"]; IN = [C.index(c) for c in NRU]; PK3 = subpackage(NRU)
def to9(n3):
    o = np.zeros(len(C)); o[IN] = n3; return o

class Book:
    def __init__(self): self.s, self.u = {}, {}
    def add(self, sid, name, frm, to, n, T, P, vf=None, phase=None):
        n = np.asarray(n, float)
        if len(n) == 3: n = to9(n)
        if vf is None and n.sum() > 0:
            try:
                nz = [i for i in range(len(C)) if n[i] > 1e-12 * n.sum()]
                pk = subpackage([C[i] for i in nz]) if len(nz) < len(C) else PKG
                vf = float(flash_TP(pk, T, P, n[nz] / n[nz].sum())[0]) if len(nz) > 1 else None
            except Exception:
                vf = None
        self.s[sid] = dict(id=sid, name=name, frm=frm, to=to, T_C=float(T - K0), P_bar=float(P / 1e5),
                           kmolh=float(n.sum() / KMOLH), tph=float(n @ MW * 3.6 / 1000), vf=vf, phase=phase,
                           x={c: float(n[i] / n.sum()) if n.sum() > 0 else 0.0 for i, c in enumerate(C)})
    def unit(self, uid, name, sect, **kw): self.u[uid] = dict(id=uid, name=name, section=sect, **kw)

def c3_flows(loads, lev, Tcond):
    """Propane circulation at each level (kg/s), same logic as lngp.propane_cycle."""
    Tev = [t + K0 for t in lev]; T_cond = Tcond + K0
    hl = lambda T: CP.PropsSI("H", "T", T, "Q", 0, "Propane"); hv = lambda T: CP.PropsSI("H", "T", T, "Q", 1, "Propane")
    m_evap = [loads[k] / (hv(Tev[k]) - hl(Tev[k])) for k in range(3)]; m_down = 0.0; vap = [0] * 3; liq_in = [0] * 3
    for k in (2, 1, 0):
        h_in = hl(Tev[k - 1]) if k > 0 else hl(T_cond)
        x = (h_in - hl(Tev[k])) / (hv(Tev[k]) - hl(Tev[k]))
        m_in = (m_evap[k] + m_down) / (1 - x); vap[k] = m_evap[k] + x * m_in; liq_in[k] = m_in; m_down = m_in
    P = lambda T: CP.PropsSI("P", "T", T, "Q", 0, "Propane")
    return dict(vap=vap, liq_in=liq_in, total=liq_in[0], P=[P(t) for t in Tev], Pc=P(T_cond))

def solve(case, x_lng, colv, cycv, hlsv):
    t0 = time.time(); B = Book()
    E.set_feed(600.0)
    cap = {}
    orig = E.mche
    def wrap(*a, **k):
        r = orig(*a, **k); cap["r"] = r; cap["kw"] = k; return r
    E.mche = wrap
    o = E.evaluate(x_lng, ua_base=12.18771779439311e6)   # main-exchanger UA of my design, W/K
    E.mche = orig
    r = cap["r"]; yN2, yC1, yC2, P_lo, P_hi, dTw, T_pre = x_lng
    mr = np.zeros(len(C)); mr[C.index("N2")], mr[C.index("CH4")], mr[C.index("C2H6")] = yN2, yC1, yC2; mr[C.index("C3H8")] = 1 - yN2 - yC1 - yC2
    F_mr = r["F_mr"]; P_ng = 49.6e5
    # ---------------- LNG section ----------------
    B.add("1", "Natural gas feed", "Battery limit", "Pre-treatment", E.z_feed * E.F_feed, 25 + K0, 50e5)
    B.add("2", "Treated gas to main exchanger", "Propane pre-coolers", "MCHE", E.z_ng * E.F_ng, T_pre + K0, P_ng)
    B.add("3", "Heavy hydrocarbons (C2+ removal, Aspen split)", "Heavy-HC separator", "Battery limit", E.z_feed * E.F_feed - E.z_ng * E.F_ng, -2.2 + K0, 49.7e5)
    B.add("4", "LNG leaving main exchanger", "MCHE", "End-flash valve", r["lng"], -152 + K0, P_ng, vf=0.0)
    B.add("5", "MR low-pressure vapour", "MCHE (warm end)", "MR compressor", mr * F_mr, r["T_warm"], P_lo * 1e5)
    B.add("6", "MR compressed and cooled", "MR compressor + after-cooler", "Propane pre-coolers", mr * F_mr, 40 + K0, P_hi * 1e5)
    B.add("7", "MR pre-cooled", "Propane pre-coolers", "MCHE (hot pass)", mr * F_mr, T_pre + K0, P_hi * 1e5 - 0.31e5)
    B.add("8", "MR subcooled", "MCHE (cold end)", "JT valve", mr * F_mr, -152 + K0, P_hi * 1e5 - 0.31e5)
    B.add("9", "MR after JT valve (cold refrigerant)", "JT valve", "MCHE (cold pass)", mr * F_mr, r["T_jt"], P_lo * 1e5)
    lev = (E.LEV_ASPEN[0], E.LEV_ASPEN[1], T_pre - E.COLD_APPROACH); steps = (lev[0] + 4, lev[1] + 4, T_pre)
    loads = np.array(precool_loads(E.z_feed, E.F_feed, 50e5, 25 + K0, T_steps=steps)) + np.array(precool_loads(mr, F_mr, P_hi * 1e5, 40 + K0, T_steps=steps))
    c3 = c3_flows(list(loads), lev, E.TCOND)
    c3n = lambda kg: np.eye(len(C))[C.index("C3H8")] * kg * 1000 / MW[C.index("C3H8")]
    B.add("C1", "Propane liquid from condenser", "Propane condenser", "Level-1 evaporator", c3n(c3["total"]), E.TCOND + K0, c3["Pc"], vf=0.0)
    for k in range(3):
        B.add(f"C{k + 2}", f"Propane vapour, level {k + 1}", f"Level-{k + 1} evaporator", "Propane compressor", c3n(c3["vap"][k]), lev[k] + K0, c3["P"][k], vf=1.0)
    B.unit("MCHE", "Main cryogenic heat exchanger", "LNG", duty_MW=r["Q"] / 1e6, dTmin_K=o["dTmin"], UA_MWK=o["UA_MWK"])
    B.unit("K-MR", "MR compressor (3 stages)", "LNG", power_MW=o["W_mr"], P_in_bar=P_lo, P_out_bar=P_hi)
    B.unit("K-C3", "Propane compressor (3 levels)", "LNG", power_MW=o["W_c3"], levels_C=[round(t, 1) for t in lev], cond_C=round(E.TCOND, 1))
    B.unit("E-C3", "Propane pre-coolers (3 levels)", "LNG", duty_MW=float(loads.sum() / 1e6), duties_MW=[float(q / 1e6) for q in loads])
    # ---------------- tank, end flash, boil-off recycle ----------------
    lng_hp = r["lng"]; z = lng_hp / lng_hp.sum(); P_t = DB.END_FLASH_P_BAR * 1e5
    T_ef, b, x, y = flash_PH(PKG, P_t, z, h_molar(-152 + K0, P_ng, z), lo=80, hi=200)
    fl_v, fl_l = b * lng_hp.sum() * y, (1 - b) * lng_hp.sum() * x
    P_bog = colv[3] * 1e5
    ret, rfv = np.zeros(len(C)), np.zeros(len(C))
    H.CB0 = copy.deepcopy(H.CB0)
    for it in range(12):
        tank = fl_l + ret
        Tb, xb, yb = bubble_T(PKG, P_t, tank / tank.sum(), lo=90, hi=200)
        tbog = yb * tank.sum() * DB.TANK_BOG_PER_HOUR
        bog = fl_v + tbog + rfv; bog3 = bog[IN].copy(); bog3[1] += bog.sum() - bog3.sum()
        H.BOG = bog3; H.CB0["bog"] = dict(H.CB0["bog"], F=float(bog3.sum()), z=(bog3 / bog3.sum()).tolist())
        cs = H.column(*colv)
        if cs is None: raise RuntimeError("column did not converge")
        bot = np.array(cs["cb"]["rich"]["z"]) * cs["cb"]["rich"]["F"]; zb = bot / bot.sum()
        Tr, br, xr, yr = flash_PH(PKG, P_t, to9(zb), h_molar(cs["cb"]["rich"]["T_out"], P_bog, to9(zb)), lo=80, hi=200)
        new_ret, new_rfv = to9(xr[IN] * (1 - br) * bot.sum()), to9(yr[IN] * br * bot.sum())
        err = (abs(new_ret - ret).sum() + abs(new_rfv - rfv).sum()) / bog.sum(); ret, rfv = new_ret, new_rfv
        if err < 1e-5: break
    cb = H.coldbox(cs, cycv); hl = H.liquefier(cs, hlsv)
    if cb is None or hl is None: raise RuntimeError("cold box or liquefier did not close")
    lng_prod = fl_l + ret - tbog
    B.add("10", "End-flash vapour", "End-flash valve", "BOG header", fl_v, T_ef, P_t, vf=1.0)
    B.add("11", "LNG to storage tank", "End-flash valve", "LNG tank", fl_l, T_ef, P_t, vf=0.0)
    B.add("12", "Tank boil-off (heat leak)", "LNG tank", "BOG header", tbog, Tb, P_t, vf=1.0)
    B.add("13", "Rich-LNG return, flash vapour", "Rich-LNG valve", "BOG header", rfv, Tr, P_t, vf=1.0)
    B.add("14", "Rich-LNG return, liquid", "Rich-LNG valve", "LNG tank", ret, Tr, P_t, vf=0.0)
    B.add("15", "LNG product to ship", "LNG tank", "Battery limit", lng_prod, Tb, P_t, vf=0.0)
    B.add("16", "Boil-off gas to compressor", "BOG header (warmed)", "BOG compressor", bog3, 303.15, P_t, vf=1.0)
    B.add("17", "Compressed BOG", "BOG compressor + cooler", "NRHR cold box", bog3, DB.COOLER_T_C + K0, P_bog, vf=1.0)
    T_feed = cs["cb"]["bog"]["T_out"]
    B.add("18", "Column feed (cooled BOG)", "NRHR cold box", "N2 rejection column", bog3, T_feed, P_bog)
    cond = cs["cb"]["cond"]
    B.add("19", "Column overhead vapour", "Column top", "Overhead condenser", np.array(cond["z"]) * cond["F"], cond["T_in"], P_bog, vf=1.0)
    B.add("20", "Column bottoms (rich LNG)", "Column bottom", "NRHR cold box", bot, cs["cb"]["rich"]["T_in"], P_bog, vf=0.0)
    B.add("21", "Rich LNG subcooled", "NRHR cold box", "Rich-LNG valve", bot, cs["cb"]["rich"]["T_out"], P_bog, vf=0.0)
    # helium separators
    s_in = None; n2liq = np.zeros(3); Tp = colv[0] + K0; temps = (90.0, 80.0, float(colv[2]))
    vap0 = np.array(cs["cb"]["hesep"][0]["z"]) * cs["cb"]["hesep"][0]["F"]
    B.add("22", "Condenser vapour (He + N2)", "Overhead condenser", "He separator 1", vap0, Tp, P_bog, vf=1.0)
    s_in = vap0
    for k, T in enumerate(temps):
        zz = s_in / s_in.sum(); b3, x3, y3, _ = flash_TP(PK3, T, P_bog, zz)
        liq = (1 - b3) * s_in.sum() * x3; vap = b3 * s_in.sum() * y3; n2liq += liq
        B.add(f"{23 + 2 * k}", f"He separator {k + 1} vapour", f"He separator {k + 1}", "He separator " + str(k + 2) if k < 2 else "Adsorber", vap, T, P_bog, vf=1.0)
        B.add(f"{24 + 2 * k}", f"He separator {k + 1} liquid (N2-rich)", f"He separator {k + 1}", "N2 flash (1.3 bar)", liq, T, P_bog, vf=0.0)
        s_in = vap
    crude = s_in; he = np.array([crude[0], 0, 0])
    Tn = CP.PropsSI("T", "P", 1.3e5, "Q", 0, "Nitrogen")
    B.add("29", "Helium product gas, ≥ 99.999 %", "Adsorber", "Helium liquefier", he, temps[2], P_bog, vf=1.0)
    B.add("30", "Adsorber tail gas", "Adsorber", "BOG header / fuel", crude - he, temps[2], 1.3e5, vf=1.0)
    B.add("31", "Liquid helium product", "Helium liquefier", "Battery limit", he, 4.22, 1.013e5, vf=0.0)
    B.add("32", "Nitrogen product", "N2 flash + cold box", "Battery limit", n2liq, Tn + 0.5, 1.3e5, vf=1.0)
    p = cb["p"]
    B.unit("K-BOG", "BOG compressor (3 stages)", "BOG", power_MW=cs["W_bog"], P_in_bar=P_t / 1e5, P_out_bar=colv[3])
    B.unit("T-NRU", "N2 rejection column (20 stages)", "NRHR", condenser_C=colv[0], boilup=colv[1], Q_cond_MW=cs["Q_cond_MW"], Q_reb_MW=cs["Q_reb_MW"])
    B.unit("K-N2", "N2 refrigeration cycle (reverse Brayton, 2 expanders)", "NRHR", power_MW=cb["W_n2"], P_hi_bar=p["P_hi"], P_lo_bar=p["P_lo"], T_exp_K=[p["T_e1"], p["T_e2"]], flow_kmolh=cb["m_n2_tph"] * 1000 / 28.0134)
    B.unit("K-He", "He refrigeration cycle", "NRHR", power_MW=cb["W_he"], P_hi_bar=p["He_P_hi"], P_lo_bar=p["He_P_lo"], T_exp_K=p["He_T_e"])
    B.unit("HLS", "Helium liquefier (Collins cycle)", "HLS", power_MW=hl["W"], dTmin_K=hl["dTmin"], **{k: float(v) for k, v in hl["q"].items()})
    W = dict(MR=o["W_mr"], C3=o["W_c3"], BOG=cs["W_bog"], N2=cb["W_n2"], He=cb["W_he"], HLS=hl["W"])
    out = dict(case=case, streams=B.s, units=B.u, power_MW=W, total_MW=sum(W.values()), he_purity_crude=float(crude[0] / crude.sum()),
               he_recovery=float(crude[0] / bog3[0]), ch4_recovery=float(bot[1] / bog3[1]), lng_tph=B.s["15"]["tph"], he_kgh=B.s["31"]["tph"] * 1000,
               seconds=time.time() - t0)
    json.dump(out, open(os.path.join(HERE, f"streams_{case}.json"), "w"), indent=1, default=float)
    print(case, f"total {out['total_MW']:.1f} MW", {k: round(v, 1) for k, v in W.items()}, f"LNG {out['lng_tph']:.1f} t/h He {out['he_kgh']:.0f} kg/h", f"{out['seconds']:.0f} s")
    return out

if __name__ == "__main__":
    solve("my_design", E.X_ASPEN, H.COL_START, H.CYC_START, H.HLS_START)
