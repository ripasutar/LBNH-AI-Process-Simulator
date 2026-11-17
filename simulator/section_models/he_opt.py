"""
Economic optimization of the helium sections (AI simulator's own layout):
N2 rejection column + helium separators, NRHR cold box (N2 and He reverse-Brayton
cycles) and the Collins helium liquefier.

Maximize the helium sections' contribution to plant profit:
    profit = helium revenue + value of methane returned to LNG
             - electricity - CRF x (compressor capital + cold-box exchanger capital)

Decision variables (18)
  column:    condenser temperature (C), boilup ratio, last helium-separator temperature (K)
  cold box:  N2 cycle high / low pressure (bar), expander inlet temperatures T1, T2 (K),
             split to expander 1; He cycle high / low pressure (bar), expander inlet (K),
             helium cycle flow (mol/s)
  liquefier: expander flow fractions x1, x2, expander inlets T_E1, T_E2 (K),
             intermediate temperature T_D (K), Joule-Thomson inlet T_JT (K)
Constraints
  helium purity: crude helium >= 98 % He before the adsorber, which polishes it
  to >= 99.999 %;  helium recovery >= 98 %;  methane recovery to LNG >= 98 %;
  cold box: no refrigeration left unmet at a 2 K minimum approach;
  liquefier: every exchanger >= 1 K approach along its whole length;
  pressure ratios and bounds below; the column and the recycle must converge.
"""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); SIM = os.path.dirname(HERE)
if SIM not in sys.path: sys.path.insert(0, SIM)
import CoolProp.CoolProp as CP
from props import subpackage
from flash import flash_TP, mixture_H
from nru import nru_with_condenser
from coldbox import process_curves, evaluate as cb_evaluate, n2_cycle, he_cycle, composites
from hls import collins
from flash import dew_T
from lngp import compress
import design_basis as DB
from plant import to9
K0 = 273.15; KMOLH = 1000 / 3600
HOURS, ELEC, CRF = 8322.0, 0.0775, 0.1 / (1 - 1.1 ** -20)
COMP_COST = 1000.0                   # $/kW installed compressor + motor (same as the LNG section)
HX_COST = 100e6 / 12.188e6           # $ per W/K of UA: same unit cost as the main LNG exchanger
P_HE, P_LNG = 100.0, 330.0           # $/kg helium, $/t LNG (paper)
MW_HE, MW_CH4 = 4.0026, 16.043
PK3 = subpackage(["He", "CH4", "N2"])

R0 = json.load(open(os.path.join(HERE, "he_design_state.json")))
CB0 = R0["coldbox"]; P_NRU = CB0["P"]
BOG = np.array(CB0["bog"]["z"]) * CB0["bog"]["F"]; T_FEED = CB0["bog"]["T_out"]
VF_FEED = flash_TP(PK3, T_FEED, P_NRU, BOG / BOG.sum())[0]

CYC_KEYS = ["P_hi", "P_lo", "T_e1", "T_e2", "f1", "He_P_hi", "He_P_lo", "He_T_e", "m_he"]
CYC_LO = np.array([20, 1.5, 150, 95, 0.1, 10, 2, 70, 0]); CYC_HI = np.array([80, 12, 300, 200, 0.9, 40, 15, 130, 1500])
HLS_KEYS = ["x1", "x2", "T_E1", "T_E2", "T_D", "T_JT"]
HLS_LO = np.array([0.05, 0.05, 25, 9, 13, 5.4]); HLS_HI = np.array([0.7, 0.6, 75, 22, 71, 8.5])
COL_KEYS = ["T_cond_C", "boilup", "T_sep3", "P_bog_bar"]
COL_LO = np.array([-162.0, 0.5, 66.0, 20.0]); COL_HI = np.array([-150.0, 1.1, 76.0, 28.0])

# starting point: the AI simulator's helium design (refrigeration sized for minimum power only)
_cb0 = json.load(open(os.path.join(SIM, "nrhr_best_sealed.json")))["params"]; _hl0 = json.load(open(os.path.join(SIM, "hls_best.json")))
CYC_START = [_cb0[k] for k in CYC_KEYS]
HLS_START = [_hl0[k] for k in HLS_KEYS]
COL_START = [-156.0, 0.7424, 70.0, DB.BOG_COMPRESSOR_P_BAR]

_warm = {}
def column(T_cond_C, boilup, T_sep3, P_bog_bar=25.0):
    """BOG compressor + N2 rejection column + 3 helium separators -> cold-box duties and products."""
    P = P_bog_bar * 1e5; z3 = BOG / BOG.sum()
    try:
        W_bog, _, _ = compress(to9(z3), BOG.sum(), 303.15, DB.END_FLASH_P_BAR * 1e5, P, 3)
        Td, *_ = dew_T(PK3, P, z3)
    except Exception:
        return None
    T_feed = Td - 8.0; vf = flash_TP(PK3, T_feed, P, z3)[0]
    st = _warm.get("st") if abs(_warm.get("P", P) - P) < 0.5e5 else None
    try:
        res = nru_with_condenser(BOG, P, T_feed, vf, T_cond_C + K0, boilup, state=st)
        if res is None and st is not None:
            res = nru_with_condenser(BOG, P, T_feed, vf, T_cond_C + K0, boilup)
    except Exception:
        res = None
    if res is None:
        return None
    _warm["st"] = dict(reflux=res["reflux"], x=res["x"]); _warm["P"] = P
    P_NRU = P
    col = res["col"].results(); bottoms = np.array(col["l"][-1]); vap = np.array(res["vapour"])
    top_v = np.array(col["v"][0]); zb = bottoms / bottoms.sum()
    cb = dict(P=P_NRU, bog=dict(CB0["bog"], T_out=float(T_feed)),
              cond=dict(F=float(top_v.sum()), z=(top_v / top_v.sum()).tolist(), T_in=float(col["T"][0]), T_out=T_cond_C + K0, Q=float(abs(res["Q_cond"]))),
              reb=dict(T=float(col["T"][-1]), T_above=float(col["T"][-2]), Q=float(res["reb"])),
              rich=dict(F=float(bottoms.sum()), z=zb.tolist(), T_in=float(col["T"][-1]), T_out=CB0["rich"]["T_out"]),
              hesep=[], n2liq=[])
    s_in, Tp, n2_liq = vap.copy(), T_cond_C + K0, 0.0
    for T in (90.0, 80.0, float(T_sep3)):
        zz = s_in / s_in.sum(); b3, x3, y3, _ = flash_TP(PK3, T, P_NRU, zz)
        cb["hesep"].append(dict(F=float(s_in.sum()), z=zz.tolist(), T_in=float(Tp), T_out=T))
        cb["n2liq"].append(dict(F=float((1 - b3) * s_in.sum()), T=T))
        n2_liq += (1 - b3) * s_in.sum(); s_in, Tp = b3 * s_in.sum() * y3, T
    crude = s_in
    return dict(cb=cb, he_prod=float(crude[0]), purity=float(crude[0] / crude.sum()), he_rec=float(crude[0] / BOG[0]),
                ch4_rec=float(bottoms[1] / BOG[1]), ch4_lng=float(bottoms[1]), n2_liq=float(n2_liq), W_bog=W_bog / 1e6,
                Q_cond_MW=abs(res["Q_cond"]) / 1e6, Q_reb_MW=res["reb"] / 1e6)

def ua_of(allc):
    c = composites(allc)
    Th, Hh = np.array(c["hot"]["T"]), np.array(c["hot"]["H"]) * 1e6
    Tc, Hc = np.array(c["cold"]["T"]), np.array(c["cold"]["H"]) * 1e6
    Q = np.linspace(0, min(Hh[-1], Hc[-1]), 400)
    dT = np.maximum(np.interp(Q, Hh, Th) - np.interp(Q, Hc, Tc), 0.5)
    return float(np.sum(np.diff(Q) / (0.5 * (dT[1:] + dT[:-1]))))

def coldbox(cs, v):
    p = dict(zip(CYC_KEYS, map(float, np.clip(v, CYC_LO, CYC_HI))))
    if p["He_P_hi"] < 1.3 * p["He_P_lo"] or p["P_hi"] < 1.5 * p["P_lo"]:
        return None
    try:
        proc = process_curves(cs["cb"]); r = cb_evaluate(proc, p)
        if r is None: return None
        n2 = n2_cycle(p["P_hi"] * 1e5, p["P_lo"] * 1e5, p["T_e1"], p["T_e2"], p["f1"]); he = he_cycle(p["He_P_hi"] * 1e5, p["He_P_lo"] * 1e5, p["He_T_e"])
        allc = proc + [c.scaled(r["m_n2"]) for c in n2["curves"]] + [c.scaled(p["m_he"]) for c in he["curves"]]
        return dict(W=r["W"] / 1e6, W_n2=r["W_n2"] / 1e6, W_he=r["W_he"] / 1e6, UA=ua_of(allc) / 1e6, m_n2_tph=r["m_n2"] * 28.0134 * 3.6 / 1000, p=p)
    except Exception:
        return None

def liquefier(cs, v):
    q = dict(zip(HLS_KEYS, map(float, np.clip(v, HLS_LO, HLS_HI))))
    if q["x1"] + q["x2"] >= 0.92: return None
    q["T_D"] = float(np.clip(q["T_D"], q["T_E2"] + 4, q["T_E1"] - 4))
    Tn = CP.PropsSI("T", "P", 1.3e5, "Q", 0, "Nitrogen")
    try:
        r = collins(cs["he_prod"], x1=q["x1"], x2=q["x2"], T_E1=q["T_E1"], T_E2=q["T_E2"], T_D=q["T_D"], T_JT=q["T_JT"], n2_cold=(cs["n2_liq"], Tn + 0.5))
    except Exception:
        return None
    if "dTmin" not in r: return None
    return dict(W=r["W"] / 1e6, dTmin=float(r["dTmin"]), q=q)

def money(cs, cbr, hlr):
    """$M/yr: the parts of plant profit the helium sections control."""
    rev_he = cs["he_prod"] * MW_HE * 3.6 * HOURS * P_HE / 1e6
    val_ch4 = cs["ch4_lng"] * MW_CH4 * 3.6 / 1000 * HOURS * P_LNG / 1e6
    W = cbr["W"] + hlr["W"] + cs["W_bog"]
    elec = W * 1e3 * HOURS * ELEC / 1e6
    capex = (COMP_COST * W * 1e3 + HX_COST * cbr["UA"] * 1e6) / 1e6
    return dict(rev_he=rev_he, val_ch4=val_ch4, elec=elec, capex=capex, W=W, profit=rev_he + val_ch4 - elec - CRF * capex)

def penalties(cs, hlr):
    pen = 0.0
    if cs["purity"] < 0.98: pen += 2000 * (0.98 - cs["purity"]) + 5
    if cs["he_rec"] < 0.98: pen += 2000 * (0.98 - cs["he_rec"]) + 5
    if cs["ch4_rec"] < 0.98: pen += 2000 * (0.98 - cs["ch4_rec"]) + 5
    if hlr is not None and hlr["dTmin"] < 1.0: pen += 30 * (1.0 - hlr["dTmin"]) + 5
    return pen

def nm(f, x0, lo, hi, maxfev, deadline=None):
    from scipy.optimize import minimize
    sc = hi - lo; u0 = (np.array(x0, float) - lo) / sc
    best = [None, 1e18]
    def g(u):
        if deadline and time.time() > deadline: raise StopIteration
        x = lo + np.clip(u, 0, 1) * sc; v = f(x)
        if v < best[1]: best[0], best[1] = x.copy(), v
        return v
    try:
        minimize(g, u0, method="Nelder-Mead", options=dict(maxfev=maxfev, xatol=1e-4, fatol=1e-3, adaptive=True))
    except StopIteration:
        pass
    return best[0], best[1]

def evaluate_all(colv, cycv, hlsv):
    cs = column(*colv)
    if cs is None: return None
    cbr = coldbox(cs, cycv); hlr = liquefier(cs, hlsv)
    if cbr is None or hlr is None: return dict(cs=cs, ok=False)
    m = money(cs, cbr, hlr)
    return dict(cs=cs, cb=cbr, hl=hlr, m=m, pen=penalties(cs, hlr), ok=True)

def optimize(log, colv=None, cycv=None, hlsv=None, budget_cyc=250, budget_hls=300, budget_col=0, inner_cyc=60, inner_hls=80, deadline=None, on_step=None, time_cyc=None, time_hls=None):
    """Run the helium-section optimization; every simulator solve is appended to log."""
    colv = list(colv or COL_START); cycv = list(cycv or CYC_START); hlsv = list(hlsv or HLS_START)
    base = evaluate_all(colv, cycv, hlsv)
    if base is None or not base["ok"]: raise RuntimeError("starting design could not be solved")
    P0 = base["m"]["profit"]
    def rec(section, cs, cbr, hlr, feasible):
        if cbr is None or hlr is None: log.append(dict(section=section, gain=None, feasible=False)); return
        m = money(cs, cbr, hlr); log.append(dict(section=section, gain=round(m["profit"] - P0, 3), feasible=bool(feasible), W=round(m["W"], 3)))
        if on_step: on_step()
    log.append(dict(section="start", gain=0.0, feasible=True, W=round(base["m"]["W"], 3)))
    cs = base["cs"]; hlr_cur = base["hl"]; cbr_cur = base["cb"]
    # 1) cold box refrigeration (N2 + He cycles), column fixed
    def f_cyc(v):
        cbr = coldbox(cs, v)
        if cbr is None: rec("cold box", cs, None, None, False); return 1e4
        m = money(cs, cbr, hlr_cur); rec("cold box", cs, cbr, hlr_cur, True); return -m["profit"]
    x, _ = nm(f_cyc, cycv, CYC_LO, CYC_HI, budget_cyc, time.time() + time_cyc if time_cyc else deadline); cycv = list(x); cbr_cur = coldbox(cs, cycv)
    # 2) helium liquefier
    def f_hls(v):
        hlr = liquefier(cs, v)
        if hlr is None: rec("liquefier", cs, None, None, False); return 1e4
        pen = penalties(cs, hlr); m = money(cs, cbr_cur, hlr); rec("liquefier", cs, cbr_cur, hlr, pen == 0); return -m["profit"] + pen
    x, _ = nm(f_hls, hlsv, HLS_LO, HLS_HI, budget_hls, time.time() + time_hls if time_hls else deadline); hlsv = list(x); hlr_cur = liquefier(cs, hlsv)
    # 3) N2 rejection column and separators, with the cycles and liquefier re-tuned each time
    if budget_col:
        state = dict(cyc=cycv, hls=hlsv)
        def f_col(v):
            c2 = column(*v)
            if c2 is None: log.append(dict(section="column", gain=None, feasible=False)); return 1e4
            xc, _ = nm(lambda w: (lambda r: 1e4 if r is None else -money(c2, r, liquefier(c2, state["hls"]) or dict(W=99, dTmin=0))["profit"])(coldbox(c2, w)), state["cyc"], CYC_LO, CYC_HI, inner_cyc, deadline)
            cb2 = coldbox(c2, xc)
            xh, _ = nm(lambda w: (lambda r: 1e4 if r is None else -money(c2, cb2, r)["profit"] + penalties(c2, r))(liquefier(c2, w)) if cb2 else 1e4, state["hls"], HLS_LO, HLS_HI, inner_hls, deadline)
            hl2 = liquefier(c2, xh)
            if cb2 is None or hl2 is None: log.append(dict(section="column", gain=None, feasible=False)); return 1e4
            pen = penalties(c2, hl2); m = money(c2, cb2, hl2)
            log.append(dict(section="column", gain=round(m["profit"] - P0, 3), feasible=pen == 0, W=round(m["W"], 3)))
            if pen == 0 and -m["profit"] < state.get("best", 1e18):
                state.update(best=-m["profit"], cyc=list(xc), hls=list(xh), col=list(v))
            if on_step: on_step()
            return -m["profit"] + pen
        nm(f_col, colv, COL_LO, COL_HI, budget_col, deadline)
        if "col" in state: colv, cycv, hlsv = state["col"], state["cyc"], state["hls"]
    best = evaluate_all(colv, cycv, hlsv)
    return dict(base=summary(base, colv=list(COL_START)), best=summary(best, colv=colv), colv=colv, cycv=cycv, hlsv=hlsv, P0=P0)

def summary(e, colv=None):
    cs, cb, hl, m = e["cs"], e["cb"], e["hl"], e["m"]
    return dict(profit=m["profit"], W=m["W"], W_bog=cs["W_bog"], W_cb=cb["W"], W_n2=cb["W_n2"], W_he=cb["W_he"], W_hls=hl["W"], UA=cb["UA"], m_n2_tph=cb["m_n2_tph"],
                elec=m["elec"], capex=m["capex"], rev_he=m["rev_he"], val_ch4=m["val_ch4"], he_kgh=cs["he_prod"] * MW_HE * 3.6,
                col=dict(T_cond_C=colv[0], boilup=colv[1], T_sep3=colv[2], P_bog_bar=colv[3]) if colv else None, purity=cs["purity"], he_rec=cs["he_rec"], ch4_rec=cs["ch4_rec"], hls_dTmin=hl["dTmin"],
                cyc=cb["p"], hls=hl["q"])
