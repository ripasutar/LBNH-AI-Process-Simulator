"""
NRHR cold box with simulated refrigeration
==========================================

Replaces the 30 %-of-Carnot estimate. Every stream that must be cooled or
heated inside the nitrogen-rejection cold box goes into one heat-exchange
network, together with the two refrigeration cycles the design uses
("N2 and He refrigeration"):

  process, to be cooled        BOG feed 313 K -> column feed temperature
                               column overhead condenser
                               rich LNG bottoms subcooled to -158 C
                               helium separators 117 -> 90 -> 80 -> 70 K
  process, to be heated        column reboiler
                               recovered liquid N2, throttled to 1.3 bar
                               and boiled (its warming above 80 K is
                               used in the helium liquefier, not here)
  N2 reverse-Brayton cycle     two turbo-expanders, warm and cold
  He reverse-Brayton cycle     one turbo-expander, for the 70-80 K end,
                               below where nitrogen can stay a gas

Feasibility is checked with the problem-table cascade at a 2 K minimum
approach: no refrigeration may be needed beyond what the two cycles
supply. For given cycle pressures and expander temperatures the smallest
nitrogen flow that closes the box is found by bisection; the pressures,
temperatures, split and helium flow are optimised for minimum compressor
power. Expander work is not recovered, the same convention as the helium
liquefier.
"""
import json

import numpy as np
import CoolProp.CoolProp as CP

import design_basis as DB
from props import subpackage
from flash import flash_TP, mixture_H

DT_MIN = 2.0
T_AMB = DB.COOLER_T_C + 273.15          # 313.15 K
T_WARM_OUT = T_AMB - DT_MIN             # cycle gas leaves the box here
NRU = ["He", "CH4", "N2"]
P_N2_BOIL = 1.3e5


# ------------------------------------------------------------------ streams
class Curve:
    """Heat content H(T) of one stream, W, over [T_lo, T_hi]; hot=True if
    it gives heat up as it cools."""
    def __init__(self, name, T, H, hot):
        o = np.argsort(T)
        self.name, self.T, self.H, self.hot = name, np.asarray(T)[o], np.asarray(H)[o], hot
        self.H = self.H - self.H[0]

    def q(self, a, b):
        """Heat exchanged between temperatures a < b (clipped to range)."""
        a, b = np.clip([a, b], self.T[0], self.T[-1])
        return np.interp(b, self.T, self.H) - np.interp(a, self.T, self.H)

    def scaled(self, k):
        c = Curve(self.name, self.T.copy(), self.H * k, self.hot)
        return c


def mixture_curve(pk, name, F, z, P, T_a, T_b, hot, n=24, Q=None):
    T = np.linspace(min(T_a, T_b), max(T_a, T_b), n)
    H = []
    for t in T:
        b, x, y, _ = flash_TP(pk, t, P, np.asarray(z))
        H.append(F * mixture_H(pk, t, P, b, x, y))
    H = np.array(H)
    if Q is not None:                    # match the column's own duty
        H = (H - H[0]) * Q / (H[-1] - H[0])
    return Curve(name, T, H, hot)


def pure_curve(fluid, name, F, P, T_a, T_b, hot, n=30, h_start=None):
    """h_start: enthalpy at the cold end when that point is two-phase
    (a wet expander exhaust sitting on the saturation line)."""
    lo, hi = min(T_a, T_b), max(T_a, T_b)
    T = np.linspace(lo, hi, n)
    if h_start is not None:
        T[0] = lo
        T = np.concatenate([[lo], np.linspace(lo + 0.05, hi, n - 1)])
    H = np.array([CP.PropsSI("Hmolar", "T", t, "P", P, fluid) for t in T])
    if h_start is not None:
        H[0] = h_start
    H = H * F
    return Curve(name, T, H, hot)


def process_curves(cb, bog_cold=None):
    """bog_cold: (T_in K, P Pa) to recover the cold of the boil-off gas
    leaving the tank - warming it against the compressed BOG before the
    compressor instead of compressing it cold-wasted from 303 K."""
    pk = subpackage(NRU)
    P = cb["P"]
    out = []
    if bog_cold:
        b = cb["bog"]
        out.append(mixture_curve(pk, "BOG from tank (cold recovery)", b["F"], b["z"],
                                 bog_cold[1], bog_cold[0], b["T_in"] - 10.0, False))
    b = cb["bog"]
    out.append(mixture_curve(pk, "BOG feed", b["F"], b["z"], P, b["T_in"], b["T_out"], True))
    c = cb["cond"]
    out.append(mixture_curve(pk, "NRU condenser", c["F"], c["z"], P, c["T_in"], c["T_out"],
                             True, Q=c["Q"]))
    r = cb["rich"]
    out.append(mixture_curve(pk, "rich LNG subcooler", r["F"], r["z"], P, r["T_in"],
                             r["T_out"], True, n=12))
    for i, s in enumerate(cb["hesep"]):
        out.append(mixture_curve(pk, f"He separator {i+1}", s["F"], s["z"], P, s["T_in"],
                                 s["T_out"], True, n=16))
    rb = cb["reb"]
    out.append(Curve("NRU reboiler", [rb["T_above"], rb["T"]], [0.0, rb["Q"]], False))
    # recovered N2: liquid at 25 bar, throttled to 1.3 bar, boiled to dew
    Tb = CP.PropsSI("T", "P", P_N2_BOIL, "Q", 1, "Nitrogen")
    hV = CP.PropsSI("Hmolar", "P", P_N2_BOIL, "Q", 1, "Nitrogen")
    Q = sum(s["F"] * (hV - CP.PropsSI("Hmolar", "T", s["T"], "P", P, "Nitrogen"))
            for s in cb["n2liq"])
    out.append(Curve("liquid N2 boiling", [Tb, Tb + 0.5], [0.0, Q], False))
    return out


# ------------------------------------------------------------------ cascade
def cascade(curves, dt=DT_MIN):
    """Problem-table cascade. Returns (hot utility, cold utility, pinch T)."""
    edges = set()
    for c in curves:
        sh = -dt / 2 if c.hot else dt / 2
        edges.update((c.T + sh).tolist())
    e = np.array(sorted(edges))[::-1]
    R, run = [0.0], 0.0
    for hi, lo in zip(e[:-1], e[1:]):
        s = 0.0
        for c in curves:
            sh = -dt / 2 if c.hot else dt / 2
            q = c.q(lo - sh, hi - sh)
            s += q if c.hot else -q
        run += s
        R.append(run)
    R = np.array(R)
    Qh = max(0.0, -R.min())
    Qc = Qh + R[-1]
    return Qh, Qc, float(e[int(np.argmin(R))])


# ------------------------------------------------------------------ cycles
def expand(fluid, T_in, P_in, P_out, eff=None):
    eff = eff or DB.EXPANDER_EFF
    h1 = CP.PropsSI("Hmolar", "T", T_in, "P", P_in, fluid)
    s1 = CP.PropsSI("Smolar", "T", T_in, "P", P_in, fluid)
    h2 = h1 - eff * (h1 - CP.PropsSI("Hmolar", "P", P_out, "Smolar", s1, fluid))
    return CP.PropsSI("T", "P", P_out, "Hmolar", h2, fluid), h1 - h2, h2


def compressor(fluid, P_lo, P_hi, T_in=T_AMB, eff=None, r_max=3.0):
    """Intercooled compression, W per mol/s."""
    eff = eff or DB.COMPRESSOR_EFF
    n = max(1, int(np.ceil(np.log(P_hi / P_lo) / np.log(r_max))))
    r = (P_hi / P_lo) ** (1 / n)
    W, P = 0.0, P_lo
    for _ in range(n):
        h1 = CP.PropsSI("Hmolar", "T", T_in, "P", P, fluid)
        s1 = CP.PropsSI("Smolar", "T", T_in, "P", P, fluid)
        W += (CP.PropsSI("Hmolar", "P", P * r, "Smolar", s1, fluid) - h1) / eff
        P *= r
    return W, n


def n2_cycle(P_hi, P_lo, T_e1, T_e2, f1):
    """Two-expander nitrogen Brayton, curves per 1 mol/s of compressor flow."""
    if T_e2 >= T_e1:
        return None
    T1o, w1, h1o = expand("Nitrogen", T_e1, P_hi, P_lo)
    T2o, w2, h2o = expand("Nitrogen", T_e2, P_hi, P_lo)
    # up to 10 % liquid is tolerated at the cold expander exhaust
    hV = CP.PropsSI("Hmolar", "P", P_lo, "Q", 1, "Nitrogen")
    hL = CP.PropsSI("Hmolar", "P", P_lo, "Q", 0, "Nitrogen")
    q2 = (h2o - hL) / (hV - hL)
    if q2 < 0.90 or (h1o - hL) / (hV - hL) < 1.0:
        return None
    wet = h2o if q2 < 1.0 else None
    c = [pure_curve("Nitrogen", "N2 cycle HP", 1.0, P_hi, T_e1, T_AMB, True),
         pure_curve("Nitrogen", "N2 cycle HP (cold leg)", 1 - f1, P_hi, T_e2, T_e1, True),
         pure_curve("Nitrogen", "N2 cycle LP from expander 1", f1, P_lo, T1o, T_WARM_OUT, False),
         pure_curve("Nitrogen", "N2 cycle LP from expander 2", 1 - f1, P_lo, T2o, T_WARM_OUT, False,
                    h_start=wet)]
    W, n = compressor("Nitrogen", P_lo, P_hi)
    return dict(curves=c, W=W, stages=n, T1o=T1o, T2o=T2o, q2=min(q2, 1.0),
                W_exp=f1 * w1 + (1 - f1) * w2)


def he_cycle(P_hi, P_lo, T_e):
    To, w, _ = expand("Helium", T_e, P_hi, P_lo)
    c = [pure_curve("Helium", "He cycle HP", 1.0, P_hi, T_e, T_AMB, True),
         pure_curve("Helium", "He cycle LP", 1.0, P_lo, To, T_WARM_OUT, False)]
    W, n = compressor("Helium", P_lo, P_hi)
    return dict(curves=c, W=W, stages=n, To=To, W_exp=w)


def close_box(proc, n2, he, m_he, m_hi=20000.0, tol=1e3):
    """Smallest N2 flow (mol/s) that leaves no refrigeration unmet."""
    hec = [c.scaled(m_he) for c in he["curves"]] if m_he > 0 else []

    def Qc(m):
        return cascade(proc + hec + [c.scaled(m) for c in n2["curves"]])[1]
    if Qc(m_hi) > tol:
        return None
    lo, hi = 0.0, m_hi
    if Qc(0.0) <= tol:
        return 0.0
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        (hi, lo) = (mid, lo) if Qc(mid) <= tol else (hi, mid)
        if hi - lo < 1e-4 * hi:
            break
    return hi


def evaluate(proc, p):
    n2 = n2_cycle(p["P_hi"] * 1e5, p["P_lo"] * 1e5, p["T_e1"], p["T_e2"], p["f1"])
    he = he_cycle(p["He_P_hi"] * 1e5, p["He_P_lo"] * 1e5, p["He_T_e"])
    if n2 is None:
        return None
    m = close_box(proc, n2, he, p["m_he"])
    if m is None:
        return None
    W_n2, W_he = m * n2["W"], p["m_he"] * he["W"]
    allc = proc + [c.scaled(m) for c in n2["curves"]] + [c.scaled(p["m_he"]) for c in he["curves"]]
    Qh, Qc, Tp = cascade(allc)
    return dict(m_n2=m, W_n2=W_n2, W_he=W_he, W=W_n2 + W_he,
                W_exp=m * n2["W_exp"] + p["m_he"] * he["W_exp"],
                n2=dict(T1o=n2["T1o"], T2o=n2["T2o"], q2=n2["q2"], stages=n2["stages"]),
                he=dict(To=he["To"], stages=he["stages"]),
                Qh=Qh, Qc=Qc, pinch_T=Tp, params=p,
                duties={c.name: float(c.H[-1]) for c in proc})


def composites(curves, n=120):
    """Hot and cold composite curves for plotting: (T, H) lists."""
    out = {}
    for hot in (True, False):
        cs = [c for c in curves if c.hot == hot]
        T = np.unique(np.concatenate([c.T for c in cs]))
        H = np.array([sum(c.q(T[0], t) for c in cs) for t in T])
        out["hot" if hot else "cold"] = dict(T=T.tolist(), H=(H / 1e6).tolist())
    return out
