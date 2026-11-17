"""
LNG production: propane precooling, heavy-hydrocarbon knockout, and the
mixed-refrigerant main cryogenic heat exchanger (C3MR, simplified to a
single mixed-refrigerant stream).

The main cryogenic heat exchanger (MCHE) is where the design lives or
dies. The natural gas has to be cooled from -35 to -155 C by a mixed
refrigerant boiling at low pressure, and everywhere along the exchanger
the refrigerant must be at least MIN_APPROACH colder than the streams it
cools. That is checked with composite curves - total heat released by the
hot streams against total heat absorbed by the cold one, both as
functions of temperature - and the smallest vertical gap between them is
the pinch.
"""
import numpy as np
import CoolProp.CoolProp as CP

import design_basis as DB
from props import PKG, C, vec, KMOLH
from flash import flash_TP, flash_PH, flash_PS, mixture_H, mixture_S

K0 = 273.15


def h_molar(T, P, z):
    b, x, y, _ = flash_TP(PKG, T, P, z)
    return mixture_H(PKG, T, P, b, x, y)


def cooling_curve(z, P, T_hi, T_lo, n=40):
    """Molar enthalpy along a stream at constant pressure, T_hi -> T_lo."""
    Ts = np.linspace(T_hi, T_lo, n)
    return Ts, np.array([h_molar(T, P, z) for T in Ts])


# ---------------------------------------------------------------------
# composite curves
# ---------------------------------------------------------------------

def composite(streams, T_grid):
    """
    Sum of heat flows of several streams over a common temperature grid.

    streams: list of (F mol/s, Ts ascending, H ascending, T_min, T_max)
    Returns cumulative heat flow (W) relative to the grid's cold end.
    """
    Q = np.zeros_like(T_grid)
    for F, Ts, Hs, tmin, tmax in streams:
        Tc = np.clip(T_grid, tmin, tmax)
        Q += F * (np.interp(Tc, Ts, Hs) - np.interp(tmin, Ts, Hs))
    return Q


def pinch(hot, cold, T_lo, T_hi, n=300):
    """
    Minimum temperature approach between hot and cold composites.

    Each composite is built only over the temperature span its streams
    actually cover. Extending one past its end leaves a flat stretch of
    constant heat, and interpolating temperature on a flat stretch is
    undefined - it produced a false 0 K pinch at the warm end.
    """
    th_lo = min(s[3] for s in hot); th_hi = max(s[4] for s in hot)
    tc_lo = min(s[3] for s in cold); tc_hi = max(s[4] for s in cold)
    Tgh = np.linspace(th_lo, th_hi, n)
    Tgc = np.linspace(tc_lo, tc_hi, n)
    Qh = composite(hot, Tgh)
    Qc = composite(cold, Tgc)
    Qmax = min(Qh[-1], Qc[-1])
    q = np.linspace(0, Qmax, n)
    Th = np.interp(q, Qh, Tgh)
    Tc = np.interp(q, Qc, Tgc)
    return float(np.min(Th - Tc)), Qh[-1], Qc[-1], (q, Th, Tc)


def stream_curve(z, P, T_lo, T_hi, n=45):
    Ts = np.linspace(T_lo, T_hi, n)
    Hs = np.array([h_molar(T, P, z) for T in Ts])
    return Ts, Hs


# ---------------------------------------------------------------------
# the main cryogenic heat exchanger
# ---------------------------------------------------------------------

def mche(ng_n, mr_z, P_ng=None, P_hi=None, P_lo=None, T_in=None,
         T_ko=None, T_out=None, dT_warm=None, n=36):
    """
    Size the mixed-refrigerant flow for the MCHE and check its pinch.

    ng_n  natural gas component flows entering the MCHE, mol/s
    mr_z  mixed refrigerant composition

    The refrigerant flow comes from the energy balance: the cold
    low-pressure refrigerant must absorb everything the natural gas and
    the high-pressure refrigerant give up. Its warm outlet is fixed at
    dT_warm below the hot inlets. Whether that flow actually works is
    then a separate question - the pinch.
    """
    P_ng = P_ng or (DB.FEED_P_BAR - 3.0) * 1e5
    P_hi = P_hi or (DB.MR_HIGH_P_BAR - 1.0) * 1e5
    P_lo = P_lo or DB.MR_LOW_P_BAR * 1e5
    T_in = T_in or DB.PRECOOL_T_C + K0
    T_ko = T_ko or DB.HHC_SEPARATOR_T_C + K0
    T_out = T_out or DB.MCHE_OUTLET_T_C + K0
    dT_warm = dT_warm if dT_warm is not None else DB.MCHE_MIN_APPROACH_K
    mr_z = np.asarray(mr_z, float); mr_z = mr_z / mr_z.sum()

    F_ng = ng_n.sum(); z_ng = ng_n / F_ng

    # heavy-hydrocarbon knockout partway down
    b, x, y, _ = flash_TP(PKG, T_ko, P_ng, z_ng)
    F_v = b * F_ng
    hhc = x * (1 - b) * F_ng              # heavy liquid, mol/s by component

    # stream enthalpy curves (ascending T)
    TA, HA = stream_curve(z_ng, P_ng, T_ko, T_in, n // 2)
    TB, HB = stream_curve(y, P_ng, T_out, T_ko, n)
    TM, HM = stream_curve(mr_z, P_hi, T_out, T_in, n)

    # JT expansion of the subcooled refrigerant
    h_jt = np.interp(T_out, TM, HM)
    T_jt, *_ = flash_PH(PKG, P_lo, mr_z, h_jt, lo=60.0, hi=T_in)
    T_warm = T_in - dT_warm
    TL, HL = stream_curve(mr_z, P_lo, T_jt, T_warm, n)

    # energy balance
    Q_ng = F_ng * (HA[-1] - HA[0]) + F_v * (HB[-1] - HB[0])
    refrig = (HL[-1] - HL[0]) - (HM[-1] - HM[0])     # per mol MR
    if refrig <= 0:
        return dict(feasible=False, reason="no net refrigeration")
    F_mr = Q_ng / refrig

    hot = [(F_ng, TA, HA, T_ko, T_in), (F_v, TB, HB, T_out, T_ko),
           (F_mr, TM, HM, T_out, T_in)]
    cold = [(F_mr, TL, HL, T_jt, T_warm)]
    dTmin, Qh, Qc, curves = pinch(hot, cold, T_jt, T_in)

    lng = y * F_v
    return dict(feasible=dTmin >= DB.MCHE_MIN_APPROACH_K - 1e-6,
                dTmin=dTmin, F_mr=F_mr, Q=Qh, T_jt=T_jt, T_warm=T_warm,
                lng=lng, hhc=hhc, curves=curves, P_lo=P_lo, P_hi=P_hi,
                mr_z=mr_z)


# ---------------------------------------------------------------------
# power
# ---------------------------------------------------------------------

def compress(z, F, T_in, P_in, P_out, stages, eff=None, T_inter=None):
    """
    Multistage compression with intercooling back to T_inter.
    Returns (shaft power W, after-cooler duty W, discharge T before the
    final cooler). Stage pressure ratios are equal.
    """
    eff = eff or DB.COMPRESSOR_EFF
    T_inter = (T_inter or DB.COOLER_T_C) + K0 if T_inter is None or T_inter < 200 else T_inter
    r = (P_out / P_in) ** (1.0 / stages)
    T, P, W, Qc = T_in, P_in, 0.0, 0.0
    for k in range(stages):
        b, x, y, _ = flash_TP(PKG, T, P, z)
        h1 = mixture_H(PKG, T, P, b, x, y); s1 = mixture_S(PKG, T, P, b, x, y)
        P2 = P * r
        Ts, bs, xs, ys = flash_PS(PKG, P2, z, s1)
        h2s = mixture_H(PKG, Ts, P2, bs, xs, ys)
        h2 = h1 + (h2s - h1) / eff
        T2, *_ = flash_PH(PKG, P2, z, h2)
        W += F * (h2 - h1)
        hc = h_molar(T_inter, P2, z)
        Qc += F * (h2 - hc)
        T_dis = T2
        T, P = T_inter, P2
    return W, Qc, T_dis


PROPANE = "Propane"


def propane_cycle(loads, T_evap_C=(5.0, -15.0, -40.0), T_cond_C=None,
                  eff=None):
    """
    Three-level propane refrigeration with flash economisers.

    loads: heat absorbed at each evaporator level, W, warmest first.
    Liquid cascades down through the levels, flashing at each; the vapour
    from every level is compressed back up in a single machine with side
    inlets. Pure propane, so CoolProp's reference equation of state is
    used directly. Returns shaft power, W.
    """
    eff = eff or DB.COMPRESSOR_EFF
    T_cond = (T_cond_C if T_cond_C is not None else DB.PROPANE_COND_T_C) + K0
    Tev = [t + K0 for t in T_evap_C]
    P = lambda T: CP.PropsSI("P", "T", T, "Q", 0, PROPANE)
    hl = lambda T: CP.PropsSI("H", "T", T, "Q", 0, PROPANE)
    hv = lambda T: CP.PropsSI("H", "T", T, "Q", 1, PROPANE)
    sv = lambda T: CP.PropsSI("S", "T", T, "Q", 1, PROPANE)

    n = len(Tev)
    # work from the coldest level upward: liquid needed at each level
    m_evap = [loads[k] / (hv(Tev[k]) - hl(Tev[k])) for k in range(n)]
    m_liq_down = 0.0
    vap = [0.0] * n
    for k in range(n - 1, -1, -1):
        h_in = hl(Tev[k - 1]) if k > 0 else hl(T_cond)
        x = (h_in - hl(Tev[k])) / (hv(Tev[k]) - hl(Tev[k]))
        need_liq = m_evap[k] + m_liq_down
        m_in = need_liq / (1.0 - x)
        vap[k] = m_evap[k] + x * m_in
        m_liq_down = m_in

    # compression: coldest vapour up to the next level, mix, repeat
    W = 0.0
    m, h = 0.0, 0.0
    levels_P = [P(t) for t in Tev] + [P(T_cond)]
    for k in range(n - 1, -1, -1):
        # mix the stream arriving from below with this level's vapour
        h = (m * h + vap[k] * hv(Tev[k])) / (m + vap[k]) if m > 0 else hv(Tev[k])
        m += vap[k]
        P_in = levels_P[k]
        P_out = levels_P[k - 1] if k > 0 else levels_P[-1]
        s_in = CP.PropsSI("S", "P", P_in, "H", h, PROPANE)
        h_s = CP.PropsSI("H", "P", P_out, "S", s_in, PROPANE)
        h_out = h + (h_s - h) / eff
        W += m * (h_out - h)
        h = h_out
    return W


def precool_loads(z, F, P, T_from, T_steps=(8.0, -12.0, -35.0)):
    """Split a stream's precooling into the three propane levels, W."""
    out, T = [], T_from
    for t in T_steps:
        out.append(F * (h_molar(T, P, z) - h_molar(t + K0, P, z)))
        T = t + K0
    return out


def liquefaction(ng_n, mr_z, P_lo, n=30):
    """Total C3MR power for a given refrigerant, with feasibility flags."""
    P_ng = (DB.FEED_P_BAR - 3.0) * 1e5
    P_hi = DB.MR_HIGH_P_BAR * 1e5
    r = mche(ng_n, mr_z, P_lo=P_lo, n=n)
    if not r.get("F_mr"):
        return dict(r, W=np.inf)
    z = r["mr_z"]; F = r["F_mr"]

    # compressor suction must be all vapour
    b, *_ = flash_TP(PKG, r["T_warm"], P_lo, z)
    W_mr, Qc_mr, T_dis = compress(z, F, r["T_warm"], P_lo, P_hi,
                                  DB.MR_STAGES)
    ng_z = ng_n / ng_n.sum()
    loads = np.array(precool_loads(ng_z, ng_n.sum(), DB.FEED_P_BAR * 1e5,
                                   DB.FEED_T_C + K0))
    loads += np.array(precool_loads(z, F, P_hi, DB.COOLER_T_C + K0))
    W_c3 = propane_cycle(list(loads))
    return dict(r, W=W_mr + W_c3, W_mr=W_mr, W_c3=W_c3,
                suction_vf=b, precool_loads=loads, Qc_mr=Qc_mr)
