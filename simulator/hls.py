"""
Helium liquefaction - a Collins-type cycle
==========================================

Interpretation of the user's "self cooling, expander cooling, expander
cooling": two turbo-expanders provide the refrigeration and a final
Joule-Thomson expansion - helium cooling itself on throttling - produces
the liquid.

Pure helium, so CoolProp's reference equation of state is used directly;
near 4 K a cubic equation of state describes helium poorly.

Layout (high-pressure stream, warm to cold):

   300 K  compressor aftercooler
     | HX-A   cooled by the returning low-pressure stream AND by the cold
     |        nitrogen recovered in NRHR (the paper's "cold energy")
   T_A    purified product helium joins here, already cold
     | HX-B
   T_E1   expander 1 takes fraction x1 of the high-pressure flow
     | HX-C
     | HX-D
   T_E2   expander 2 takes fraction x2
     | HX-E
   T_JT   Joule-Thomson valve to the storage pressure -> liquid product

The expander outlets rejoin the low-pressure return at the matching
temperature. The yield of the JT step fixes how much must pass the valve
to make the product; the energy balance of each exchanger then fixes the
low-pressure temperatures, and every exchanger must hold its minimum
approach along its whole length, not just at its ends - helium's heat
capacity changes sharply near 5-10 K, which is where a pinch hides.
"""

import numpy as np
import CoolProp.CoolProp as CP

import design_basis as DB

HE = "Helium"
N2 = "Nitrogen"


def h(T, P, f=HE):
    return CP.PropsSI("Hmolar", "T", T, "P", P, f)


def s(T, P, f=HE):
    return CP.PropsSI("Smolar", "T", T, "P", P, f)


def T_ph(P, H, f=HE):
    return CP.PropsSI("T", "P", P, "Hmolar", H, f)


def expander(T_in, P_in, P_out, eff):
    h1 = h(T_in, P_in)
    s1 = s(T_in, P_in)
    h2s = CP.PropsSI("Hmolar", "P", P_out, "Smolar", s1, HE)
    h2 = h1 - eff * (h1 - h2s)
    return CP.PropsSI("T", "P", P_out, "Hmolar", h2, HE), h1 - h2


def jt_yield(T_in, P_in, P_out):
    """Liquid fraction from throttling, and the flash temperature."""
    H = h(T_in, P_in)
    q = CP.PropsSI("Q", "P", P_out, "Hmolar", H, HE)
    T = CP.PropsSI("T", "P", P_out, "Hmolar", H, HE)
    if q < 0 or q > 1:
        return 0.0, T
    return 1.0 - q, T


def compressor_power(m, P_lo, P_hi, T=300.0, stages=2, eff=None):
    """Intercooled compression of helium, W for m mol/s."""
    eff = eff or DB.COMPRESSOR_EFF
    r = (P_hi / P_lo) ** (1 / stages)
    W, P = 0.0, P_lo
    for _ in range(stages):
        P2 = P * r
        h1, s1 = h(T, P), s(T, P)
        h2s = CP.PropsSI("Hmolar", "P", P2, "Smolar", s1, HE)
        W += m * (h2s - h1) / eff
        P = P2
    return W


def exchanger_pinch(m_hot, P_hot, Th_in, Th_out, m_cold, P_cold, hc_in,
                    extra_cold=None, n=40):
    """
    Counter-current exchanger, hot side fully specified. Returns the cold
    outlet temperature from the energy balance and the smallest
    temperature difference along the exchanger.

    extra_cold: optional (m, P, T_in, T_out, fluid) second cold stream,
    for the nitrogen precooling in HX-A.
    """
    Q = m_hot * (h(Th_in, P_hot) - h(Th_out, P_hot))
    Q_extra = 0.0
    if extra_cold:
        me, Pe, te_in, te_out, fe = extra_cold
        Q_extra = me * (h(te_out, Pe, fe) - h(te_in, Pe, fe))
    # The cold inlet is carried as an enthalpy: the JT flash vapour sits
    # exactly on the saturation line, where T and P do not fix the state.
    hc_out = hc_in + (Q - Q_extra) / m_cold
    try:
        Tc_out = T_ph(P_cold, hc_out)
    except Exception:
        return None, -1.0, Q
    # walk along the exchanger from the cold end
    dmin = np.inf
    Ths = np.linspace(Th_out, Th_in, n)
    for Th in Ths:
        q = m_hot * (h(Th, P_hot) - h(Th_out, P_hot))       # heat from cold end
        qe = 0.0
        if extra_cold:
            me, Pe, te_in, te_out, fe = extra_cold
            # extra stream spans the whole exchanger linearly in duty
            qe = Q_extra * q / Q if Q > 0 else 0.0
        hc = hc_in + (q - qe) / m_cold
        try:
            Tc = T_ph(P_cold, hc)
        except Exception:
            return Tc_out, -1.0, Q
        dmin = min(dmin, Th - Tc)
    return hc_out, dmin, Q


def collins(F_prod, T_feed=70.0, P_feed=25e5, x1=0.35, x2=0.30,
            T_A=80.0, T_E1=60.0, T_E2=15.0, T_JT=6.5, T_D=None,
            P_hi=None, P_lo=None, n2_cold=None):
    """
    Size a Collins liquefier to produce F_prod mol/s of liquid helium.

    n2_cold: (mol/s, T_in K) of recovered nitrogen available to precool
    HX-A, from NRHR.
    Returns a dict with flows, power, and the pinch in each exchanger.
    """
    P_hi = P_hi or DB.HE_CYCLE_HIGH_P_BAR * 1e5
    P_lo = P_lo or DB.HE_CYCLE_LOW_P_BAR * 1e5
    T_D = T_D or 30.0
    dT = DB.HLS_MIN_APPROACH_K

    y, T_flash = jt_yield(T_JT, P_hi, P_lo)
    if y <= 0:
        return dict(ok=False, reason="no liquid from JT")
    m_jt = F_prod / y
    m_cold = m_jt / max(1.0 - x1 - x2, 1e-6)       # HP flow below T_A
    m_c = m_cold - F_prod                            # cycle compressor flow
    if m_c <= 0:
        return dict(ok=False, reason="negative cycle flow")

    # expander outlets
    T_e1o, w1 = expander(T_E1, P_hi, P_lo, DB.EXPANDER_EFF)
    T_e2o, w2 = expander(T_E2, P_hi, P_lo, DB.EXPANDER_EFF)

    # Low-pressure return, cold end upward, carried as enthalpy.
    # Each expander exhaust joins where the return stream is at about its
    # own temperature: expander 2 leaves at ~6.6 K, so it joins at the
    # COLD end of the JT exchanger and supplies the cold flow that
    # exchanger needs; expander 1 leaves near 28 K and joins at the cold
    # end of HX-C, placed there by choosing T_D just above it. Joining
    # either exhaust higher up starves the exchanger below it of cold flow.
    m2 = x2 * m_cold
    m1 = x1 * m_cold
    m_lp = m_jt - F_prod
    h_flash = CP.PropsSI("Hmolar", "P", P_lo, "Q", 1, HE)
    h_lp = (m_lp * h_flash + m2 * h(T_e2o, P_lo)) / (m_lp + m2)
    m_lp += m2
    pinches = {}

    hl, d, _ = exchanger_pinch(m_jt, P_hi, T_E2, T_JT, m_lp, P_lo, h_lp)
    pinches["HX-E"] = d
    if hl is None:
        return dict(ok=False, reason="HX-E energy balance")
    h_lp = hl

    hl, d, _ = exchanger_pinch(m_jt + m2, P_hi, T_D, T_E2, m_lp, P_lo, h_lp)
    pinches["HX-D"] = d
    if hl is None:
        return dict(ok=False, reason="HX-D energy balance")
    h_lp = (m_lp * hl + m1 * h(T_e1o, P_lo)) / (m_lp + m1)
    m_lp += m1

    hl, d, _ = exchanger_pinch(m_jt + m2, P_hi, T_E1, T_D, m_lp, P_lo, h_lp)
    pinches["HX-C"] = d
    if hl is None:
        return dict(ok=False, reason="HX-C energy balance")
    h_lp = hl

    hl, d, _ = exchanger_pinch(m_cold, P_hi, T_A, T_E1, m_lp, P_lo, h_lp)
    pinches["HX-B"] = d
    if hl is None:
        return dict(ok=False, reason="HX-B energy balance")
    h_lp = hl

    extra = None
    if n2_cold and DB.HLS_PRECOOL_BY_N2:
        mn, Tn = n2_cold
        extra = (mn, 1.3e5, Tn, 295.0, N2)
    hl, d, QA = exchanger_pinch(m_c, P_hi, 300.0, T_A, m_lp, P_lo, h_lp,
                                extra_cold=extra)
    pinches["HX-A"] = d
    if hl is None:
        return dict(ok=False, reason="HX-A energy balance")
    Tl = T_ph(P_lo, hl)

    # product helium: 70 K, 25 bar, throttled to P_hi and mixed in at T_A
    # (its own small cooling from T_feed to T_A is folded into HX-B)
    W = compressor_power(m_c, P_lo, P_hi)
    W_exp = m1 * w1 + m2 * w2            # not recovered - braked
    dmin = min(pinches.values())
    return dict(ok=dmin >= dT - 1e-6, dTmin=dmin, pinches=pinches,
                y_jt=y, T_flash=T_flash, m_c=m_c, m_jt=m_jt, x1=x1, x2=x2,
                T_e1o=T_e1o, T_e2o=T_e2o, W=W, W_exp=W_exp,
                T_lp_warm=Tl, Q_A=QA,
                specific_kWh_per_L=W / 1000.0 / (F_prod * 4.0026e-3
                                                 / 0.125 * 3600.0))
