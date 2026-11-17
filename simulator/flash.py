"""
Phase equilibrium flashes
=========================

Every unit operation outside a column reduces to one of these. A valve
is an isenthalpic flash, a compressor is an isentropic flash followed by
an efficiency correction, a flash drum is a temperature-pressure or a
duty flash, and a heat exchanger is a pair of enthalpy flashes.

The temperature-pressure flash is the primitive: successive substitution
on the K-values with Rachford-Rice for the phase split, then a few
Newton steps once it is close. Successive substitution is slow but very
hard to break; Newton is fast but needs a good starting point. Using
one to feed the other gets both properties.
"""

import numpy as np
from scipy.optimize import brentq

R = 8.31446261815324


def wilson_K(pkg, T, P):
    """Wilson correlation - a starting estimate, never an answer."""
    Tr = T / pkg.Tc
    return (pkg.Pc / P) * np.exp(5.373 * (1.0 + pkg.omega) * (1.0 - 1.0 / Tr))


def rachford_rice(z, K, tol=1e-13):
    """
    Solve sum z(K-1)/(1+beta(K-1)) = 0 for the vapour fraction beta.

    Bracketed on the interval where every liquid and vapour mole
    fraction stays positive, which keeps the solve well behaved even
    when some K-values are extreme.
    """
    z = np.asarray(z, float)
    K = np.asarray(K, float)
    if np.all(K <= 1.0):
        return 0.0
    if np.all(K >= 1.0):
        return 1.0

    def f(b):
        return float(np.sum(z * (K - 1.0) / (1.0 + b * (K - 1.0))))

    lo = 1.0 / (1.0 - K.max()) + 1e-12
    hi = 1.0 / (1.0 - K.min()) - 1e-12
    lo = max(lo, 0.0)
    hi = min(hi, 1.0)
    if hi <= lo:
        return 0.0 if f(0.0) < 0 else 1.0
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return 0.0 if flo < 0 else 1.0
    return brentq(f, lo, hi, xtol=tol)


def stability(pkg, T, P, z, max_iter=200):
    """
    Michelsen tangent-plane stability test.

    Tries a liquid-like and a vapour-like trial phase. If either lowers
    the Gibbs energy of the feed (sum of trial mole numbers > 1 at the
    stationary point), the single phase is unstable and the returned
    K-values point toward the second phase.

    Needed because successive substitution alone can settle on the
    trivial solution, K = 1, even where a second phase is overwhelmingly
    favoured - which is what happened for the helium-rich vapour at
    70-80 K, where liquid nitrogen plainly has to condense.
    """
    Ta, Pa = np.array([T]), np.array([P])
    z = np.asarray(z, float)
    phz_v, _ = pkg.phi(Ta, Pa, z[None, :], True)
    phz_l, _ = pkg.phi(Ta, Pa, z[None, :], False)
    # the feed's own reference: take the phase with lower Gibbs energy
    gv = np.sum(z * np.log(np.maximum(phz_v[0], 1e-300)))
    gl = np.sum(z * np.log(np.maximum(phz_l[0], 1e-300)))
    lnphz = np.log(phz_v[0] if gv <= gl else phz_l[0])
    d = np.log(np.maximum(z, 1e-300)) + lnphz

    Kw = np.clip(wilson_K(pkg, T, P), 1e-12, 1e12)
    trials = [("liquid", z / Kw, False), ("vapour", z * Kw, True)]
    # a heavy-component-rich liquid seed catches cases Wilson misses
    heavy = np.full_like(z, 1e-6); heavy[np.argmax(pkg.Tc * (z > 0))] = 1.0
    trials.append(("liquid-heavy", heavy, False))

    best = None
    for name, W, vap in trials:
        W = np.maximum(W, 1e-300)
        for _ in range(max_iter):
            x = W / W.sum()
            ph, _ = pkg.phi(Ta, Pa, x[None, :], vap)
            Wn = np.exp(d - np.log(np.maximum(ph[0], 1e-300)))
            if np.max(np.abs(np.log(Wn / W))) < 1e-10:
                W = Wn; break
            W = Wn
        tm = 1.0 - W.sum()                 # < 0 means unstable
        x = W / W.sum()
        if not np.all(np.isfinite(x)):
            continue
        if np.max(np.abs(x - z)) < 1e-6:
            continue                       # collapsed onto the feed
        if tm < -1e-9 and (best is None or tm < best[0]):
            K = (z / x) if not vap else (x / z)
            best = (tm, np.clip(K, 1e-14, 1e14))
    return (best is not None), (best[1] if best else None)


def flash_TP(pkg, T, P, z, K0=None, max_iter=300, tol=1e-11, _retry=True):
    """
    Isothermal flash.

    Returns (beta, x, y, K) where beta is the vapour mole fraction.
    A single-phase result comes back with beta at exactly 0 or 1 and
    both compositions equal to the feed.
    """
    z = np.asarray(z, float)
    z = z / z.sum()
    Ta, Pa = np.array([T]), np.array([P])

    K = wilson_K(pkg, T, P) if K0 is None else np.array(K0, float)
    K = np.clip(K, 1e-12, 1e12)

    beta = 0.5
    for it in range(max_iter):
        beta = rachford_rice(z, K)
        den = 1.0 + beta * (K - 1.0)
        x = z / den
        y = K * x
        sx, sy = x.sum(), y.sum()
        x = x / sx
        y = y / sy

        phil, Zl = pkg.phi(Ta, Pa, x[None, :], vapour=False)
        phiv, Zv = pkg.phi(Ta, Pa, y[None, :], vapour=True)
        Knew = np.clip(phil[0] / phiv[0], 1e-14, 1e14)

        err = np.max(np.abs(np.log(Knew / K)))
        K = Knew
        if err < tol:
            break

    beta = rachford_rice(z, K)
    den = 1.0 + beta * (K - 1.0)
    x = z / den; x /= x.sum()
    y = K * z / den; y /= y.sum()

    # A trivial solution - both phases identical - may mean one phase, or
    # may mean successive substitution collapsed. Test stability first.
    if np.max(np.abs(np.log(np.clip(K, 1e-300, None)))) < 1e-7 and _retry:
        unstable, Ks = stability(pkg, T, P, z)
        if unstable:
            return flash_TP(pkg, T, P, z, K0=Ks, max_iter=max_iter,
                            tol=tol, _retry=False)
    if np.max(np.abs(np.log(np.clip(K, 1e-300, None)))) < 1e-7:
        phil, Zl = pkg.phi(Ta, Pa, z[None, :], vapour=False)
        phiv, Zv = pkg.phi(Ta, Pa, z[None, :], vapour=True)
        beta = 1.0 if Zv[0] > 0.35 else 0.0
        return beta, z.copy(), z.copy(), K

    if beta <= 1e-10:
        return 0.0, z.copy(), y, K
    if beta >= 1.0 - 1e-10:
        return 1.0, x, z.copy(), K
    return beta, x, y, K


def mixture_H(pkg, T, P, beta, x, y):
    Ta, Pa = np.array([T]), np.array([P])
    h = 0.0
    if beta < 1.0:
        h += (1.0 - beta) * pkg.H(Ta, Pa, x[None, :], vapour=False)[0]
    if beta > 0.0:
        h += beta * pkg.H(Ta, Pa, y[None, :], vapour=True)[0]
    return h


def mixture_S(pkg, T, P, beta, x, y):
    Ta, Pa = np.array([T]), np.array([P])
    s = 0.0
    if beta < 1.0:
        s += (1.0 - beta) * pkg.S(Ta, Pa, x[None, :], vapour=False)[0]
    if beta > 0.0:
        s += beta * pkg.S(Ta, Pa, y[None, :], vapour=True)[0]
    return s


def _solve_T(pkg, P, z, target, prop, lo=55.0, hi=400.0, tol=1e-8):
    """Find the temperature at which a flashed property hits a target."""
    def f(T):
        beta, x, y, _ = flash_TP(pkg, T, P, z)
        val = (mixture_H if prop == "H" else mixture_S)(pkg, T, P, beta, x, y)
        return val - target

    flo, fhi = f(lo), f(hi)
    n = 0
    while flo * fhi > 0 and n < 12:
        lo = max(40.0, lo - 40.0); hi = min(900.0, hi + 120.0)
        flo, fhi = f(lo), f(hi); n += 1
    if flo * fhi > 0:
        raise ValueError(f"cannot bracket {prop} = {target:.6g} "
                         f"at P = {P/1e5:.3f} bar")
    T = brentq(f, lo, hi, xtol=tol)
    beta, x, y, _ = flash_TP(pkg, T, P, z)
    return T, beta, x, y


def flash_PH(pkg, P, z, H, **kw):
    """Isenthalpic flash: a valve, or any adiabatic pressure change."""
    return _solve_T(pkg, P, z, H, "H", **kw)


def flash_PS(pkg, P, z, S, **kw):
    """Isentropic flash: the reversible step inside a compressor."""
    return _solve_T(pkg, P, z, S, "S", **kw)


def _saturation_T(pkg, P, z, dew, lo=55.0, hi=None):
    """
    Bubble (dew=False) or dew (dew=True) temperature.

    Solved on sum(z K) = 1 or sum(z / K) = 1 rather than on the vapour
    fraction. At a vapour fraction of exactly zero or one the residual
    is flat over a whole temperature range, so a root finder lands on
    whichever bracket edge it started from instead of the real answer.
    These two equations cross zero cleanly.
    """
    z = np.asarray(z, float); z = z / z.sum()
    Pa = np.array([P])
    hi = min(float(pkg.Tc.max()) * 0.999, 400.0) if hi is None else hi

    def resid(T):
        Ta = np.array([T])
        K = wilson_K(pkg, T, P)
        xx = z.copy() if not dew else z / np.clip(K, 1e-300, None)
        yy = z * K if not dew else z.copy()
        for _ in range(400):
            xx = xx / xx.sum(); yy = yy / yy.sum()
            phil, Zl = pkg.phi(Ta, Pa, xx[None, :], vapour=False)
            phiv, Zv = pkg.phi(Ta, Pa, yy[None, :], vapour=True)
            Kn = np.clip(phil[0] / phiv[0], 1e-300, 1e300)
            if dew:
                xn = z / Kn
                if np.max(np.abs(xn / xn.sum() - xx)) < 1e-13:
                    xx = xn / xn.sum(); K = Kn; break
                xx = xn / xn.sum()
            else:
                yn = z * Kn
                if np.max(np.abs(yn / yn.sum() - yy)) < 1e-13:
                    yy = yn / yn.sum(); K = Kn; break
                yy = yn / yn.sum()
            K = Kn
        phil, Zl = pkg.phi(Ta, Pa, (z if not dew else xx)[None, :], vapour=False)
        phiv, Zv = pkg.phi(Ta, Pa, (yy if not dew else z)[None, :], vapour=True)
        if abs(Zv[0] - Zl[0]) < 1e-8:
            # Single root - no phase split here. The two residuals run in
            # OPPOSITE directions, so the sign has to follow the one
            # being solved: sum(zK)-1 is negative below the bubble point,
            # while sum(z/K)-1 is POSITIVE below the dew point.
            liquid_like = Zl[0] < 0.35
            if dew:
                return 1e3 if liquid_like else -1e3
            return -1e3 if liquid_like else 1e3
        K = phil[0] / phiv[0]
        return float(np.sum(z / K) - 1.0) if dew else float(np.sum(z * K) - 1.0)

    flo, fhi = resid(lo), resid(hi)
    if flo * fhi > 0:
        raise ValueError(f"no {'dew' if dew else 'bubble'} point between "
                         f"{lo:.1f} and {hi:.1f} K at {P/1e5:.3f} bar")
    return brentq(resid, lo, hi, xtol=1e-9)


def bubble_T(pkg, P, z, **kw):
    T = _saturation_T(pkg, P, z, dew=False, **kw)
    beta, x, y, _ = flash_TP(pkg, T, P, z)
    return T, x, y


def dew_T(pkg, P, z, **kw):
    T = _saturation_T(pkg, P, z, dew=True, **kw)
    beta, x, y, _ = flash_TP(pkg, T, P, z)
    return T, x, y


def flash_PQ(pkg, P, z, vapour_frac, lo=55.0, hi=None):
    """Flash to a specified vapour fraction."""
    hi = min(float(pkg.Tc.max()) * 0.999, 400.0) if hi is None else hi
    if vapour_frac <= 1e-12:
        T = _saturation_T(pkg, P, z, dew=False, lo=lo, hi=hi)
    elif vapour_frac >= 1.0 - 1e-12:
        T = _saturation_T(pkg, P, z, dew=True, lo=lo, hi=hi)
    else:
        Tb = _saturation_T(pkg, P, z, dew=False, lo=lo, hi=hi)
        Td = _saturation_T(pkg, P, z, dew=True, lo=lo, hi=hi)

        def f(T):
            beta, x, y, _ = flash_TP(pkg, T, P, z)
            return beta - vapour_frac

        T = brentq(f, Tb, Td, xtol=1e-9)
    beta, x, y, _ = flash_TP(pkg, T, P, z)
    return T, beta, x, y
