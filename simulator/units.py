"""
Unit operations
===============

Everything in the flowsheet that is not a distillation column. Each unit
takes Stream objects in and returns Stream objects out, so a flowsheet
is just these called in order with a loop around the recycles.

Conventions: component molar flows in mol/s, temperature in K, pressure
in Pa. Aspen states the ASU in kmol/hr, degC and bar, so the flowsheet
layer converts at the boundary and nothing in here has to think about it.
"""

import numpy as np
from flash import (flash_TP, flash_PH, flash_PS, flash_PQ,
                   bubble_T, dew_T, mixture_H, mixture_S)


class Stream:
    """A material stream. Holds component flows and a state."""

    __slots__ = ("pkg", "n", "T", "P", "_beta", "_x", "_y", "name",
                 "forced_beta")

    def __init__(self, pkg, n, T, P, name="", forced_beta=None):
        """
        forced_beta pins the vapour fraction instead of deriving it from
        the temperature and pressure.

        A stream sitting exactly at its bubble or dew point is ambiguous
        from T and P alone: a fraction of a percent of error in the
        vapour pressure flips it between liquid and vapour, and its
        enthalpy then differs by the whole latent heat. Several ASU
        streams sit exactly there - the nitrogen leaving the high
        pressure column is saturated liquid by construction - so when
        the phase is known it is stated rather than rediscovered.
        """
        self.pkg = pkg
        self.n = np.asarray(n, float)
        self.T = float(T)
        self.P = float(P)
        self.name = name
        self.forced_beta = forced_beta
        self._beta = None
        self._x = None
        self._y = None

    # ---- basic quantities -------------------------------------------

    @property
    def F(self):
        """Total molar flow, mol/s."""
        return float(self.n.sum())

    @property
    def z(self):
        """Overall mole fractions."""
        t = self.F
        return self.n / t if t > 0 else self.n.copy()

    @property
    def mass(self):
        """Mass flow, kg/s."""
        return float(self.n @ self.pkg.MW) / 1000.0

    def _flash(self):
        if self._beta is not None:
            return
        if self.F <= 1e-30:
            # An empty stream is legitimate - HPCTSEP produces no vapour
            # because its feed is subcooled, and Aspen's own stream table
            # has no SEPT in it either. Flashing a zero flow would
            # normalise 0/0 and put NaN through the whole flowsheet.
            nc = len(self.n)
            self._beta = 0.0
            self._x = np.zeros(nc)
            self._y = np.zeros(nc)
            return
        b, x, y, _ = flash_TP(self.pkg, self.T, self.P, self.z)
        if self.forced_beta is not None:
            fb = float(self.forced_beta)
            if fb <= 1e-9:
                b, x, y = 0.0, self.z.copy(), y
            elif fb >= 1.0 - 1e-9:
                b, x, y = 1.0, x, self.z.copy()
            else:
                b = fb
        self._beta, self._x, self._y = b, x, y

    @property
    def beta(self):
        self._flash(); return self._beta

    @property
    def x(self):
        self._flash(); return self._x

    @property
    def y(self):
        self._flash(); return self._y

    @property
    def H(self):
        """Total enthalpy flow, W."""
        if self.F <= 1e-30:
            return 0.0
        self._flash()
        return self.F * mixture_H(self.pkg, self.T, self.P,
                                  self._beta, self._x, self._y)

    @property
    def S(self):
        if self.F <= 1e-30:
            return 0.0
        self._flash()
        return self.F * mixture_S(self.pkg, self.T, self.P,
                                  self._beta, self._x, self._y)

    @property
    def Hmol(self):
        if self.F <= 1e-30:
            return 0.0
        self._flash()
        return mixture_H(self.pkg, self.T, self.P, self._beta,
                         self._x, self._y)

    def copy(self, **kw):
        return Stream(self.pkg, kw.get("n", self.n.copy()),
                      kw.get("T", self.T), kw.get("P", self.P),
                      kw.get("name", self.name),
                      kw.get("forced_beta", self.forced_beta))

    def __repr__(self):
        return (f"<{self.name or 'stream'} {self.F*3.6:.1f} kmol/hr "
                f"{self.T-273.15:.2f} C {self.P/1e5:.3f} bar "
                f"VF={self.beta:.3f}>")


# ---------------------------------------------------------------------
# simple units
# ---------------------------------------------------------------------

def mixer(streams, P=None, name=""):
    """Adiabatic mix. Outlet pressure defaults to the lowest inlet."""
    streams = [s for s in streams if s.F > 0]
    if not streams:
        raise ValueError("mixer with no feed")
    pkg = streams[0].pkg
    n = sum(s.n for s in streams)
    Pout = min(s.P for s in streams) if P is None else P
    H = sum(s.H for s in streams)
    T, b, x, y = flash_PH(pkg, Pout, n / n.sum(), H / n.sum())
    return Stream(pkg, n, T, Pout, name)


def splitter(stream, fractions, basis="mole", names=None):
    """Split one stream into several with identical composition."""
    fr = np.asarray(fractions, float)
    if basis == "mass":
        # a mass split is a mole split too when composition is unchanged
        fr = fr / fr.sum()
    else:
        fr = fr / fr.sum()
    out = []
    for i, f in enumerate(fr):
        nm = names[i] if names else ""
        out.append(Stream(stream.pkg, stream.n * f, stream.T, stream.P, nm))
    return out


def separator(stream, keep_mask, names=("keep", "reject")):
    """
    Ideal component separator - the PPU.

    Components flagged False leave in the reject stream entirely. Used
    for the pre-purification unit, which takes CO2 and water out of the
    air before the cold box; Aspen models it the same way.
    """
    keep = np.asarray(keep_mask, bool)
    n_keep = np.where(keep, stream.n, 0.0)
    n_rej = np.where(keep, 0.0, stream.n)
    a = Stream(stream.pkg, n_keep, stream.T, stream.P, names[0])
    b = Stream(stream.pkg, n_rej, stream.T, stream.P, names[1])
    return a, b


def valve(stream, P_out, name=""):
    """Isenthalpic expansion."""
    if stream.F <= 1e-30:
        return Stream(stream.pkg, stream.n.copy(), stream.T, P_out, name)
    T, b, x, y = flash_PH(stream.pkg, P_out, stream.z, stream.Hmol)
    return Stream(stream.pkg, stream.n.copy(), T, P_out, name)


def pump(stream, P_out, efficiency=0.75, name=""):
    """Liquid pump. The temperature rise is negligible here but kept."""
    s = Stream(stream.pkg, stream.n.copy(), stream.T, P_out, name)
    return s


def heater(stream, P_out=None, T_out=None, duty=None, vapour_frac=None,
           name=""):
    """
    Heat exchange against an unspecified utility.

    Give exactly one of outlet temperature, duty, or vapour fraction.
    Returns (stream, duty in W).
    """
    pkg = stream.pkg
    P = stream.P if P_out is None else P_out
    if stream.F <= 1e-30:
        return Stream(pkg, stream.n.copy(),
                      stream.T if T_out is None else T_out, P, name), 0.0
    H_in = stream.H

    if T_out is not None:
        out = Stream(pkg, stream.n.copy(), T_out, P, name)
    elif vapour_frac is not None:
        T, b, x, y = flash_PQ(pkg, P, stream.z, vapour_frac)
        out = Stream(pkg, stream.n.copy(), T, P, name)
    elif duty is not None:
        Hm = (H_in + duty) / stream.F
        T, b, x, y = flash_PH(pkg, P, stream.z, Hm)
        out = Stream(pkg, stream.n.copy(), T, P, name)
    else:
        raise ValueError("heater needs T_out, duty or vapour_frac")
    return out, out.H - H_in


def flash_drum(stream, T=None, P=None, duty=None,
               names=("vapour", "liquid")):
    """
    Two-phase separator.

    Specify temperature and pressure, or pressure and duty. Returns
    (vapour, liquid, duty).
    """
    pkg = stream.pkg
    P = stream.P if P is None else P
    if stream.F <= 1e-30:
        z = Stream(pkg, stream.n.copy(), stream.T, P, names[0])
        z2 = Stream(pkg, stream.n.copy(), stream.T, P, names[1])
        return z, z2, 0.0
    if T is None:
        if duty is None:
            duty = 0.0
        Hm = (stream.H + duty) / stream.F
        T, b, x, y = flash_PH(pkg, P, stream.z, Hm)
    else:
        b, x, y, _ = flash_TP(pkg, T, P, stream.z)

    F = stream.F
    v = Stream(pkg, y * b * F, T, P, names[0])
    l = Stream(pkg, x * (1.0 - b) * F, T, P, names[1])
    Q = (v.H + l.H) - stream.H
    return v, l, Q


def compressor(stream, P_out, isentropic_eff=0.80, mechanical_eff=1.0,
               name=""):
    """
    Adiabatic compression.

    The reversible outlet is found by an isentropic flash, then the real
    outlet by adding the ideal enthalpy rise divided by the efficiency.
    Returns (stream, shaft power in W).
    """
    pkg = stream.pkg
    H_in = stream.Hmol
    S_in = stream.S / stream.F

    Ts, bs, xs, ys = flash_PS(pkg, P_out, stream.z, S_in)
    H_ideal = mixture_H(pkg, Ts, P_out, bs, xs, ys)
    dH_ideal = H_ideal - H_in
    dH_real = dH_ideal / isentropic_eff

    T, b, x, y = flash_PH(pkg, P_out, stream.z, H_in + dH_real)
    out = Stream(pkg, stream.n.copy(), T, P_out, name)
    power = stream.F * dH_real / mechanical_eff
    return out, power


def expander(stream, P_out, isentropic_eff=0.90, mechanical_eff=1.0,
             name=""):
    """
    Adiabatic expansion through a turbine.

    Same structure as the compressor, but the efficiency MULTIPLIES the
    ideal enthalpy drop instead of dividing it - a real turbine recovers
    less work than the reversible one, where a real compressor consumes
    more. Getting that the wrong way round is a silent error: the
    flowsheet still converges, to a plant that is slightly too cold.
    """
    pkg = stream.pkg
    H_in = stream.Hmol
    S_in = stream.S / stream.F

    Ts, bs, xs, ys = flash_PS(pkg, P_out, stream.z, S_in)
    H_ideal = mixture_H(pkg, Ts, P_out, bs, xs, ys)
    dH_ideal = H_ideal - H_in              # negative
    dH_real = dH_ideal * isentropic_eff

    T, b, x, y = flash_PH(pkg, P_out, stream.z, H_in + dH_real)
    out = Stream(pkg, stream.n.copy(), T, P_out, name)
    power = stream.F * dH_real * mechanical_eff   # negative = produced
    return out, power


def multistage_compressor(stream, P_out, n_stages=4, poly_eff=0.80,
                          mech_eff=0.98, intercool_T=303.15, name=""):
    """
    Multistage compression with intercooling, as AIRCOMP.

    Equal pressure ratio per stage, cooled back to the same temperature
    between stages, which is how the block is specified in the .bkp.
    Returns (stream, total shaft power, total cooling duty).
    """
    ratio = (P_out / stream.P) ** (1.0 / n_stages)
    s = stream
    W = 0.0
    Qc = 0.0
    for k in range(n_stages):
        P_k = s.P * ratio
        s, w = compressor(s, P_k, isentropic_eff=poly_eff,
                          mechanical_eff=mech_eff)
        W += w
        s, q = heater(s, T_out=intercool_T)
        Qc += q
    s.name = name
    return s, W, Qc


# ---------------------------------------------------------------------
# heat exchangers
# ---------------------------------------------------------------------

def heatx(hot_in, cold_in, spec, value, P_hot=None, P_cold=None,
          names=("hot_out", "cold_out")):
    """
    Two-stream exchanger, shortcut model.

    Supported specifications, matching the .bkp:
      "DELT-COLD"  : cold outlet approaches the hot inlet within `value` K
      "DEGSUB-HOT" : hot outlet is `value` K below its bubble point
      "T-HOT"      : hot outlet temperature given directly

    The duty follows from whichever stream is specified, and the other
    stream takes the equal and opposite duty.
    """
    pkgh, pkgc = hot_in.pkg, cold_in.pkg
    Ph = hot_in.P if P_hot is None else P_hot
    Pc = cold_in.P if P_cold is None else P_cold

    if spec == "DELT-COLD":
        T_cold_out = hot_in.T - value
        cold_out = Stream(pkgc, cold_in.n.copy(), T_cold_out, Pc, names[1])
        Q = cold_out.H - cold_in.H
        hot_out, _ = heater(hot_in, P_out=Ph, duty=-Q, name=names[0])
    elif spec == "DEGSUB-HOT":
        Tb, _, _ = bubble_T(pkgh, Ph, hot_in.z)
        T_hot_out = Tb - value
        hot_out = Stream(pkgh, hot_in.n.copy(), T_hot_out, Ph, names[0])
        Q = hot_in.H - hot_out.H
        cold_out, _ = heater(cold_in, P_out=Pc, duty=Q, name=names[1])
    elif spec == "T-HOT":
        hot_out = Stream(pkgh, hot_in.n.copy(), value, Ph, names[0])
        Q = hot_in.H - hot_out.H
        cold_out, _ = heater(cold_in, P_out=Pc, duty=Q, name=names[1])
    else:
        raise ValueError(f"unknown HeatX specification {spec!r}")

    return hot_out, cold_out, abs(Q)


def multistream_hx(hot_specs, cold_specs, name=""):
    """
    Multi-stream exchanger, as HX1.

    Each hot side carries its own outlet temperature specification, so
    the hot duties are known without reference to the cold side. Their
    total is then distributed over the cold streams in proportion to
    each one's capacity to absorb it, and the block reports whether the
    overall balance closes.

    hot_specs  : list of (stream, T_out, P_out)
    cold_specs : list of (stream, P_out)
    Returns (hot_outs, cold_outs, duty, imbalance_fraction).
    """
    hot_outs = []
    Q_hot = 0.0
    for s, T_out, P_out in hot_specs:
        P = s.P if P_out is None else P_out
        out = Stream(s.pkg, s.n.copy(), T_out, P, s.name + "_out")
        Q_hot += s.H - out.H
        hot_outs.append(out)

    # Warm every cold stream to the warmest hot inlet, then scale back so
    # the duties match: the exchanger cannot deliver more heat than the
    # hot side gives up.
    T_target = max(s.T for s, _, _ in hot_specs)
    trial = []
    Q_cold_max = 0.0
    for s, P_out in cold_specs:
        P = s.P if P_out is None else P_out
        out = Stream(s.pkg, s.n.copy(), T_target, P, s.name + "_out")
        q = out.H - s.H
        trial.append((s, P, q))
        Q_cold_max += q

    scale = Q_hot / Q_cold_max if Q_cold_max > 0 else 0.0
    cold_outs = []
    for s, P, q in trial:
        out, _ = heater(s, P_out=P, duty=q * scale, name=s.name + "_out")
        cold_outs.append(out)

    Q_cold = sum(o.H for o in cold_outs) - sum(s.H for s, _ in cold_specs)
    imbalance = abs(Q_hot - Q_cold) / max(Q_hot, 1.0)
    return hot_outs, cold_outs, Q_hot, imbalance
