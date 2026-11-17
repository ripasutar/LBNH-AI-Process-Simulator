"""
General cubic equation of state, vectorised over stages
=======================================================

Peng-Robinson or Soave-Redlich-Kwong, for any set of components, with a
choice of temperature function for the attraction term:

  soave            the standard form
  boston-mathias   an exponential extrapolation for components above
                   their critical temperature

Why the choice matters here. This plant carries helium, whose critical
temperature is 5.2 K, through streams at 80 to 300 K - reduced
temperatures of 15 to 60. Helium also has a NEGATIVE acentric factor,
and with a negative kappa the standard Soave function increases with
temperature instead of decaying, so helium would be predicted to grow
stickier the hotter it gets. The Boston-Mathias form fixes the shape
above Tc, but it still inherits the sign of kappa. So for helium the
package pins kappa at zero - the quantum gas is treated as a simple
fluid, which is the common practice in cryogenic process work - and the
binary parameters that decide helium's solubility are fitted to the
GERG-2008 reference equation instead of taken from a table that has
none.

Shared structure with the air separation package: every method takes
T (n,), P (n,), z (n, nc) and returns arrays with a leading stage axis.
"""

import numpy as np

R = 8.31446261815324

EOS_CONST = {
    # u, w in  P = RT/(v-b) - a/(v^2 + u b v + w b^2)
    "PR":  dict(u=2.0, w=-1.0, oa=0.45723552892138218938,
                ob=0.07779607390388846387),
    "SRK": dict(u=1.0, w=0.0, oa=0.42748023354034140439,
                ob=0.08664034996495772158),
}


class Cubic:
    def __init__(self, names, Tc, Pc, omega, MW, cp_poly, kij=None,
                 eos="PR", alpha="boston-mathias", kappa_override=None):
        self.names = list(names)
        self.nc = len(names)
        self.Tc = np.asarray(Tc, float)
        self.Pc = np.asarray(Pc, float)
        self.omega = np.asarray(omega, float)
        self.MW = np.asarray(MW, float)
        self.cp = np.asarray(cp_poly, float)          # (nc, 4)
        self.eos = eos
        self.alpha_form = alpha
        k = EOS_CONST[eos]
        self.u, self.w = k["u"], k["w"]
        self.delta = np.sqrt(self.u ** 2 - 4 * self.w)
        self.d1 = (self.u + self.delta) / 2.0
        self.d2 = (self.u - self.delta) / 2.0

        w_ = self.omega
        if eos == "PR":
            self.kappa = 0.37464 + 1.54226 * w_ - 0.26992 * w_ * w_
        else:
            self.kappa = 0.480 + 1.574 * w_ - 0.176 * w_ * w_
        if kappa_override:
            for nm, val in kappa_override.items():
                self.kappa[self.names.index(nm)] = val

        self.ac = k["oa"] * R ** 2 * self.Tc ** 2 / self.Pc
        self.b_i = k["ob"] * R * self.Tc / self.Pc
        self.kij = (np.zeros((self.nc, self.nc)) if kij is None
                    else np.asarray(kij, float).copy())

    # ---------------------------------------------------------------

    def _alpha(self, T):
        Tr = T[:, None] / self.Tc[None, :]
        m = self.kappa[None, :]
        sq = np.sqrt(Tr)
        s = 1.0 + m * (1.0 - sq)
        a_sub = s * s
        da_sub = -m * s / (sq * self.Tc[None, :])

        if self.alpha_form != "boston-mathias":
            return a_sub, da_sub

        # Boston-Mathias above Tc: alpha = exp[2c(1 - Tr^d)], d = 1 + m/2,
        # c = 1 - 1/d. It matches the Soave value and slope at Tc and decays
        # smoothly beyond it instead of turning upward.
        d = 1.0 + m / 2.0
        c = np.where(np.abs(d) > 1e-12, 1.0 - 1.0 / d, 0.0)
        Trd = Tr ** d
        a_sup = np.exp(2.0 * c * (1.0 - Trd))
        da_sup = a_sup * (-2.0 * c * d * Trd / T[:, None])
        sup = Tr > 1.0
        return np.where(sup, a_sup, a_sub), np.where(sup, da_sup, da_sub)

    def _mix(self, T, z):
        alpha, dalpha = self._alpha(T)
        a_i = self.ac[None, :] * alpha
        da_i = self.ac[None, :] * dalpha
        sa = np.sqrt(np.maximum(a_i, 1e-300))
        one_k = 1.0 - self.kij[None, :, :]
        a_ij = sa[:, :, None] * sa[:, None, :] * one_k
        ratio = sa[:, None, :] / sa[:, :, None]
        da_ij = 0.5 * one_k * (ratio * da_i[:, :, None]
                               + da_i[:, None, :] / ratio)
        a_part = np.einsum("nj,nij->ni", z, a_ij)
        a = np.einsum("ni,ni->n", z, a_part)
        da = np.einsum("ni,nj,nij->n", z, z, da_ij)
        b = z @ self.b_i
        return a, b, da, a_part

    def _roots(self, A, B, vapour):
        u, w = self.u, self.w
        c2 = -(1.0 + B - u * B)
        c1 = A + w * B * B - u * B - u * B * B
        c0 = -(A * B + w * B * B + w * B ** 3)
        shift = c2 / 3.0
        p = c1 - c2 * c2 / 3.0
        q = 2.0 * c2 ** 3 / 27.0 - c2 * c1 / 3.0 + c0
        disc = (q / 2.0) ** 2 + (p / 3.0) ** 3
        Z = np.empty_like(A)
        three = disc < 0
        if np.any(three):
            pt, qt = p[three], q[three]
            mm = 2.0 * np.sqrt(np.maximum(-pt / 3.0, 0.0))
            arg = np.clip(3.0 * qt / (pt * mm + 1e-300), -1.0, 1.0)
            th = np.arccos(arg) / 3.0
            r = np.stack([mm * np.cos(th - 2 * np.pi * k / 3) for k in range(3)],
                         -1) - shift[three, None]
            if vapour:
                Z[three] = r.max(-1)
            else:
                Bt = B[three, None]
                v = np.where(r > Bt * 1.0000001, r, np.inf).min(-1)
                Z[three] = np.where(np.isfinite(v), v, r.max(-1))
        one = ~three
        if np.any(one):
            s = np.sqrt(np.maximum(disc[one], 0.0))
            Z[one] = np.cbrt(-q[one] / 2 + s) + np.cbrt(-q[one] / 2 - s) - shift[one]
        return np.maximum(Z, B * 1.0000001)

    # ---------------------------------------------------------------

    def phi(self, T, P, z, vapour):
        T = np.atleast_1d(np.asarray(T, float))
        P = np.atleast_1d(np.asarray(P, float))
        z = np.atleast_2d(np.asarray(z, float))
        a, b, _, ap = self._mix(T, z)
        A = a * P / (R * T) ** 2
        B = b * P / (R * T)
        Z = self._roots(A, B, vapour)
        bb = self.b_i[None, :] / b[:, None]
        L = np.log((Z + self.d1 * B) / (Z + self.d2 * B))[:, None]
        lnphi = (bb * (Z - 1.0)[:, None]
                 - np.log(np.maximum(Z - B, 1e-300))[:, None]
                 - (A / (self.delta * B))[:, None] * (2.0 * ap / a[:, None] - bb) * L)
        return np.exp(np.clip(lnphi, -700, 700)), Z

    def K(self, T, P, x, y):
        pl, _ = self.phi(T, P, x, False)
        pv, _ = self.phi(T, P, y, True)
        return pl / pv

    def _h_ig(self, T):
        c = self.cp
        Tr = 298.15
        def I(t):
            return (c[:, 0] * t + c[:, 1] * t ** 2 / 2
                    + c[:, 2] * t ** 3 / 3 + c[:, 3] * t ** 4 / 4)
        return I(T[:, None]) - I(np.array([[Tr]]))

    def _s_ig(self, T):
        c = self.cp
        Tr = 298.15
        def J(t):
            return (c[:, 0] * np.log(t) + c[:, 1] * t
                    + c[:, 2] * t ** 2 / 2 + c[:, 3] * t ** 3 / 3)
        return J(T[:, None]) - J(np.array([[Tr]]))

    def H(self, T, P, z, vapour):
        T = np.atleast_1d(np.asarray(T, float))
        P = np.atleast_1d(np.asarray(P, float))
        z = np.atleast_2d(np.asarray(z, float))
        h_ig = np.einsum("ni,ni->n", z, self._h_ig(T))
        a, b, da, _ = self._mix(T, z)
        A = a * P / (R * T) ** 2
        B = b * P / (R * T)
        Z = self._roots(A, B, vapour)
        h_dep = (R * T * (Z - 1.0)
                 + (T * da - a) / (self.delta * b)
                 * np.log((Z + self.d1 * B) / (Z + self.d2 * B)))
        return h_ig + h_dep

    def S(self, T, P, z, vapour):
        T = np.atleast_1d(np.asarray(T, float))
        P = np.atleast_1d(np.asarray(P, float))
        z = np.atleast_2d(np.asarray(z, float))
        s_ig = np.einsum("ni,ni->n", z, self._s_ig(T)) - R * np.log(P / 1e5)
        zz = np.clip(z, 1e-300, 1.0)
        s_ig = s_ig - R * np.einsum("ni,ni->n", z, np.log(zz))
        a, b, da, _ = self._mix(T, z)
        A = a * P / (R * T) ** 2
        B = b * P / (R * T)
        Z = self._roots(A, B, vapour)
        s_dep = (R * np.log(np.maximum(Z - B, 1e-300))
                 + da / (self.delta * b)
                 * np.log((Z + self.d1 * B) / (Z + self.d2 * B)))
        return s_ig + s_dep


# ---------------------------------------------------------------------
# the LBNH component set
# ---------------------------------------------------------------------

COMPONENTS = ["He", "CH4", "N2", "C2H6", "C3H8", "iC4", "nC4", "iC5", "nC5"]
_NAMES = {"He": "helium", "CH4": "methane", "N2": "nitrogen",
          "C2H6": "ethane", "C3H8": "propane", "iC4": "isobutane",
          "nC4": "n-butane", "iC5": "isopentane", "nC5": "n-pentane",
          "H2O": "water"}


def _props(comps):
    from chemicals import CAS_from_any, Tc, Pc, omega, MW
    from thermo import HeatCapacityGas
    cas = [CAS_from_any(_NAMES[c]) for c in comps]
    Ts = np.linspace(60, 340, 40)
    cp = []
    for c in cas:
        hc = HeatCapacityGas(CASRN=c)
        v = np.array([hc(t) for t in Ts])
        cp.append(np.polyfit(Ts, v, 3)[::-1])       # c0 + c1 T + c2 T^2 + c3 T^3
    return (cas, [Tc(c) for c in cas], [Pc(c) for c in cas],
            [omega(c) for c in cas], [MW(c) for c in cas], np.array(cp))


def build(comps=COMPONENTS, eos="PR", alpha="boston-mathias",
          he_kij=None):
    """
    The LBNH property package.

    Water is left out on purpose: it is removed in pretreatment ahead of
    the cold section and appears in the cold box only as a trace, where
    a cubic equation of state would describe it badly anyway.
    """
    from thermo.interaction_parameters import IPDB
    cas, Tc, Pc, om, MW, cp = _props(comps)
    k = np.array(IPDB.get_ip_asymmetric_matrix("ChemSep PR", cas, "kij"),
                 float)
    k = 0.5 * (k + k.T)
    if he_kij and "He" in comps:
        ih = comps.index("He")
        for other, v in he_kij.items():
            if other in comps:
                j = comps.index(other)
                k[ih, j] = k[j, ih] = v
    return Cubic(comps, Tc, Pc, om, MW, cp, kij=k, eos=eos, alpha=alpha,
                 kappa_override={"He": 0.0} if "He" in comps else None)
