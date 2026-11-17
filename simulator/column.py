"""
Rigorous distillation by simultaneous correction (Naphtali-Sandholm)
====================================================================

pyflowsheet has only sum-rates absorbers, which will not converge a
column at high reflux. This solves all stages at once by Newton's
method, which is what a 69-stage double column and a 150-stage argon
column need.

Formulation
-----------
Unknowns per stage j:  the component vapour flows v_ij, the component
liquid flows l_ij, and the stage temperature T_j.  For C components and
N stages that is N(2C+1) unknowns, matched by N(2C+1) equations:

  M_ij   l_{i,j-1} + v_{i,j+1} + f_ij
           - (1 + U_j/L_j) l_ij - (1 + W_j/V_j) v_ij  = 0
  E_ij   K_ij (l_ij / L_j) V_j - v_ij                 = 0
  H_j    energy balance over the stage                = 0

Stage 1 is the top. Liquid entering stage 1 and vapour entering stage N
from outside are zero unless supplied as feeds.

Why the Jacobian is built the way it is
---------------------------------------
Stage j couples only to j-1 and j+1, so the Jacobian is block
tridiagonal. Rather than deriving every analytic derivative - easy to
get subtly wrong, and a wrong derivative shows up as slow convergence
rather than an error - the Jacobian is built by finite differences with
graph colouring: stages 1, 4, 7, ... are perturbed simultaneously
because they cannot influence one another. That costs 3(2C+1) residual
evaluations per iteration regardless of column height, which for three
components is 21, and each residual evaluation is a fraction of a
millisecond.
"""

import numpy as np
from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import splu


class Stage:
    """One equilibrium stage."""

    __slots__ = ("feed_l", "feed_v", "H_feed", "U", "W", "Q", "P")

    def __init__(self, nc):
        self.feed_l = np.zeros(nc)   # liquid feed component flows, mol/s
        self.feed_v = np.zeros(nc)   # vapour feed component flows, mol/s
        self.H_feed = 0.0            # enthalpy carried in by the feed, W
        self.U = 0.0                 # liquid side draw, mol/s
        self.W = 0.0                 # vapour side draw, mol/s
        self.Q = 0.0                 # heat added, W
        self.P = 1e5                 # stage pressure, Pa


class Column:
    """
    A rigorous column.

    Parameters
    ----------
    pkg       property package exposing K(T,P,x,y) and H(T,P,z,vapour)
    n_stages  number of equilibrium stages, stage 1 at the top
    P_top,
    P_bottom  linear pressure profile between them, Pa
    """

    def __init__(self, pkg, n_stages, P_top, P_bottom=None):
        self.pkg = pkg
        self.N = int(n_stages)
        self.nc = pkg.nc
        P_bottom = P_top if P_bottom is None else P_bottom
        self.P = np.linspace(P_top, P_bottom, self.N)
        self.stages = [Stage(self.nc) for _ in range(self.N)]
        for j, s in enumerate(self.stages):
            s.P = self.P[j]

        # Stages whose energy-balance row is replaced by a specification.
        # The duty at such a stage is then whatever the energy balance
        # implies, computed after the solve rather than solved for.
        self.spec_rows = {}      # stage index -> callable -> residual
        # spec row index -> the stage its residual actually depends on.
        # Colouring assumes a row responds to at most one perturbed stage
        # per pass; a non-local row breaks that, so its derivatives are
        # taken separately.
        self.nonlocal_rows = {}

    # -----------------------------------------------------------------
    # specification
    # -----------------------------------------------------------------

    @property
    def active(self):
        """
        Which components are actually fed to this column.

        A component with no feed anywhere cannot appear in any product,
        and trying to solve its equilibrium is meaningless: both its
        flows collapse toward zero, their ratio is undetermined, and in
        logarithmic form the residual is stuck at ln(K) forever. Such
        components are parked at a fixed negligible value and their two
        rows replaced by that trivial statement, which keeps the system
        square and well conditioned.
        """
        f = sum((s.feed_l + s.feed_v) for s in self.stages)
        return np.asarray(f, float) > 0.0

    @property
    def component_scale(self):
        """
        Per-component reference flow for row scaling.

        Scaling every row by the TOTAL feed leaves a trace component's
        equations with residuals around 1e-13, which is below the noise
        floor of the finite-difference Jacobian - so Newton simply
        cannot see them, and the solve stalls exactly when a product
        approaches high purity. Scaling each component by its own feed
        rate keeps all of them visible.
        """
        f = sum((s.feed_l + s.feed_v) for s in self.stages)
        f = np.asarray(f, float)
        # a component with no feed still needs a finite scale
        return np.maximum(f, max(f.sum(), 1e-9) * 1e-6)

    @property
    def feed_scale(self):
        """Total feed rate, used to normalise specification residuals.

        Normalising by the live V or L is unstable: early in the solve
        those can approach zero and the residual blows up, which is what
        wrecks the Newton direction."""
        f = sum((s.feed_l + s.feed_v).sum() for s in self.stages)
        return max(f, 1e-9)

    def add_feed(self, stage, flows, T=None, P=None, vapour_frac=None,
                 H=None):
        """
        Feed component molar flows (mol/s) into a stage, 1-indexed.

        The enthalpy is evaluated at the FEED's own temperature and
        pressure, not the stage's. A feed arrives with a state it
        already has; what happens to it on the stage is what the stage
        equations are there to determine. Evaluating it at the stage
        pressure instead silently changes how much vapour it is assumed
        to carry - the nitrogen feed entering the low pressure column
        arrives at 1.3 bar with 17% vapour but is fully vapour at the
        1.1 bar stage pressure, a 40% error in its enthalpy.

        The two phases also get their own compositions from a real
        flash, rather than both being given the overall composition.
        """
        s = self.stages[stage - 1]
        flows = np.asarray(flows, float)
        tot = flows.sum()
        if tot <= 0:
            return
        z = flows / tot

        P_feed = s.P if P is None else P

        if vapour_frac is None and T is not None:
            from flash import flash_TP
            vapour_frac, xf, yf, _ = flash_TP(self.pkg, T, P_feed, z)
        else:
            vapour_frac = 0.0 if vapour_frac is None else float(vapour_frac)
            xf = yf = z
            if 1e-9 < vapour_frac < 1.0 - 1e-9 and T is not None:
                from flash import flash_TP
                _, xf, yf, _ = flash_TP(self.pkg, T, P_feed, z)

        s.feed_v += flows * vapour_frac
        s.feed_l += flows * (1.0 - vapour_frac)

        if H is not None:
            s.H_feed += H
        elif T is not None:
            Ta, Pa = np.array([T]), np.array([P_feed])
            h = 0.0
            if vapour_frac > 1e-12:
                h += vapour_frac * self.pkg.H(Ta, Pa, yf[None, :],
                                              vapour=True)[0]
            if vapour_frac < 1.0 - 1e-12:
                h += (1.0 - vapour_frac) * self.pkg.H(Ta, Pa, xf[None, :],
                                                      vapour=False)[0]
            s.H_feed += tot * h

    def add_heat(self, stage, Q):
        """Heat added to a stage, W. Negative removes heat."""
        self.stages[stage - 1].Q += Q

    # ----- specifications -------------------------------------------
    #
    # Each replaces one stage's energy balance. All of them involve only
    # that stage's own unknowns, so the Jacobian stays block tridiagonal.

    def set_reflux_ratio(self, value, stage=1):
        """
        Reflux ratio L/D at a condenser stage.

        Here the vapour leaving stage 1 is the overhead product and the
        liquid leaving stage 1 is the internal reflux, so the ratio is a
        plain relation between two unknowns on the same stage.
        """
        j = stage - 1
        self._reflux_hint = float(value)
        self.spec_rows[j] = (
            lambda v, l, T, V, L, j=j, r=float(value), s=self:
                (L[j] - r * V[j]) / s.feed_scale)

    def set_boilup_ratio(self, value, stage=None):
        """Boilup ratio V/B at a reboiler stage."""
        j = (self.N - 1) if stage is None else stage - 1
        self._boilup_hint = float(value)
        self.spec_rows[j] = (
            lambda v, l, T, V, L, j=j, r=float(value), s=self:
                (V[j] - r * L[j]) / s.feed_scale)

    def set_stage_temperature(self, stage, T_value, duty_stage=None):
        """Drive a stage to a temperature, freeing a duty to achieve it."""
        j = stage - 1
        k = j if duty_stage is None else duty_stage - 1
        self.spec_rows[k] = (
            lambda v, l, T, V, L, j=j, t=float(T_value): (T[j] - t) / 100.0)

    def set_vapour_product(self, flow, stage=1, duty_stage=None):
        """
        Fix the overhead vapour rate.

        The constraint is on stage `stage`, but it is achieved by the
        duty at `duty_stage` - the reboiler, if the condenser is already
        carrying a reflux specification. That makes the row non-local:
        it reaches a stage other than its own, so the Jacobian cannot
        assume it stays inside the tridiagonal band.
        """
        j = stage - 1
        k = j if duty_stage is None else duty_stage - 1
        self.spec_rows[k] = (
            lambda v, l, T, V, L, j=j, f=float(flow), s=self:
                (V[j] - f) / s.feed_scale)
        if k != j:
            self.nonlocal_rows[k] = j

    def set_distillate_to_feed(self, ratio, duty_stage=None):
        """D:F specification, as Aspen states it for the argon column."""
        self.set_vapour_product(ratio * self.feed_scale, stage=1,
                                duty_stage=duty_stage or self.N)

    def set_liquid_product(self, flow, stage=None, duty_stage=None):
        """Fix the bottoms liquid rate."""
        j = (self.N - 1) if stage is None else stage - 1
        k = j if duty_stage is None else duty_stage - 1
        self.spec_rows[k] = (
            lambda v, l, T, V, L, j=j, f=float(flow), s=self:
                (L[j] - f) / s.feed_scale)
        if k != j:
            self.nonlocal_rows[k] = j

    def add_liquid_draw(self, stage, flow):
        self.stages[stage - 1].U += flow

    def add_vapour_draw(self, stage, flow):
        self.stages[stage - 1].W += flow

    # -----------------------------------------------------------------
    # unknown vector packing
    #
    # The component flows are carried as their LOGARITHMS. In a column
    # at high reflux the trace component at the product end falls to
    # 1e-7 of the bulk flow while the bulk itself is order 1e2, and
    # Newton working on the raw flows cannot resolve the small end: an
    # absolute step that is sensible for the bulk obliterates the trace.
    # Working in logarithms gives every flow the same RELATIVE
    # resolution, and has the side benefit that a flow can never go
    # negative, so no clipping is needed to keep the thermodynamics
    # defined.
    # -----------------------------------------------------------------

    @property
    def nvar(self):
        return self.N * (2 * self.nc + 1)

    nstage_var = nvar

    LOG_MIN, LOG_MAX = -120.0, 25.0
    # Where an unfed component is parked. This must be exactly where the
    # initial guess puts it, otherwise those rows start with a large
    # residual, dominate the first Newton step, and the line search
    # accepts a direction that fixes an absent component while wrecking
    # the real ones.
    INERT_LOG = -80.0

    def _pack(self, v, l, T):
        nc, N = self.nc, self.N
        blk = 2 * nc + 1
        x = np.empty(self.nvar)
        floor = np.exp(self.INERT_LOG)
        lv = np.log(np.maximum(v, floor))
        ll = np.log(np.maximum(l, floor))
        for j in range(N):
            x[j*blk:j*blk+nc] = lv[j]
            x[j*blk+nc:j*blk+2*nc] = ll[j]
            x[j*blk+2*nc] = T[j]
        return x

    def _unpack(self, x):
        nc, N = self.nc, self.N
        blk = 2 * nc + 1
        xr = x.reshape(N, blk)
        lv = np.clip(xr[:, :nc], self.LOG_MIN, self.LOG_MAX)
        ll = np.clip(xr[:, nc:2*nc], self.LOG_MIN, self.LOG_MAX)
        return np.exp(lv), np.exp(ll), xr[:, 2*nc].copy()

    # -----------------------------------------------------------------
    # residuals
    # -----------------------------------------------------------------

    def residuals(self, x):
        nc, N = self.nc, self.N
        v, l, T = self._unpack(x)

        # flows are exponentials, so already strictly positive
        T = np.clip(T, 40.0, 400.0)

        V = v.sum(axis=1)
        L = l.sum(axis=1)
        y = v / V[:, None]
        xm = l / L[:, None]

        K = self.pkg.K(T, self.P, xm, y)
        HV = self.pkg.H(T, self.P, y, vapour=True)
        HL = self.pkg.H(T, self.P, xm, vapour=False)

        U = np.array([s.U for s in self.stages])
        W = np.array([s.W for s in self.stages])
        Q = np.array([s.Q for s in self.stages], float)
        FL = np.array([s.feed_l for s in self.stages])
        FV = np.array([s.feed_v for s in self.stages])
        HF = np.array([s.H_feed for s in self.stages])

        # draw ratios
        ru = 1.0 + U / L
        rw = 1.0 + W / V

        R = np.empty((N, 2 * nc + 1))

        F = self.feed_scale

        # --- M: component material balance ---
        l_above = np.zeros_like(l); l_above[1:] = l[:-1]
        v_below = np.zeros_like(v); v_below[:-1] = v[1:]
        R[:, :nc] = (l_above + v_below + FL + FV
                     - ru[:, None] * l - rw[:, None] * v)

        # --- E: equilibrium ---
        R[:, nc:2*nc] = K * xm * V[:, None] - v

        # --- H: energy balance ---
        L_above = np.zeros(N); L_above[1:] = L[:-1]
        HL_above = np.zeros(N); HL_above[1:] = HL[:-1]
        V_below = np.zeros(N); V_below[:-1] = V[1:]
        HV_below = np.zeros(N); HV_below[:-1] = HV[1:]

        R[:, 2*nc] = (L_above * HL_above + V_below * HV_below + HF + Q
                      - ru * L * HL - rw * V * HV)

        # scale the energy row so Newton sees all rows on a similar
        # footing; enthalpies are ~1e4 J/mol against flows of ~1e2 mol/s
        scale = np.maximum(np.abs(ru * L * HL) + np.abs(rw * V * HV), 1.0)
        R[:, 2*nc] /= scale

        # Scale the material and equilibrium rows by the feed rate.
        # Without this a trace component contributes a residual of order
        # 1e-20 while a bulk component contributes order 1e2, and Newton
        # produces directions the line search has to cut to nothing.
        #
        # This single global scale is a compromise, and it is the reason
        # the argon column stalls near RR=24: at that point the top-stage
        # oxygen residual falls below the finite-difference noise floor.
        # Per-component scaling and local row equilibration were both
        # tried and were worse - the first lets a 3 ppm nitrogen
        # impurity dominate the system, the second makes the residual
        # non-smooth where the largest term switches.
        if getattr(self, "scale_mode", "total") == "component":
            # Each component's balances scaled by its own feed. Only safe
            # when every component is a significant fraction of the feed:
            # with a trace impurity the trace rows dominate instead.
            cs = sum((s.feed_l + s.feed_v) for s in self.stages)
            cs = np.maximum(np.asarray(cs, float), F * 1e-6)
            R[:, :nc] /= cs[None, :]
            R[:, nc:2*nc] /= cs[None, :]
        else:
            R[:, :2*nc] /= F

        # a specified stage trades its energy balance for its spec
        for j, f in self.spec_rows.items():
            R[j, 2*nc] = f(v, l, T, V, L)

        return R.reshape(-1)

    # -----------------------------------------------------------------
    # Jacobian by coloured finite differences
    # -----------------------------------------------------------------

    def jacobian(self, x, R0=None, eps=1e-7):
        """
        Block-tridiagonal Jacobian by coloured finite differences.

        Stage j couples only to j-1 and j+1, so stages three apart can be
        perturbed in the same pass. That is 3(2C+1) residual evaluations
        per iteration whatever the column height - 21 for three
        components, whether the column has 20 stages or 150.
        """
        nc, N = self.nc, self.N
        blk = 2 * nc + 1
        n = self.nvar
        if R0 is None:
            R0 = self.residuals(x)

        J = lil_matrix((n, n))
        for colour in range(3):
            stages = list(range(colour, N, 3))
            if not stages:
                continue
            for k in range(blk):
                dx = np.zeros(n)
                idx = []
                for j in stages:
                    i = j * blk + k
                    # log variables take an absolute step (which is a
                    # relative step in the flow); temperature takes one
                    # scaled to its own magnitude
                    h = 1e-6 if k < 2 * nc else eps * max(abs(x[i]), 1.0)
                    dx[i] = h
                    idx.append((j, i, h))
                d = self.residuals(x + dx) - R0
                for j, i, h in idx:
                    lo = max(0, (j - 1) * blk)
                    hi = min(n, (j + 2) * blk)
                    J[lo:hi, i] = (d[lo:hi] / h).reshape(-1, 1)

        # Non-local specification rows, one column at a time.
        # Perturbing a whole colour at once would fold several stages'
        # effects into a single difference, and assigning that sum to
        # each stage's column multiply-counts the derivative - which
        # leaves the specification uncorrectable and the solve stuck.
        for row_stage, dep_stage in self.nonlocal_rows.items():
            row = row_stage * blk + 2 * nc
            for k in range(blk):
                i = dep_stage * blk + k
                h = 1e-6 if k < 2 * nc else eps * max(abs(x[i]), 1.0)
                dx = np.zeros(n); dx[i] = h
                J[row, i] = (self.residuals(x + dx)[row] - R0[row]) / h

        return csc_matrix(J)

    # -----------------------------------------------------------------
    # solve
    # -----------------------------------------------------------------

    def solve(self, x0=None, tol=1e-10, max_iter=150, verbose=False):
        x = self.initial_guess() if x0 is None else x0.copy()
        R = self.residuals(x)
        norm = np.linalg.norm(R)
        history = [norm]

        for it in range(max_iter):
            if norm < tol:
                break
            J = self.jacobian(x, R)
            try:
                dx = splu(J).solve(-R)
            except Exception:
                # singular Jacobian: nudge the diagonal and retry once
                J = J + 1e-10 * csc_matrix(np.eye(self.nvar))
                try:
                    dx = splu(J).solve(-R)
                except Exception:
                    break

            # Damped line search, accepting a step that merely does not
            # make things much worse. Strict descent stalls in the
            # shallow valleys these systems have; allowing a small
            # transient increase lets Newton through to the basin where
            # it converges quadratically.
            step = 1.0
            best = (None, None, np.inf)
            for _ in range(30):
                xt = self._clip(x + step * dx)
                Rt = self.residuals(xt)
                nt = np.linalg.norm(Rt)
                if nt < best[2]:
                    best = (xt, Rt, nt)
                if nt < norm:
                    break
                step *= 0.5
            else:
                xt, Rt, nt = best
                if xt is None or nt > norm * 1.5:
                    break

            x, R, norm = xt, Rt, nt
            history.append(norm)
            if verbose:
                print(f"  iter {it+1:>3}  |R| = {norm:.4e}  step = {step:.3f}")

        self.x = x
        self.converged = norm < tol * 1e3
        self.residual_norm = norm
        self.history = history
        self.iterations = len(history) - 1
        return self.results()

    def _clip(self, x):
        nc, N = self.nc, self.N
        blk = 2 * nc + 1
        x = x.copy()
        xr = x.reshape(N, blk)
        xr[:, :2*nc] = np.clip(xr[:, :2*nc], self.LOG_MIN, self.LOG_MAX)
        xr[:, 2*nc] = np.clip(xr[:, 2*nc], 50.0, 400.0)
        return xr.reshape(-1)

    def initial_guess(self):
        """
        A starting point that already knows roughly what the column does.

        Initialising every stage at feed composition is the usual cause
        of a stalled Newton solve: for a column that actually separates,
        that point is far from the answer and the first Jacobian points
        somewhere useless. Instead:

          * split the feed into a light-enriched top and a heavy-enriched
            bottom using K-values at the feed bubble point,
          * interpolate composition down the column,
          * take each stage temperature from its own bubble point,
          * set flows by constant molar overflow, arranged so the
            reflux and boilup specifications are already satisfied.
        """
        nc, N = self.nc, self.N
        total = sum(s.feed_l + s.feed_v for s in self.stages)
        F = max(total.sum(), 1e-9)
        z = total / F

        # --- composition profile ---------------------------------------
        Tf = self._bubble_T(z, float(self.P.mean()))
        Ta, Pa = np.array([Tf]), np.array([float(self.P.mean())])
        zz = z[None, :]
        y = zz.copy()
        for _ in range(60):
            K = self.pkg.K(Ta, Pa, zz, y)
            yn = zz * K
            yn /= yn.sum(axis=1, keepdims=True)
            if np.max(np.abs(yn - y)) < 1e-12:
                y = yn; break
            y = yn
        K = self.pkg.K(Ta, Pa, zz, y)[0]
        K = np.clip(K, 1e-6, 1e6)

        x_top = z * K ** 1.5
        x_top /= x_top.sum()
        x_bot = z / K ** 1.5
        x_bot /= x_bot.sum()

        w = np.linspace(0.0, 1.0, N)[:, None]
        xprof = (1.0 - w) * x_top[None, :] + w * x_bot[None, :]
        xprof /= xprof.sum(axis=1, keepdims=True)

        # --- temperature from each stage's own bubble point -------------
        # A caller can supply the end points instead. With a
        # non-condensable in the feed - helium here - bubble points are
        # undefined, the search falls back to an arbitrary midpoint, and
        # the solver starts from an upside-down profile it cannot recover
        # from.
        g = getattr(self, "guess", None)
        if g:
            if g.get("x_top") is not None:
                x_top = np.asarray(g["x_top"], float); x_top /= x_top.sum()
            if g.get("x_bot") is not None:
                x_bot = np.asarray(g["x_bot"], float); x_bot /= x_bot.sum()
            xprof = (1.0 - w) * x_top[None, :] + w * x_bot[None, :]
            xprof /= xprof.sum(axis=1, keepdims=True)
            T = np.linspace(g["T_top"], g["T_bot"], N)
        else:
            T = np.array([self._bubble_T(xprof[j], self.P[j]) for j in range(N)])
            # keep it monotone; a ragged profile makes the first step erratic
            T = np.maximum.accumulate(T)

        # --- flows by constant molar overflow, honouring the specs ------
        V = np.full(N, F)
        L = np.full(N, F)
        r = getattr(self, "_reflux_hint", None)
        if r is not None:
            V[0] = F / (1.0 + r)
            L[0] = r * F / (1.0 + r)
        b = getattr(self, "_boilup_hint", None)
        if b is not None:
            L[-1] = F / (1.0 + b)
            V[-1] = b * F / (1.0 + b)

        # Vapour composition from equilibrium with the guessed liquid, not
        # a copy of it. Giving both phases the same composition starts
        # Newton exactly on the trivial solution (y = x, every K = 1),
        # which is itself a solution of the equilibrium equations - and
        # with a non-condensable in the column the solver stays there.
        yprof = xprof.copy()
        for _ in range(8):
            K = self.pkg.K(T, self.P, xprof, yprof)
            yn = np.clip(K, 1e-12, 1e12) * xprof
            yprof = yn / yn.sum(axis=1, keepdims=True)
        v = V[:, None] * yprof
        l = L[:, None] * xprof
        return self._pack(v, l, T)

    def _bubble_T(self, z, P, lo=60.0, hi=None):
        from scipy.optimize import brentq
        pkg = self.pkg
        hi = min(float(pkg.Tc.max()) * 0.99, 400.0) if hi is None else hi
        zz = np.atleast_2d(z)
        Pa = np.array([P])

        def f(T):
            Ta = np.array([T])
            y = zz.copy()
            for _ in range(80):
                phil, Zl = pkg.phi(Ta, Pa, zz, vapour=False)
                phiv, Zv = pkg.phi(Ta, Pa, y, vapour=True)
                K = phil / phiv
                yn = zz * K
                yn /= yn.sum(axis=1, keepdims=True)
                if np.max(np.abs(yn - y)) < 1e-12:
                    y = yn; break
                y = yn
            phil, Zl = pkg.phi(Ta, Pa, zz, vapour=False)
            phiv, Zv = pkg.phi(Ta, Pa, y, vapour=True)
            if abs(Zv[0] - Zl[0]) < 1e-8:
                return -1e3 if Zl[0] < 0.35 else 1e3
            return float((zz * (phil / phiv)).sum() - 1.0)

        try:
            if f(lo) * f(hi) > 0:
                return 0.5 * (lo + hi)
            return brentq(f, lo, hi, xtol=1e-8)
        except Exception:
            return 0.5 * (lo + hi)

    # -----------------------------------------------------------------

    def results(self):
        v, l, T = self._unpack(self.x)
        V, L = v.sum(axis=1), l.sum(axis=1)
        return dict(
            duties=self.implied_duties(),
            converged=self.converged,
            iterations=self.iterations,
            residual=self.residual_norm,
            T=T, P=self.P.copy(),
            V=V, L=L,
            y=v / V[:, None], x=l / L[:, None],
            v=v, l=l,
        )

    # -----------------------------------------------------------------

    def implied_duties(self):
        """
        Heat duty on every stage, from the converged profiles.

        On a specified stage the energy balance was replaced by the
        specification, so its duty is whatever the balance implies -
        that is the condenser or reboiler duty, and it is recovered here
        rather than solved for.
        """
        nc, N = self.nc, self.N
        v, l, T = self._unpack(self.x)
        V, L = v.sum(axis=1), l.sum(axis=1)
        y, xm = v / V[:, None], l / L[:, None]

        HV = self.pkg.H(T, self.P, y, vapour=True)
        HL = self.pkg.H(T, self.P, xm, vapour=False)
        U = np.array([s.U for s in self.stages])
        W = np.array([s.W for s in self.stages])
        HF = np.array([s.H_feed for s in self.stages])
        ru, rw = 1.0 + U / L, 1.0 + W / V

        L_above = np.zeros(N); L_above[1:] = L[:-1]
        HL_above = np.zeros(N); HL_above[1:] = HL[:-1]
        V_below = np.zeros(N); V_below[:-1] = V[1:]
        HV_below = np.zeros(N); HV_below[:-1] = HV[1:]

        return (ru * L * HL + rw * V * HV
                - L_above * HL_above - V_below * HV_below - HF)

    # -----------------------------------------------------------------

    def solve_continuation(self, spec_setter, targets, tol=1e-10,
                           max_iter=150, verbose=False):
        """
        Reach a hard specification by walking up to it.

        A column at very high reflux is a poor Newton problem from a
        cold start: the product end spans seven decades of composition
        and the initial guess is nowhere near it. The same column at
        moderate reflux is easy. So the specification is approached in
        steps, each solve starting from the last converged profile,
        which is the standard way these columns are made to converge.

        `spec_setter(value)` applies one value of the specification;
        `targets` is the ladder of values ending at the real one.
        """
        x = None
        last = None
        for t in targets:
            spec_setter(t)
            res = self.solve(x0=x, tol=tol, max_iter=max_iter)
            if verbose:
                print(f"    spec={t:<8g} iters={self.iterations:>4}"
                      f"  |R|={self.residual_norm:.3e}"
                      f"  {'ok' if self.converged else 'STALLED'}")
            if not self.converged and last is not None:
                # this rung failed; keep the last good profile
                self.x = last
                break
            x = self.x.copy()
            last = x
            res_final = res
        return self.results()

    # -----------------------------------------------------------------

    def stretch_profile(self, source, n_source):
        """
        Interpolate a converged profile from a shorter column onto this one.

        A tall column started from a generic guess is a poor Newton
        problem - the guess is simply too far from the answer, and the
        solve stalls even though the same system converges cleanly when
        it is short. A short column, though, converges easily, and its
        profile has the right SHAPE: composition and temperature vary
        smoothly from one end to the other. Stretching that shape over
        more stages gives the tall column a starting point that is
        already qualitatively correct, which is all Newton needs.
        """
        nc, N = self.nc, self.N
        blk = 2 * nc + 1
        src = source.reshape(n_source, blk)

        u_src = np.linspace(0.0, 1.0, n_source)
        u_dst = np.linspace(0.0, 1.0, N)
        out = np.empty((N, blk))
        for k in range(blk):
            out[:, k] = np.interp(u_dst, u_src, src[:, k])

        # Flows scale with column height only through the feeds, which
        # are already set on this column, so the interpolated logarithms
        # carry over directly.
        return out.reshape(-1)

    def solve_by_growth(self, ladder=(8, 16, 32), tol=1e-10,
                        max_iter=250, verbose=False):
        """
        Solve by growing the column: short first, then stretch upward.

        Feed and draw stages are placed proportionally on each short
        column so its profile has the same shape as the real one.
        """
        import copy
        x = None
        for n in list(ladder) + [self.N]:
            if n >= self.N:
                short = self
            else:
                short = Column(self.pkg, n, self.P[0], self.P[-1])
                if getattr(self, "guess", None):
                    short.guess = self.guess
                scale = (n - 1) / max(self.N - 1, 1)
                for j, s in enumerate(self.stages):
                    tgt = int(round(j * scale))
                    t = short.stages[tgt]
                    t.feed_l = t.feed_l + s.feed_l
                    t.feed_v = t.feed_v + s.feed_v
                    t.H_feed += s.H_feed
                    t.U += s.U
                    t.W += s.W
                    t.Q += s.Q
                for j, f in self.spec_rows.items():
                    tgt = int(round(j * scale))
                    if hasattr(self, "_boilup_hint") and j == self.N - 1:
                        short.set_boilup_ratio(self._boilup_hint)
                    elif hasattr(self, "_reflux_hint") and j == 0:
                        short.set_reflux_ratio(self._reflux_hint)

            x0 = None if x is None else short.stretch_profile(x, prev)
            short.solve(x0=x0, tol=tol, max_iter=max_iter)
            if verbose:
                print(f"    {n:>3} stages  iters={short.iterations:>4}"
                      f"  |R|={short.residual_norm:.2e}"
                      f"  {'ok' if short.converged else 'STALL'}")
            if not short.converged and n < self.N:
                continue
            x = short.x.copy()
            prev = n
        return self.results()

    # -----------------------------------------------------------------

    def solve_lm(self, x0=None, tol=1e-10, max_nfev=400, verbose=False):
        """
        Trust-region solve (scipy least_squares, 'trf').

        Newton with a backtracking line search can only shorten a bad
        step. Near a trivial solution or across a region where a phase
        disappears, the Newton direction itself is wrong, and shortening
        it just stalls. A trust region bends the step toward steepest
        descent when the model stops being trustworthy, which is what
        gets the helium-bearing nitrogen rejection column through.
        Uses the same coloured sparse Jacobian as the Newton solver.
        """
        from scipy.optimize import least_squares
        x = self.initial_guess() if x0 is None else np.asarray(x0, float)
        nc, N = self.nc, self.N
        blk = 2 * nc + 1
        lo = np.full(self.nvar, self.LOG_MIN)
        hi = np.full(self.nvar, self.LOG_MAX)
        lo[2 * nc::blk] = 50.0
        hi[2 * nc::blk] = 400.0
        x = np.clip(x, lo + 1e-9, hi - 1e-9)
        res = least_squares(self.residuals, x,
                            jac=lambda z: self.jacobian(z),
                            bounds=(lo, hi), method="trf",
                            x_scale="jac", ftol=1e-15, xtol=1e-15,
                            gtol=1e-15, max_nfev=max_nfev,
                            tr_solver="lsmr", verbose=2 if verbose else 0)
        self.x = res.x
        self.residual_norm = float(np.linalg.norm(res.fun))
        self.converged = self.residual_norm < tol * 1e3
        self.iterations = int(res.nfev)
        self.history = []
        return self.results()
