"""
Nitrogen rejection column.

Getting this column to converge took four things together, each fixing a
distinct failure:

  1. solve only in the components actually fed - components with zero
     feed make the Newton matrix singular;
  2. scale each component's balances by its own feed - otherwise
     helium's small flows are swamped by methane's and helium quietly
     disappears from the material balance;
  3. start the vapour at equilibrium with the liquid - an initial guess
     with identical phases sits exactly on the trivial solution, which
     is itself a solution of the equilibrium equations;
  4. grow the column in small steps from a short one that converges.

Reflux enters as a liquid feed from an external partial condenser, so the
top stage is held inside the two-phase region, as in the ASU work.
"""
import numpy as np
from props import subpackage, KMOLH
from column import Column
from flash import flash_TP, dew_T, bubble_T

NRU_COMPS = ["He", "CH4", "N2"]


def solve_nru(feed, reflux, P, T_feed, vf_feed, T_reflux, n_stages=20,
              feed_stage=8, boilup=1.2, warm=None):
    pk = subpackage(NRU_COMPS)
    TbN, *_ = bubble_T(pk, P, np.array([0, 0, 1.0]))
    TbC, *_ = bubble_T(pk, P, np.array([0, 1.0, 0]))

    def make(N):
        c = Column(pk, N, P, P + 0.2e5)
        c.scale_mode = "component"
        c.guess = dict(T_top=TbN, T_bot=TbC, x_top=[0.002, 0.05, 0.948],
                       x_bot=[1e-7, 0.985, 0.015])
        c.add_feed(1, reflux, T=T_reflux, P=P, vapour_frac=0.0)
        c.add_feed(max(2, round(feed_stage * N / n_stages)), feed,
                   T=T_feed, P=P, vapour_frac=vf_feed)
        c.set_boilup_ratio(boilup)
        return c

    if warm is not None:
        c = make(n_stages)
        c.solve(x0=warm, max_iter=120)
        if c.converged:
            return c

    x, prev = None, None
    ladder = [n for n in (6, 8, 10, 12, 14, 16, 18) if n < n_stages] + [n_stages]
    for N in ladder:
        c = make(N)
        x0 = None if x is None else c.stretch_profile(x, prev)
        c.solve(x0=x0, max_iter=120)
        if not c.converged:
            c.solve_lm(x0=x0, max_nfev=200)
            c.solve(x0=c.x, max_iter=80)
        if not c.converged:
            break
        x, prev = c.x.copy(), N
    return c


def nru_with_condenser(feed, P, T_feed, vf_feed, T_cond, boilup,
                       n_stages=20, feed_stage=8, tol=1e-6, max_outer=60,
                       damping=0.6, state=None):
    """
    Column plus external partial condenser, with the reflux loop closed.

    The overhead vapour is partially condensed at T_cond; the liquid is
    the reflux, the vapour goes on to helium recovery. Reflux flow and
    composition are therefore results, set by the condenser temperature.
    Returns (column, reflux flows, condenser vapour flows, condenser duty
    in W).
    """
    pk = subpackage(NRU_COMPS)
    TbN, *_ = bubble_T(pk, P, np.array([0, 0, 1.0]))
    reflux = (state["reflux"] if state else
              np.array([1e-6, 0.005, 0.995]) * 0.8 * feed.sum())
    warm = state["x"] if state else None
    from units import Stream, flash_drum
    for it in range(max_outer):
        c = solve_nru(feed, reflux, P, T_feed, vf_feed, T_cond,
                      n_stages=n_stages, feed_stage=feed_stage,
                      boilup=boilup, warm=warm)
        if not c.converged:
            return None
        warm = c.x.copy()
        r = c.results()
        ov = Stream(pk, r["v"][0], r["T"][0], P, "NRU_OV", forced_beta=1.0)
        vap, liq, Qc = flash_drum(ov, T=T_cond, P=P)
        new = damping * liq.n + (1 - damping) * reflux
        err = np.max(np.abs(new - reflux)) / feed.sum()
        reflux = np.maximum(new, 1e-12)
        if err < tol:
            return dict(col=c, reflux=reflux, vapour=vap.n, T_cond=T_cond,
                        Q_cond=Qc, reb=r["duties"][-1], outer=it + 1,
                        x=c.x.copy())
    return None
