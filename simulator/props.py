"""The LBNH property package, built once with the fitted helium parameters."""
import json, os
import numpy as np
from eos import build, COMPONENTS
import design_basis as DB

_here = os.path.dirname(os.path.abspath(__file__))
_k = json.load(open(os.path.join(_here, "helium_kij.json")))
# He-CH4 and He-N2 fitted to GERG-2008. Heavier hydrocarbons take the
# He-CH4 value: helium is essentially insoluble in all of them, and GERG
# gives no better target that the cubic flash can reach.
HE_KIJ = {"CH4": _k["CH4"], "N2": _k["N2"], "C2H6": _k["CH4"],
          "C3H8": _k["CH4"], "iC4": _k["CH4"], "nC4": _k["CH4"],
          "iC5": _k["CH4"], "nC5": _k["CH4"]}

PKG = build(COMPONENTS, eos=DB.EOS, alpha=DB.ALPHA, he_kij=HE_KIJ)
C = COMPONENTS
IDX = {c: i for i, c in enumerate(C)}
KMOLH = 1000.0 / 3600.0          # kmol/h -> mol/s


def vec(d):
    """Component dict -> array in package order."""
    return np.array([d.get(c, 0.0) for c in C], float)


def feed_kmolh():
    from chemicals import MW, CAS_from_any
    from eos import _NAMES
    out = {}
    for c, tph in DB.FEED_TPH.items():
        if c == "H2O":
            continue                     # removed in dehydration
        out[c] = tph * 1000.0 / MW(CAS_from_any(_NAMES[c]))
    return out


_SUB = {}


def subpackage(comps):
    """
    The same thermodynamics restricted to a subset of components.

    A column fed with only some of the nine components must be solved in
    only those. Components with zero feed sit at the floor of the log
    variables on every stage; their equations reduce to 0 = 0 and their
    Jacobian columns vanish, so the Newton matrix is singular - which is
    what stalled the nitrogen rejection column.
    """
    key = tuple(comps)
    if key not in _SUB:
        _SUB[key] = build(list(comps), eos=DB.EOS, alpha=DB.ALPHA,
                          he_kij=HE_KIJ)
    return _SUB[key]
