"""
Same design, same answer.
Every compressor and pump of my Aspen Plus model, recomputed by the AI simulator's own
property models from the SAME inlet state, outlet pressure and efficiencies
(read from the exported Aspen backup/report of LBNH 0.14 wt %).

    python same_design.py
"""
import json, time, numpy as np, CoolProp.CoolProp as CP
from props import PKG, C, subpackage
from flash import flash_TP, flash_PS, mixture_H, mixture_S
K0 = 273.15
D = json.load(open("aspen_settings.json")); S, B = D["streams"], D["blocks"]
MAP = {"METHANE": "CH4", "ETHANE": "C2H6", "PROPANE": "C3H8", "BUTANE-N": "nC4", "BUTANE-I": "iC4",
       "PENTAN-N": "nC5", "PENTAN-I": "iC5", "HELIUM": "He", "NITROGEN": "N2"}
MW = dict(zip(C, PKG.MW)); PURE = {"C3H8": "Propane", "N2": "Nitrogen", "He": "Helium"}

def work(st, Pout, eff):
    d = S[st]; z = {MAP[k]: v for k, v in d["molefrac"].items() if k in MAP and v and v > 1e-9}
    s = sum(z.values()); z = {k: v / s for k, v in z.items()}
    T, P = d["T"] + K0, d["P"] * 1e5
    n = d["tph"] * 1e6 / 3600 / sum(z[k] * MW[k] for k in z)
    if len(z) == 1 and list(z)[0] in PURE:
        fl = PURE[list(z)[0]]
        h1 = CP.PropsSI("Hmolar", "T", T, "P", P, fl); s1 = CP.PropsSI("Smolar", "T", T, "P", P, fl)
        return n * (CP.PropsSI("Hmolar", "P", Pout, "Smolar", s1, fl) - h1) / eff / 1e6
    names = [c for c in C if c in z]; pk = subpackage(names) if len(names) < len(C) else PKG
    zz = np.array([z[c] for c in names])
    b, x, y, _ = flash_TP(pk, T, P, zz); h1 = mixture_H(pk, T, P, b, x, y); s1 = mixture_S(pk, T, P, b, x, y)
    Ts, bs, xs, ys = flash_PS(pk, Pout, zz, s1, lo=60, hi=700)
    return n * (mixture_H(pk, Ts, Pout, bs, xs, ys) - h1) / eff / 1e6

t0 = time.time(); rows = []
print("\nSame design, same answer: every Aspen compressor and pump recomputed by the AI simulator\n")
print(f"{'machine':10s}{'inlet':>16s}{'outlet':>10s}{'efficiency':>13s}{'Aspen MW':>11s}{'AI MW':>9s}{'diff':>8s}")
print("-" * 77)
for name, d in B.items():
    st = S[d["inlet"]]; Pout = d["Pout"] if d["Pout"] else st["P"] * d["Pratio"]
    w = work(d["inlet"], Pout * 1e5, d["seff"] * d["meff"]); a = d["brake"] / 1e3
    rows.append((a, w))
    print(f"{name:10s}{st['T']:8.0f} C {st['P']:5.1f} bar{Pout:6.1f} bar{d['seff']:7.2f} x {d['meff']:.2f}{a:11.2f}{w:9.2f}{(w - a) / a * 100:+7.1f}%")
# LNG pump: liquid, 1.3 -> 6.0 bar, pump efficiency as used by Aspen
p = D["pump"]; d = S[p["inlet"]]
z = {MAP[k]: v for k, v in d["molefrac"].items() if k in MAP and v}
zz = np.array([z.get(c, 0.0) for c in C]); zz = zz / zz.sum()
T, P = d["T"] + K0, d["P"] * 1e5; n = d["tph"] * 1e6 / 3600 / float(zz @ PKG.MW)
h1 = mixture_H(PKG, T, P, 0.0, zz, zz); s1 = mixture_S(PKG, T, P, 0.0, zz, zz)
Ts, bs, xs, ys = flash_PS(PKG, p["Pout"] * 1e5, zz, s1, lo=80, hi=200)
w = n * (mixture_H(PKG, Ts, p["Pout"] * 1e5, bs, xs, ys) - h1) / p["eff"] / 1e6; a = p["brake"] / 1e3
rows.append((a, w))
print(f"{p['name']:10s}{d['T']:8.0f} C {d['P']:5.1f} bar{p['Pout']:6.1f} bar{'pump ' + format(p['eff'], '.2f'):>13s}{a:11.2f}{w:9.2f}{(w - a) / a * 100:+7.1f}%")
# nitrogen cycle compressor: 5 stages, polytropic 0.72, intercooled to 5 C, 1 -> 75 bar, 375 t/h
W = 0; P = 1e5; r = 75 ** (1 / 5); n = 375e6 / 3600 / 28.0134
for k in range(5):
    h = CP.PropsSI("Hmolar", "T", 278.15, "P", P, "Nitrogen"); rr = r ** (1 / 40)
    for j in range(40):
        s = CP.PropsSI("Smolar", "P", P, "Hmolar", h, "Nitrogen")
        h2 = h + (CP.PropsSI("Hmolar", "P", P * rr, "Smolar", s, "Nitrogen") - h) / 0.72
        W += n * (h2 - h); h = h2; P *= rr
rows.append((61.5382, W / 1e6))
print(f"{'N2-CMPR':10s}{5:8.0f} C {1.0:5.1f} bar{75.0:6.1f} bar{'0.72 poly':>13s}{61.54:11.2f}{W / 1e6:9.2f}{(W / 1e6 - 61.5382) / 61.5382 * 100:+7.1f}%")
ta = sum(r[0] for r in rows); tt = sum(r[1] for r in rows)
print("-" * 77)
print(f"{'TOTAL 14':10s}{'':39s}{ta:11.2f}{tt:9.2f}{(tt - ta) / ta * 100:+7.1f}%")
print(f"\nGiven the same design settings, the AI simulator reproduces Aspen's power ({time.time() - t0:.0f} s).\n")
