"""
Cash-flow model on the paper's financial basis (Table 8):
TCI = CAPEX + 5 % working capital, sales growing 5 %/y, straight-line
depreciation over 20 y to a 10 % salvage value, 25 % tax, 10 % discount,
20-year life, 8322 operating hours.
"""
import numpy as np

BASIS = dict(wc=0.05, growth=0.05, life=20, salvage=0.10, tax=0.25,
             rate=0.10, hours=8322,
             price_ng_t=180.0, price_lng_t=330.0, price_he_t=100000.0,
             price_n2_t=130.0, elec_kwh=0.0775)


def cashflow(capex, sales1, opex, b=BASIS, opex_growth=0.0):
    wc = b["wc"] * capex
    tci = capex + wc
    dep = capex * (1 - b["salvage"]) / b["life"]
    cf = [-tci]
    for y in range(1, b["life"] + 1):
        s = sales1 * (1 + b["growth"]) ** (y - 1)
        o = opex * (1 + opex_growth) ** (y - 1)
        ebt = s - o - dep
        c = ebt * (1 - b["tax"]) + dep
        if y == b["life"]:
            c += b["salvage"] * capex + wc
        cf.append(c)
    cf = np.array(cf); t = np.arange(len(cf))
    npv = float(np.sum(cf / (1 + b["rate"]) ** t))
    f = lambda r: np.sum(cf / (1 + r) ** t)
    lo, hi = -0.5, 2.0
    for _ in range(200):
        m = 0.5 * (lo + hi)
        (lo, hi) = (lo, m) if f(lo) * f(m) <= 0 else (m, hi)
    cum = np.cumsum(cf)
    k = int(np.argmax(cum >= 0))
    pbp = (k - 1) + (-cum[k - 1] / cf[k]) if k > 0 else None
    pi = 1 + npv / tci
    return dict(npv=npv, irr=0.5 * (lo + hi), pbp=pbp, pi=pi, tci=tci,
                cf=cf.tolist())


if __name__ == "__main__":
    # the paper's own numbers through this engine
    r = cashflow(4499.61, 2204.92, 1398.84)
    print("paper's CAPEX/sales/OPEX through this cash-flow engine:")
    print(f"  NPV {r['npv']:,.2f} M  (paper 5,958.67)")
    print(f"  IRR {r['irr']*100:.2f} %  (paper 21.5)")
    print(f"  payback {r['pbp']:.2f} y  (paper 7.7)")
    print(f"  PI {r['pi']:.2f}  (paper 2.26)")
    print(f"  TCI {r['tci']:,.2f} M  (paper 4,724.59)")
