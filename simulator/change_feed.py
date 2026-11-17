"""
LBNH Process Twin - change the feed and solve the whole plant again.

    python change_feed.py                 asks three questions
    python change_feed.py 660 0.28 9      feed t/h, helium wt %, nitrogen wt %

The feed flow, the helium content and the nitrogen content can be changed.
The hydrocarbons keep their proportions and fill the rest of the feed.
The whole plant is then solved from scratch, with no Aspen Plus.
"""
import sys, time
import design_basis as DB

BASE_FEED = dict(DB.FEED_TPH)
F0 = sum(BASE_FEED.values())
HE0 = BASE_FEED["He"] / F0 * 100
N20 = BASE_FEED["N2"] / F0 * 100
# the design case of the presentation (python demo.py)
BASE = dict(lng=555.81, he=807.35, n2=31.51, LNGP=155.39, BOGC=10.98, NRHR=40.66, HLS=9.35)
PRICE = dict(lng=330.0, he=100000.0, n2=130.0)      # $/t, as in the paper
HOURS = 8322.0
LIMITS = dict(flow=(300.0, 900.0), he=(0.05, 0.40), n2=(4.0, 13.0))
TIME_LIMIT = 150      # seconds


def ask(text, default, lo, hi, show):
    while True:
        s = input(f"{text} [{show}] (allowed {lo:g} to {hi:g}): ").strip()
        if not s:
            return default
        try:
            v = float(s)
        except ValueError:
            print("   please type a number"); continue
        if lo <= v <= hi:
            return v
        print(f"   outside the range the model was tested for ({lo:g} to {hi:g})")


def set_feed(flow, he_wt, n2_wt):
    hc = {k: v for k, v in BASE_FEED.items() if k not in ("He", "N2")}
    s = sum(hc.values()); rest = flow * (1 - he_wt / 100 - n2_wt / 100)
    f = {k: v / s * rest for k, v in hc.items()}
    f["He"] = flow * he_wt / 100; f["N2"] = flow * n2_wt / 100
    DB.FEED_TPH.clear(); DB.FEED_TPH.update(f)


def revenue(d):
    r = {k: d[k] * PRICE[k] * HOURS / 1e6 for k in ("lng", "n2")}
    r["he"] = d["he"] / 1000 * PRICE["he"] * HOURS / 1e6
    r["total"] = r["lng"] + r["he"] + r["n2"]
    return r


def main():
    worker = "--worker" in sys.argv
    if not worker:
        print("\nLBNH Process Twin: change the feed\n")
        print(f"Design feed: {F0:.0f} t/h, helium {HE0:.2f} wt %, nitrogen {N20:.1f} wt %\n")
    if len([a for a in sys.argv[1:] if a != "--worker"]) >= 3:
        flow, he, n2 = (float(x) for x in sys.argv[1:4])
        for name, v in (("flow", flow), ("he", he), ("n2", n2)):
            lo, hi = LIMITS[name]
            if not lo <= v <= hi:
                print(f"{name} = {v:g} is outside the tested range {lo:g} to {hi:g}"); return
    else:
        print("Press Enter to keep the design value.\n")
        flow = ask("Feed flow, t/h", F0, *LIMITS["flow"], f"{F0:.0f}")
        he = ask("Helium in the feed, wt %", HE0, *LIMITS["he"], f"{HE0:.2f}")
        n2 = ask("Nitrogen in the feed, wt %", N20, *LIMITS["n2"], f"{N20:.1f}")
    if "--worker" not in sys.argv:
        # Solve in a separate process so a feed the model cannot handle stops
        # after TIME_LIMIT seconds instead of running forever.
        import subprocess
        try:
            subprocess.run([sys.executable, __file__, str(flow), str(he), str(n2), "--worker"], timeout=TIME_LIMIT)
        except subprocess.TimeoutExpired:
            print(f"\nThe model did not finish this feed in {TIME_LIMIT} s and was stopped.")
            print("This is an honest limit of the tool: there are feeds it cannot solve.\n")
        return
    set_feed(flow, he, n2)
    print(f"Solving the whole plant for {flow:.0f} t/h, helium {he:.2f} wt %, nitrogen {n2:.1f} wt % ...\n")
    from plant import run
    t0 = time.time()
    try:
        R = run(nrhr="cycle")
    except Exception as e:
        print(f"\nThe model could not solve this feed: {str(e)[:120]}")
        print("This is the honest limit of the tool: this feed is outside what it can handle.\n")
        return
    dt = time.time() - t0
    p, w = R["products"], R["power_MW"]
    new = dict(lng=p["lng_tph"], he=p["he_kgh"], n2=p["n2_kmolh"] * 28.0 / 1000, LNGP=w["c3mr"],
               BOGC=w["bog_compressor"], NRHR=w["nru_condenser"] + w["he_separation"], HLS=w["helium_liquefier"])
    for d in (BASE, new):
        d["total"] = d["LNGP"] + d["BOGC"] + d["NRHR"] + d["HLS"]
    rb, rn = revenue(BASE), revenue(new)
    rows = [("LNG, t/h", "lng", 1), ("Liquid helium, kg/h", "he", 0), ("Nitrogen, t/h", "n2", 1),
            ("Power LNG train, MW", "LNGP", 1), ("Power BOG compression, MW", "BOGC", 1),
            ("Power N2 rejection + He recovery, MW", "NRHR", 1), ("Power He liquefier, MW", "HLS", 1),
            ("Total power, MW", "total", 1)]
    print("\n" + "-" * 74)
    print(f"{'':38s}{'Design feed':>12s}{'New feed':>12s}{'change':>10s}")
    print("-" * 74)
    for label, k, d in rows:
        a, b = BASE[k], new[k]
        print(f"{label:38s}{a:12.{d}f}{b:12.{d}f}{(b - a) / a * 100:+9.1f}%")
    print("-" * 74)
    for label, k in (("Revenue from LNG, $M/year", "lng"), ("Revenue from helium, $M/year", "he"),
                     ("Revenue from nitrogen, $M/year", "n2"), ("Total revenue, $M/year", "total")):
        a, b = rb[k], rn[k]
        print(f"{label:38s}{a:12.0f}{b:12.0f}{(b - a) / a * 100:+9.1f}%")
    print(f"{'Helium share of revenue':38s}{rb['he'] / rb['total'] * 100:11.0f}%{rn['he'] / rn['total'] * 100:11.0f}%")
    print("-" * 74)
    dT = float(R["lngp"]["dTmin"])
    if dT < 2.95:
        print(f"Note: the main exchanger approach is {dT:.1f} K, below the 3 K design value.")
        print("      The refrigerant recipe was designed for the design feed and should be re-tuned for this one.")
    print(f"Whole plant solved in {dt:.0f} s, no Aspen Plus.\n")


if __name__ == "__main__":
    main()
