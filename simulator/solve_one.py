"""
Solve the whole LBNH plant for one feed and print the result as JSON.
Used by live_app.py (the live simulator page).

    python solve_one.py FLOW_TPH HELIUM_WT_PCT NITROGEN_WT_PCT
"""
import sys, json, time, io, contextlib
import change_feed as CF


# Operating cost of the AI-optimized design ($ million per year, the paper's cost basis):
# feed gas scales with the feed flow, electricity with the total power, the rest is fixed.
RAW0, POWER_COST0, OTHER0, POWER0 = 898.776, 139.554, 289.291, 216.378


def money(d, flow):
    r = CF.revenue(d)
    r["opex"] = RAW0 * flow / CF.F0 + POWER_COST0 * d["total"] / POWER0 + OTHER0
    r["profit"] = r["total"] - r["opex"]
    return r


def main():
    flow, he, n2 = (float(x) for x in sys.argv[1:4])
    CF.set_feed(flow, he, n2)
    from plant import run
    t0 = time.time()
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        R = run(nrhr="cycle")
    p, w = R["products"], R["power_MW"]
    new = dict(lng=p["lng_tph"], he=p["he_kgh"], n2=p["n2_kmolh"] * 28.0 / 1000, LNGP=w["c3mr"],
               BOGC=w["bog_compressor"], NRHR=w["nru_condenser"] + w["he_separation"], HLS=w["helium_liquefier"])
    new["total"] = new["LNGP"] + new["BOGC"] + new["NRHR"] + new["HLS"]
    # boil-off gas: flow and helium content (helium concentrates here)
    from props import PKG, C
    MW = dict(zip(C, PKG.MW)); bc = R["bogc"]["bog_comp"]; F = R["bogc"]["bog_total_kmolh"]
    new["bog"] = sum(F * x * MW[c] for c, x in bc.items()) / 1000.0
    new["bog_he"] = F * bc.get("He", 0.0) * MW["He"] / 1000.0 / new["bog"] * 100.0
    new["feed"] = flow
    base = dict(CF.BASE); base["total"] = base["LNGP"] + base["BOGC"] + base["NRHR"] + base["HLS"]
    base.update(bog=68.12, bog_he=1.203, feed=CF.F0)
    out = dict(ok=True, flow=flow, he_wt=he, n2_wt=n2, seconds=time.time() - t0,
               base=base, new=new, rev_base=money(base, CF.F0), rev_new=money(new, flow),
               dTmin=float(R["lngp"]["dTmin"]), log=log.getvalue().splitlines()[-12:])
    print("@@JSON@@" + json.dumps(out, default=float))


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("@@JSON@@" + json.dumps(dict(ok=False, error=f"{type(e).__name__}: {str(e)[:200]}")))
