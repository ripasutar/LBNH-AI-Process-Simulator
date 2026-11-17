"""
LBNH Process Twin - live simulator page.

    python live_app.py

Starts a small web page on this computer (http://localhost:8765) where you can
change the feed, press Solve, and the whole plant is solved again by the Python
model, with no Aspen Plus. Close the black window to stop it.
"""
import json, os, subprocess, sys, threading, time, webbrowser
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import change_feed as CF

PORT = 8765
LOCK = threading.Lock()


def exact(v, shown, exact_value):
    """The page shows rounded design values; treat them as the exact design value."""
    return exact_value if abs(v - shown) < 1e-9 else v


def solve(flow, he, n2):
    flow = exact(flow, 600.0, CF.F0)
    he = exact(he, round(CF.HE0, 3), CF.HE0)
    n2 = exact(n2, round(CF.N20, 2), CF.N20)
    for name, v in (("flow", flow), ("he", he), ("n2", n2)):
        lo, hi = CF.LIMITS[name]
        if not lo - 1e-9 <= v <= hi + 1e-9:
            return dict(ok=False, error=f"{name} = {v:g} is outside the tested range {lo:g} to {hi:g}")
    try:
        r = subprocess.run([sys.executable, os.path.join(HERE, "solve_one.py"), str(flow), str(he), str(n2)],
                           cwd=HERE, capture_output=True, text=True, timeout=CF.TIME_LIMIT)
    except subprocess.TimeoutExpired:
        return dict(ok=False, error=f"The model did not finish this feed in {CF.TIME_LIMIT} s. "
                                    "This is an honest limit of the tool: there are feeds it cannot solve.")
    for line in r.stdout.splitlines():
        if line.startswith("@@JSON@@"):
            return json.loads(line[8:])
    return dict(ok=False, error="The solver stopped without a result. " + (r.stderr or "")[-300:])


# ---- My Aspen Plus base design (LBNH 0.14 wt %, from the Aspen file and report) ----
ASPEN = dict(feed=600.0, lng=541.36, he=820.0, n2=30.91, bog=64.23, bog_he=1.28,
             LNGP=261.082, BOGC=0.7545, NRHR=67.1935, HLS=3.5131)
ASPEN["total"] = ASPEN["LNGP"] + ASPEN["BOGC"] + ASPEN["NRHR"] + ASPEN["HLS"]          # 332.5 MW
ASPEN_MONEY = dict(lng=1486.71, he=684.77, n2=33.44)                                      # $ M / yr
ASPEN_MONEY["total"] = sum(ASPEN_MONEY.values())                                          # 2,204.9
RAW_A, POWER_A, OTHER_A = 898.78, 214.61, 1398.84 - 898.78 - 214.61                       # opex 1,398.8


def aspen_money(new, flow, ai_rb, ai_rn):
    m = {k: ASPEN_MONEY[k] * ai_rn[k] / ai_rb[k] for k in ("lng", "he", "n2")}
    m["total"] = m["lng"] + m["he"] + m["n2"]
    m["raw"] = RAW_A * flow / ASPEN["feed"]; m["elec"] = POWER_A * new["total"] / ASPEN["total"]; m["other"] = OTHER_A
    m["opex"] = m["raw"] + m["elec"] + m["other"]
    m["profit"] = m["total"] - m["opex"]
    return m


def aspen_basis(res):
    """Show results on my Aspen base design: start from the Aspen numbers and apply the change
    that the AI model finds when it re-solves the whole plant for the new feed."""
    if not res.get("ok"):
        return res
    B, N = res["base"], res["new"]
    new = {k: ASPEN[k] * N[k] / B[k] for k in ("lng", "he", "n2", "bog", "bog_he", "LNGP", "BOGC", "NRHR", "HLS")}
    new["total"] = new["LNGP"] + new["BOGC"] + new["NRHR"] + new["HLS"]
    new["feed"] = N["feed"]
    base = dict(ASPEN)
    rb = aspen_money(base, ASPEN["feed"], res["rev_base"], res["rev_base"])
    rn = aspen_money(new, N["feed"], res["rev_base"], res["rev_new"])
    out = dict(res, base=base, new=new, rev_base=rb, rev_new=rn,
               ai=dict(base=B, new=N, rev_base=res["rev_base"], rev_new=res["rev_new"]))
    return out


def base_result():
    """The design-feed result, shown when the page opens (no solving needed)."""
    import solve_one as SO
    base = dict(CF.BASE); base["total"] = base["LNGP"] + base["BOGC"] + base["NRHR"] + base["HLS"]
    base.update(bog=68.12, bog_he=1.203, feed=CF.F0)
    m = SO.money(base, CF.F0)
    r = dict(ok=True, he_wt=CF.HE0, n2_wt=CF.N20, seconds=None, dTmin=3.0, log=None,
                base=base, new=dict(base), rev_base=m, rev_new=dict(m))
    return aspen_basis(r)


CODE = [
    ("Basis", "design_basis.py", "Design basis: feed, operating conditions, prices and every assumption, in one place"),
    ("Basis", "props.py", "Components and the property package"),
    ("Thermodynamics", "eos.py", "Peng-Robinson / SRK cubic equation of state, with the helium treatment"),
    ("Thermodynamics", "flash.py", "Phase equilibrium: flash with stability test, bubble and dew points"),
    ("Units", "units.py", "Compressors, pumps, valves, mixers and other unit models"),
    ("1 LNG liquefaction", "lngp.py", "C3MR liquefaction: main cryogenic exchanger, MR compressors, propane cycle"),
    ("3 N2 rejection", "column.py", "Distillation column solver (Naphtali-Sandholm Newton, trust region, homotopy)"),
    ("3 N2 rejection", "nru.py", "Nitrogen rejection unit around the column"),
    ("3 N2 rejection", "coldbox.py", "Cold box with the N2 and He refrigeration cycles"),
    ("4 Helium", "hls.py", "Helium liquefier, Collins-type cycle"),
    ("Whole plant", "plant.py", "The whole plant: LNGP -> BOGC -> NRHR -> HLS with the boil-off recycle"),
    ("Economics", "economics.py", "Cost correlations, cash flow, NPV, IRR and payback"),
    ("Economics", "econ_aspen_basis.py", "Cash flow on the published paper's basis"),
    ("Checks", "same_design.py", "Same design as my Aspen file: every compressor recomputed"),
    ("Live simulator", "change_feed.py", "Change the feed and re-solve the plant"),
    ("Live simulator", "solve_one.py", "One full solve for the live page"),
    ("Live simulator", "live_app.py", "The live simulator web page server"),
]


def code_files():
    out = []
    for group, name, desc in CODE:
        path = os.path.join(HERE, *name.split("/"))
        if os.path.exists(path):
            txt = open(path, encoding="utf-8", errors="replace").read()
            out.append(dict(group=group, name=name, desc=desc, lines=txt.count("\n") + 1, text=txt))
    return out


def page():
    html = open(os.path.join(HERE, "live_page.html"), encoding="utf-8").read()
    return (html.replace("__HE__", f"{CF.HE0:.3f}").replace("__N2__", f"{CF.N20:.2f}")
                .replace("__BASE__", json.dumps(base_result()))
                .replace('"__MODE__"', '"laptop"'))


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def _send(self, code, body, ctype):
        b = body.encode("utf-8")
        self.send_response(code); self.send_header("Content-Type", ctype); self.send_header("Content-Length", str(len(b)))
        self.end_headers(); self.wfile.write(b)

    def do_GET(self):
        if self.path == "/code":
            return self._send(200, json.dumps(code_files()), "application/json")
        if self.path in ("/three.min.js", "/OrbitControls.js", "/plant3d.js"):
            b = open(os.path.join(HERE, self.path[1:]), "rb").read()
            self.send_response(200); self.send_header("Content-Type", "application/javascript"); self.send_header("Content-Length", str(len(b)))
            self.end_headers(); self.wfile.write(b); return
        if self.path == "/flowsheet/PFD_my_design.jpg":
            b = open(os.path.join(HERE, self.path[1:]), "rb").read()
            self.send_response(200); self.send_header("Content-Type", "image/jpeg"); self.send_header("Content-Length", str(len(b)))
            self.end_headers(); self.wfile.write(b); return
        if self.path.startswith("/aspen_flowsheet.jpg"):
            b = open(os.path.join(HERE, "aspen_flowsheet.jpg"), "rb").read()
            self.send_response(200); self.send_header("Content-Type", "image/jpeg"); self.send_header("Content-Length", str(len(b)))
            self.end_headers(); self.wfile.write(b); return
        self._send(200, page(), "text/html; charset=utf-8")

    def do_POST(self):
        if self.path != "/solve":
            return self._send(404, "{}", "application/json")
        n = int(self.headers.get("Content-Length", 0))
        try:
            q = json.loads(self.rfile.read(n) or b"{}")
            flow, he, n2 = float(q["flow"]), float(q["he"]), float(q["n2"])
        except Exception:
            return self._send(400, json.dumps(dict(ok=False, error="Please type numbers in all three boxes.")), "application/json")
        with LOCK:
            print(time.strftime("%H:%M:%S"), f"solving feed {flow:g} t/h, helium {he:g} wt %, nitrogen {n2:g} wt % ...", flush=True)
            res = aspen_basis(solve(flow, he, n2))
            print("          ", "done in %.0f s" % res["seconds"] if res.get("ok") else "not solved: " + res.get("error", "")[:80], flush=True)
        self._send(200, json.dumps(res, default=float), "application/json")


if __name__ == "__main__":
    srv = ThreadingHTTPServer(("127.0.0.1", PORT), Handler)
    url = f"http://localhost:{PORT}"
    print(f"LBNH live simulator running at {url}")
    print("Your browser should open by itself. Keep this window open; close it to stop.")
    threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
