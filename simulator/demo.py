"""
LBNH Process Twin - class demo
Runs the whole LNG + helium complex and prints it next to the Aspen Plus results
from Sutar, Modi & Xu (I&EC Research).

    python demo.py            the design compared in the presentation
    python demo.py strong     the designs found by the strong optimisers
"""
import sys, time
from plant import run

ASPEN = dict(lng=541.36, he=820.0, n2=30.91, LNGP=261.08, BOGC=0.75, NRHR=67.19, HLS=3.72)

strong = len(sys.argv) > 1 and sys.argv[1] == "strong"
print("\nLBNH Process Twin: solving the whole plant (no Aspen) ...\n")
t0 = time.time()
kw = dict(nrhr="cycle")
if strong:
    kw.update(mr_file="mr_best_strong.json", hls_file="hls_best_strong.json")
R = run(**kw)
dt = time.time() - t0
p, w = R["products"], R["power_MW"]
ai = dict(lng=p["lng_tph"], he=p["he_kgh"], n2=p["n2_kmolh"] * 28.0 / 1000,
          LNGP=w["c3mr"], BOGC=w["bog_compressor"],
          NRHR=w["nru_condenser"] + w["he_separation"], HLS=w["helium_liquefier"])
rows = [("LNG, t/h", "lng"), ("Liquid helium, kg/h", "he"), ("Nitrogen, t/h", "n2"),
        ("Power LNG train, MW", "LNGP"), ("Power BOG compression, MW", "BOGC"),
        ("Power N2 rejection + He recovery, MW", "NRHR"), ("Power He liquefier, MW", "HLS")]
print("\n" + "-" * 72)
print(f"{'':38s}{'Aspen Plus':>12s}{'AI':>10s}{'diff':>10s}")
print("-" * 72)
for label, k in rows:
    a, t = ASPEN[k], ai[k]
    print(f"{label:38s}{a:12.2f}{t:10.2f}{(t - a) / a * 100:+9.1f}%")
ta = sum(ASPEN[k] for k in ("LNGP", "BOGC", "NRHR", "HLS"))
tt = sum(ai[k] for k in ("LNGP", "BOGC", "NRHR", "HLS"))
print(f"{'Total power, MW':38s}{ta:12.1f}{tt:10.1f}{(tt - ta) / ta * 100:+9.1f}%")
print("-" * 72)
print(f"Whole plant solved in {dt:.0f} s on one CPU core, free open-source Python.\n")
