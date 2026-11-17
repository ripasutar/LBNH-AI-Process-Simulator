"""Aspen vs Python stream comparison -> compare.json (used by the simulator page)."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
A = json.load(open(os.path.join(HERE, "..", "section_models", "aspen", "streams.json")))
P = json.load(open(os.path.join(HERE, "streams_my_design.json")))["streams"]
MAP = [("Natural gas feed", "NG", ["1"], "LNG"), ("Heavy hydrocarbons removed", "HEAVY-HC", ["3"], "LNG"), ("Treated gas to main exchanger", "PURE-NG", ["2"], "LNG"),
       ("LNG leaving main exchanger", "LNG", ["4"], "LNG"), ("MR vapour to MR compressor", "MR-O-MCX", ["5"], "LNG"), ("MR after compressor + cooler", "MR-O-CL3", ["6"], "LNG"),
       ("MR after propane pre-cooling", "MR-O-X22", ["7"], "LNG"), ("End-flash vapour", "FBOG", ["10"], "Storage"), ("LNG after end flash", "LNG--FLS", ["11"], "Storage"),
       ("Tank boil-off", "TBOG", ["12"], "Storage"), ("Total BOG to the NRHR", "BOG-FTJ1", ["16"], "Storage"), ("Methane returned to LNG", "METHANE", ["14"], "NRHR"),
       ("Nitrogen product", "N2PROD", ["32"], "NRHR"), ("Helium product gas", "HELIUM", ["29"], "NRHR"), ("Liquid helium product", "LIQ-HE", ["31"], "HLS")]
rows = []
for name, a, ps, sec in MAP:
    s = A[a]; p = [P[i] for i in ps]; tph = sum(x["tph"] for x in p)
    rows.append(dict(name=name, sec=sec, aspen=a, py="+".join("S" + i for i in ps), aT=s["T"], aP=s["P"], aF=s["tph"], pT=p[0]["T_C"], pP=p[0]["P_bar"], pF=tph,
                     aVF=s.get("vf"), pVF=p[0]["vf"]))
json.dump(rows, open(os.path.join(HERE, "compare.json"), "w"), indent=1)
for r in rows: print(f"{r['name']:32s} {r['aspen']:9s} {r['aT']:8.1f} {r['pT']:8.1f} | {r['aP']:6.2f} {r['pP']:6.2f} | {r['aF']:8.2f} {r['pF']:8.2f}")
