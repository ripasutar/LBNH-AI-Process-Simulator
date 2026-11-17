"""
AI tool vs Aspen Plus: side-by-side, with each gap classified.

Economics are compared ENGINE-NORMALISED: the paper's CAPEX/sales/OPEX
and the tool's simulated quantities both go through the same cash-flow
engine on the paper's Table 8 rules, so differences reflect the plants,
not the spreadsheets.
"""
import json
import numpy as np
from props import PKG, C
from econ_aspen_basis import cashflow, BASIS

MW = dict(zip(C, PKG.MW))
H = BASIS["hours"]

# ---------------- the paper (Aspen Plus V14) ----------------------------
A = dict(
    lng_tph=508.86 + 32.50, he_kgh=820.0, n2_tph=30.91,
    bog_tph=64.23, bog_He=1.28, bog_CH4=50.12, bog_N2=48.60,
    hhc_tph=600.00 - 570.07, rich_tph=32.50, rich_N2_wt=7.75,
    nru_ch4_rec=32.50 * 0.9224 / (64.23 * 0.5012),
    he_purity=100.00, mr_tph=1408.0,
    mr_wt=dict(N2=10.63, CH4=24.35, C2H6=39.93, C3H8=25.10),
    power=dict(LNGP=939.90e6 / 3.6e6, BOGC=(565215 + 1979334 + 171532) / 3.6e6,
               NRHR=(20359256 + 221537439) / 3.6e6,
               HLS=(2473775 + 4120625 + 6809231) / 3.6e6),
    capex=4499.61, capex_lng=4253.0, sales=2204.92, opex=1398.84,
    rev_lng=1486.71, rev_he=684.77, rev_n2=33.44, raw=898.78,
    util=217.23, power_cost=214.61)
A["power"]["total"] = sum(A["power"].values())
A["capex_he"] = A["capex"] - A["capex_lng"]

# ---------------- the tool ---------------------------------------------
def tool_case(tag):
  sfx = "_" + tag if tag else ""
  R = json.load(open(f"plant_result{sfx}.json"))
  E_blind = json.load(open(f"economics_result{sfx}.json"))
  p, w = R["products"], R["power_MW"]
  bc = R["bogc"]["bog_comp"]
  mwb = sum(bc[c] * MW[c] for c in bc)
  mr = R["inputs"]["mr"]; mwmr = sum(mr[c] * MW[c] for c in mr)
  rich_mw = 0.0096 * MW["N2"] + 0.9904 * MW["CH4"]
  T = dict(
      lng_tph=p["lng_tph"], he_kgh=p["he_kgh"], n2_tph=p["n2_kmolh"] * 28.0 / 1000,
      bog_tph=R["bogc"]["bog_total_kmolh"] * mwb / 1000,
      bog_He=bc["He"] * MW["He"] / mwb * 100, bog_CH4=bc["CH4"] * MW["CH4"] / mwb * 100,
      bog_N2=bc["N2"] * MW["N2"] / mwb * 100,
      hhc_tph=sum(v * MW[c] for c, v in R["lngp"]["hhc_comp"].items()) / 1000,
      rich_tph=R["nru"]["bottoms_kmolh"] * rich_mw / 1000,
      rich_N2_wt=R["nru"]["bottoms_N2"] * MW["N2"] / rich_mw * 100,
      nru_ch4_rec=R["nru"]["CH4_recovery"], he_purity=99.999,
      mr_tph=R["lngp"]["F_mr_kmolh"] * mwmr / 1000,
      mr_wt={c: mr[c] * MW[c] / mwmr * 100 for c in mr},
      power=dict(LNGP=w["c3mr"], BOGC=w["bog_compressor"],
                 NRHR=w["nru_condenser"] + w["he_separation"],
                 HLS=w["helium_liquefier"]))
  T["power"]["total"] = sum(T["power"].values())

  # ---------------- economics, engine-normalised -------------------------
  T["capex_he"] = E_blind["addon"]["capex"] / 1e6
  T["capex_lng"] = A["capex_lng"]              # same cited benchmark (Chen et al.)
  T["capex"] = T["capex_lng"] + T["capex_he"]
  rev_lng = T["lng_tph"] * H * BASIS["price_lng_t"] / 1e6
  rev_he = T["he_kgh"] / 1000 * H * BASIS["price_he_t"] / 1e6
  rev_n2 = T["n2_tph"] * H * BASIS["price_n2_t"] / 1e6
  raw = 600.0 * H * BASIS["price_ng_t"] / 1e6
  pw = T["power"]["total"] * 1e3 * H * BASIS["elec_kwh"] / 1e6
  cw = 1.47 + 1.15                                  # not modelled: paper's value
  sales = rev_lng + rev_he + rev_n2
  # the paper's cost structure, as ratios
  lab = 82.43 / 4499.61 * T["capex"]
  over = 0.5 * lab
  opc = 130.11 / 2204.92 * sales
  ga = 29.08 / 2204.92 * sales
  opex = raw + pw + cw + lab + over + opc + ga
  T.update(sales=sales, opex=opex, rev_lng=rev_lng, rev_he=rev_he, rev_n2=rev_n2,
           raw=raw, util=pw + cw, power_cost=pw)
  return T, cashflow(T["capex"], T["sales"], T["opex"]), R

eA = cashflow(A["capex"], A["sales"], A["opex"])
CASES = {}
for tag in ("", "cycle", "bogcold", "strong"):
    try:
        CASES[tag or "sealed"] = tool_case(tag)
    except FileNotFoundError:
        pass
PRIMARY = "cycle" if "cycle" in CASES else "sealed"
T, eT, R = CASES[PRIMARY]

out = dict(aspen=A, tool=T, econ_aspen=eA, econ_tool=eT, primary=PRIMARY,
           cases={k: dict(tool=v[0], econ=v[1]) for k, v in CASES.items()},
           econ_paper_reported=dict(npv=5958.67, irr=0.215, pbp=7.7, pi=2.26))
json.dump(out, open("comparison.json", "w"), indent=1, default=float)


def row(label, a, t, unit="", fmt="{:,.2f}", cls=""):
    d = (t - a) / a * 100 if a else float("nan")
    print(f"  {label:<34}{fmt.format(a):>12}{fmt.format(t):>12}{d:>+9.1f}%  {unit:<8}{cls}")


print(f"{'':36}{'Aspen':>12}{'AI tool':>12}{'diff':>10}")
print("PRODUCTS")
row("LNG incl. rich LNG, t/h", A["lng_tph"], T["lng_tph"], "", cls="design: heavy-HC removal")
row("liquid helium, kg/h", A["he_kgh"], T["he_kgh"])
row("nitrogen product, t/h", A["n2_tph"], T["n2_tph"])
print("STREAMS")
row("heavy HC removed, t/h", A["hhc_tph"], T["hhc_tph"], cls="design: scrubber vs drum")
row("BOG to NRHR, t/h", A["bog_tph"], T["bog_tph"])
row("  He in BOG, wt%", A["bog_He"], T["bog_He"])
row("  N2 in BOG, wt%", A["bog_N2"], T["bog_N2"])
row("rich LNG, t/h", A["rich_tph"], T["rich_tph"])
row("  N2 in rich LNG, wt%", A["rich_N2_wt"], T["rich_N2_wt"], cls="design: NRU pressure/spec")
row("NRU methane recovery, %", A["nru_ch4_rec"] * 100, T["nru_ch4_rec"] * 100)
row("mixed refrigerant, t/h", A["mr_tph"], T["mr_tph"], cls="design: MR optimised")
print("POWER, MW")
for k in ("LNGP", "BOGC", "NRHR", "HLS", "total"):
    row(k, A["power"][k], T["power"][k])
print("ECONOMICS, $M  (paper's prices and cost structure; same engine)")
row("CAPEX", A["capex"], T["capex"])
row("  helium section CAPEX", A["capex_he"], T["capex_he"], cls="independent estimates")
row("sales", A["sales"], T["sales"])
row("OPEX", A["opex"], T["opex"])
row("NPV (engine)", eA["npv"], eT["npv"])
row("IRR, %", eA["irr"] * 100, eT["irr"] * 100)
row("payback, y (simple)", eA["pbp"], eT["pbp"])
