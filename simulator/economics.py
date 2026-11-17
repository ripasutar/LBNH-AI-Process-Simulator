"""
Economics of the LBNH complex
=============================

Two cases:

  WHOLE COMPLEX   LNGP + BOGC + NRHR + HLS built together.
  HELIUM ADD-ON   NRHR + HLS retrofitted to an existing LNG plant - the
                  case the paper describes as "retrofit-friendly". Its
                  benefit is the helium sold plus the boil-off methane that
                  is now re-liquefied and sold as LNG instead of being burnt
                  as fuel gas.

Capital, by the method that suits each part:

  LNG production and storage - a mature technology with published unit
  costs: $/tpa from Songhurst, OIES NG137 (2018), US lean gas, escalated
  by CEPCI.

  NRHR and HLS - bespoke equipment with no public benchmark at this size:
  bottom-up. Every compressor, motor and expander is costed with the
  Turton correlations from its SIMULATED power, split into parallel
  machines inside the correlation's validity range. Total equipment is
  then taken as rotating equipment / 0.40 (rotating machinery is about
  40 % of equipment cost in cryogenic plants), and fixed capital as
  equipment / 0.30, using the same equipment share of total capital that
  the OIES review reports for LNG plants.

Everything above that is a factor rather than a simulation result is in
design_basis.py and open to change from the dashboard.
"""
import json

import numpy as np
import CoolProp.CoolProp as CP

import design_basis as DB
from props import C

CEPCI_NOW = DB.CEPCI[2024]

# Turton, Analysis, Synthesis and Design of Chemical Processes, App. A
# (CEPCI 397, 2001):  log10 Cp = K1 + K2 log10 A + K3 (log10 A)^2
TURTON = {
    "compressor": dict(K=(2.2897, 1.3604, -0.1027), A=(450, 3000), FBM=2.8),
    "motor":      dict(K=(2.4604, 1.4191, -0.1798), A=(75, 2600),  FBM=1.5),
    "expander":   dict(K=(2.2476, 1.4965, -0.1618), A=(100, 1500), FBM=3.5),
}

HHV_MJ_KG = {"CH4": 55.50, "C2H6": 51.90, "C3H8": 50.35, "iC4": 49.36,
             "nC4": 49.50, "iC5": 48.95, "nC5": 49.01, "N2": 0.0, "He": 0.0}
MJ_PER_MMBTU = 1055.06


def turton(kind, A_kW):
    """Purchased cost (2024 $) and bare-module cost, splitting into
    parallel machines where the duty exceeds the correlation's range."""
    t = TURTON[kind]
    lo, hi = t["A"]
    n = max(1, int(np.ceil(A_kW / hi)))
    a = max(A_kW / n, lo)
    K1, K2, K3 = t["K"]
    la = np.log10(a)
    cp = 10 ** (K1 + K2 * la + K3 * la * la) * CEPCI_NOW / DB.CEPCI[2001]
    return n * cp, n * cp * t["FBM"], n


def hhv_mmbtu_per_t(comp_molefrac, MW):
    w = np.array([comp_molefrac[c] * MW[i] for i, c in enumerate(C)])
    w = w / w.sum()
    mj_kg = sum(w[i] * HHV_MJ_KG[c] for i, c in enumerate(C))
    return mj_kg * 1000.0 / MJ_PER_MMBTU


def cashflow(capex, annual_net, years=None, build=None, rate=None, tax=None,
             dep_y=None):
    """After-tax cash flows, NPV, IRR and payback."""
    years = years or DB.PLANT_LIFE_Y
    build = build or DB.CONSTRUCTION_Y
    rate = DB.DISCOUNT_RATE if rate is None else rate
    tax = DB.TAX_RATE if tax is None else tax
    dep_y = dep_y or DB.DEPRECIATION_Y
    spend = {3: [0.3, 0.4, 0.3], 2: [0.5, 0.5], 1: [1.0]}.get(build, [1 / build] * build)
    cf = [-capex * s for s in spend]
    dep = capex / dep_y
    for y in range(years):
        d = dep if y < dep_y else 0.0
        taxable = annual_net - d
        cf.append(annual_net - max(taxable, 0.0) * tax)
    cf = np.array(cf)
    t = np.arange(len(cf))
    npv = float(np.sum(cf / (1 + rate) ** t))
    # IRR by bisection
    def f(r):
        return np.sum(cf / (1 + r) ** t)
    irr = None
    lo, hi = -0.9, 5.0
    if f(lo) * f(hi) < 0:
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f(lo) * f(mid) <= 0:
                hi = mid
            else:
                lo = mid
        irr = 0.5 * (lo + hi)
    cum = np.cumsum(cf)
    pb = next((i - build + 1 for i, v in enumerate(cum) if v >= 0 and i >= build), None)
    return dict(npv=npv, irr=irr, payback_y=pb, cashflow=cf.tolist())


def evaluate(R, lng_price=None, he_price=None, gas_price=None,
             elec_price=None, lng_basis="greenfield", hours=None):
    from props import PKG
    lng_price = DB.LNG_PRICE_USD_MMBTU if lng_price is None else lng_price
    he_price = DB.HELIUM_PRICE_USD_KG if he_price is None else he_price
    gas_price = DB.FEED_GAS_USD_MMBTU if gas_price is None else gas_price
    elec_price = DB.ELECTRICITY_USD_KWH if elec_price is None else elec_price
    hours = hours or DB.OPERATING_HOURS
    MW = PKG.MW
    P, W = R["products"], R["power_MW"]

    # ---------------- capital: NRHR + HLS, bottom-up -----------------
    items = []
    def machine(label, kind, kW):
        cp, cbm, n = turton(kind, kW)
        items.append(dict(item=label, kind=kind, kW=kW, units=n,
                          purchased=cp, bare_module=cbm))
    machine("BOG compressor", "compressor", W["bog_compressor"] * 1e3)
    machine("BOG compressor motors", "motor", W["bog_compressor"] * 1e3)
    machine("NRU refrigeration compressor", "compressor", W["nru_condenser"] * 1e3)
    machine("NRU refrigeration motors", "motor", W["nru_condenser"] * 1e3)
    if W.get("he_separation", 0) > 0:
        machine("NRHR helium refrigeration compressor", "compressor", W["he_separation"] * 1e3)
        machine("NRHR helium refrigeration motors", "motor", W["he_separation"] * 1e3)
    nr = R.get("nrhr_refrig", {})
    if "W_exp_MW" in nr:
        machine("NRHR turbo-expanders", "expander", max(nr["W_exp_MW"] * 1e3, 100.0))
    machine("Helium cycle compressor", "compressor", W["helium_liquefier"] * 1e3)
    machine("Helium cycle motors", "motor", W["helium_liquefier"] * 1e3)
    # helium expanders: absorbed power from the cycle
    hl = R["hls"]
    exp_kW = max(0.15 * W["helium_liquefier"] * 1e3, 100.0)
    machine("Helium turbo-expanders", "expander", exp_kW)
    rotating = sum(i["purchased"] for i in items)
    equipment_he = rotating / 0.40
    fci_he = equipment_he / 0.30

    # ---------------- capital: LNG, benchmark ------------------------
    tpa = P["lng_tph"] * hours
    usd_tpa = DB.LNG_USD_PER_TPA_2018[lng_basis] * CEPCI_NOW / DB.CEPCI[2018]
    fci_lng = tpa * usd_tpa

    # ---------------- revenue and costs -------------------------------
    lng_hhv = hhv_mmbtu_per_t(P["lng_comp"], MW)
    lng_mmbtu = P["lng_tph"] * hours * lng_hhv
    rev_lng = lng_mmbtu * lng_price
    rev_he = P["he_kgh"] * hours * he_price
    feed_tph = {c: v for c, v in DB.FEED_TPH.items() if c in HHV_MJ_KG}
    feed_mmbtu = sum(v * HHV_MJ_KG[c] for c, v in feed_tph.items()) * 1000 / MJ_PER_MMBTU * hours
    cost_gas = feed_mmbtu * gas_price
    power_total = W["total"]
    cost_elec_total = power_total * 1e3 * hours * elec_price
    he_power = W["bog_compressor"] + W["nru_condenser"] + W["he_separation"] + W["helium_liquefier"]
    cost_elec_he = he_power * 1e3 * hours * elec_price

    # ---------------- WHOLE COMPLEX ----------------------------------
    fci_all = fci_lng + fci_he
    fixed_all = fci_all * (DB.MAINTENANCE_FRAC + DB.FIXED_OPEX_FRAC)
    net_all = rev_lng + rev_he - cost_gas - cost_elec_total - fixed_all
    whole = dict(capex=fci_all, capex_lng=fci_lng, capex_he=fci_he,
                 revenue_lng=rev_lng, revenue_he=rev_he, cost_feed=cost_gas,
                 cost_power=cost_elec_total, cost_fixed=fixed_all,
                 annual_net=net_all, **cashflow(fci_all, net_all))

    # ---------------- HELIUM ADD-ON ----------------------------------
    # Without NRHR the boil-off gas is fuel. With it, the methane in the
    # rich LNG is sold as LNG instead: its uplift is (LNG - fuel) price.
    rich_ch4_tph = R["bogc"]["bog_total_kmolh"] * R["bogc"]["bog_comp"]["CH4"] \
        * R["nru"]["CH4_recovery"] * 16.043 / 1000
    ch4_mmbtu = rich_ch4_tph * hours * HHV_MJ_KG["CH4"] * 1000 / MJ_PER_MMBTU
    uplift = ch4_mmbtu * (lng_price - gas_price)
    fixed_he = fci_he * (DB.MAINTENANCE_FRAC + DB.FIXED_OPEX_FRAC)
    net_he = rev_he + uplift - cost_elec_he - fixed_he
    addon = dict(capex=fci_he, revenue_he=rev_he, methane_uplift=uplift,
                 cost_power=cost_elec_he, cost_fixed=fixed_he,
                 annual_net=net_he, **cashflow(fci_he, net_he))

    return dict(
        assumptions=dict(lng_price=lng_price, he_price=he_price,
                         gas_price=gas_price, elec_price=elec_price,
                         lng_basis=lng_basis, usd_per_tpa_2024=usd_tpa,
                         hours=hours, discount=DB.DISCOUNT_RATE,
                         tax=DB.TAX_RATE, life=DB.PLANT_LIFE_Y),
        lng_hhv_mmbtu_t=lng_hhv, lng_mtpa=tpa / 1e6,
        equipment=items, rotating_purchased=rotating,
        equipment_he=equipment_he,
        whole=whole, addon=addon)


if __name__ == "__main__":
    import sys
    tag = sys.argv[1] if len(sys.argv) > 1 else ""
    R = json.load(open(f"plant_result{'_' + tag if tag else ''}.json"))
    E = evaluate(R)
    json.dump(E, open(f"economics_result{'_' + tag if tag else ''}.json", "w"), indent=1, default=float)
    a, w, h = E["assumptions"], E["whole"], E["addon"]
    M = 1e6
    print(f"LNG {E['lng_mtpa']:.2f} Mtpa, HHV {E['lng_hhv_mmbtu_t']:.2f} MMBtu/t")
    print("\nNRHR + HLS equipment (bottom-up, 2024 $)")
    for i in E["equipment"]:
        print(f"  {i['item']:<32}{i['kW']:>9.0f} kW  {i['units']:>2} units  purchased ${i['purchased']/M:6.2f} M")
    print(f"  rotating ${E['rotating_purchased']/M:.1f} M -> equipment ${E['equipment_he']/M:.1f} M -> fixed capital ${h['capex']/M:.0f} M")
    print(f"\nWHOLE COMPLEX  capex ${w['capex']/M:,.0f} M  (LNG ${w['capex_lng']/M:,.0f} M @ ${a['usd_per_tpa_2024']:.0f}/tpa, He ${w['capex_he']/M:.0f} M)")
    print(f"  revenue LNG ${w['revenue_lng']/M:,.0f} M/y  He ${w['revenue_he']/M:,.0f} M/y")
    print(f"  feed gas ${w['cost_feed']/M:,.0f} M/y  power ${w['cost_power']/M:,.0f} M/y  fixed ${w['cost_fixed']/M:,.0f} M/y")
    print(f"  NPV ${w['npv']/M:,.0f} M   IRR {w['irr']*100 if w['irr'] else float('nan'):.1f}%   payback {w['payback_y']} y")
    print(f"\nHELIUM ADD-ON  capex ${h['capex']/M:,.0f} M")
    print(f"  helium ${h['revenue_he']/M:,.0f} M/y  methane uplift ${h['methane_uplift']/M:,.0f} M/y  power ${h['cost_power']/M:,.0f} M/y  fixed ${h['cost_fixed']/M:,.0f} M/y")
    print(f"  NPV ${h['npv']/M:,.0f} M   IRR {h['irr']*100 if h['irr'] else float('nan'):.1f}%   payback {h['payback_y']} y")
