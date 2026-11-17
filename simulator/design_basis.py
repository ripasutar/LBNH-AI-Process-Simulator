"""
LBNH complex - design basis
===========================

Every number the simulation uses that is not computed by it. Written
BEFORE any Aspen result has been seen, from three sources only:

  GIVEN   stated by the user: Table 1 of the paper, the block diagram,
          and the technologies (C3MR; N2 and He refrigeration for the
          nitrogen rejection unit; self-cooling and expander cooling for
          the helium liquefier)
  CHOSEN  a design decision made here, with the reason beside it
  TYPICAL a standard literature value for this kind of plant

Anything marked CHOSEN or TYPICAL is a place where this design may differ
from the Aspen design, and a difference in the results there is a
difference in DESIGN rather than an error in the simulation. The
comparison has to keep those two apart.
"""

# ---------------------------------------------------------------------
# feed                                                        GIVEN
# ---------------------------------------------------------------------
FEED_TPH = {                     # tonne/h, Table 1
    "He": 0.82, "CH4": 499.73, "H2O": 0.01, "N2": 45.78, "C2H6": 41.68,
    "C3H8": 2.69, "nC4": 0.03, "iC4": 0.08, "nC5": 4.72, "iC5": 4.47,
}
FEED_T_C = 35.0          # TYPICAL   pretreated gas leaving the acid-gas unit
FEED_P_BAR = 60.0        # TYPICAL   pipeline / plant inlet pressure

# ---------------------------------------------------------------------
# thermodynamics                                              CHOSEN
# ---------------------------------------------------------------------
# Peng-Robinson with Boston-Mathias alpha everywhere. It matched the
# reference equations of state better than SRK for every component
# tested, so there is no reason to switch methods between subsystems.
# Helium: kappa = 0 (negative acentric factor would invert the alpha
# function), binary parameters fitted to GERG-2008.
# Pure-fluid refrigerant loops (propane, nitrogen, helium) use CoolProp's
# reference equations of state directly - no mixture there to model.
EOS = "PR"
ALPHA = "boston-mathias"

# ---------------------------------------------------------------------
# LNGP - pretreatment and C3MR liquefaction
# ---------------------------------------------------------------------
DEHY_WATER_REMOVAL = 1.0          # CHOSEN  molecular sieve: all water out

PRECOOL_T_C = -35.0               # TYPICAL propane precooling outlet
PROPANE_EVAP_T_C = -40.0          # TYPICAL coldest propane level
PROPANE_COND_T_C = 40.0           # TYPICAL air-cooled propane condenser
PROPANE_LEVELS = 3                # TYPICAL three pressure levels in C3MR

HHC_SEPARATOR_T_C = -60.0         # CHOSEN  heavy-hydrocarbon knockout part-
                                  #         way down the MCHE. At -35 C nothing
                                  #         condenses at 58 bar (tested), and
                                  #         propane cannot reach colder.

MCHE_OUTLET_T_C = -155.0          # TYPICAL subcooled LNG leaving the MCHE
MCHE_MIN_APPROACH_K = 3.0         # TYPICAL spiral-wound exchanger pinch
NG_PRESSURE_DROP_BAR = 5.0        # TYPICAL across precooling + MCHE

MR_COMPOSITION = {                # TYPICAL C3MR mixed refrigerant, mol
    "N2": 0.08, "CH4": 0.40, "C2H6": 0.40, "C3H8": 0.12,
}
MR_HIGH_P_BAR = 45.0              # TYPICAL
MR_LOW_P_BAR = 3.5                # TYPICAL
MR_STAGES = 2                     # TYPICAL two-stage MR compression

COMPRESSOR_EFF = 0.80             # TYPICAL isentropic, all compressors
EXPANDER_EFF = 0.85               # TYPICAL isentropic, all turbo-expanders
COOLER_T_C = 40.0                 # TYPICAL air cooler outlet

# ---------------------------------------------------------------------
# BOGC - end flash and tank boil-off
# ---------------------------------------------------------------------
END_FLASH_P_BAR = 1.2             # TYPICAL storage-tank pressure
TANK_BOG_RATE = 0.0005            # TYPICAL 0.05 % of LNG per day ...
TANK_BOG_PER_HOUR = TANK_BOG_RATE / 24.0   # ... expressed per hour
LNG_N2_SPEC = 0.01                # TYPICAL max 1 mol % N2 in sales LNG

# ---------------------------------------------------------------------
# NRHR - nitrogen rejection and helium recovery
# ---------------------------------------------------------------------
BOG_COMPRESSOR_P_BAR = 25.0       # CHOSEN  NRU operating pressure
NRU_STAGES = 20                   # CHOSEN
NRU_FEED_STAGE = 8                # CHOSEN
NRU_RICH_LNG_N2_MAX = 0.01        # CHOSEN  bottoms meets the LNG spec
HE_SEP_T_C = -193.0               # CHOSEN  first He/N2 partial condenser
                                  #         (80 K), cooled by the N2 loop
HE_SEP2_T_C = -203.0              # CHOSEN  second stage (70 K), cooled by
                                  #         the He refrigeration loop -
                                  #         just above N2's triple point
HE_PURIFIER = "adsorber"          # CHOSEN  final polish to 99.999 %
N2_LOOP_HIGH_P_BAR = 30.0         # TYPICAL nitrogen refrigeration loop
N2_LOOP_LOW_P_BAR = 5.0

# ---------------------------------------------------------------------
# HLS - helium liquefaction and storage
# ---------------------------------------------------------------------
# "self cooling, expander cooling, expander cooling": interpreted as a
# Collins-type cycle - two turbo-expander stages for refrigeration and a
# final Joule-Thomson (self-cooling) expansion to produce liquid.
HE_CYCLE_HIGH_P_BAR = 15.0        # TYPICAL
HE_CYCLE_LOW_P_BAR = 1.3          # TYPICAL
HE_PRODUCT_T_K = 4.22             # GIVEN   -269 C: liquid He at ~1 atm
HE_PRODUCT_P_BAR = 1.013
HLS_PRECOOL_BY_N2 = True          # GIVEN   recovered N2 supplies cold
HLS_MIN_APPROACH_K = 1.0          # TYPICAL plate-fin at helium temperatures

# ---------------------------------------------------------------------
# economics                                    TYPICAL unless noted
# ---------------------------------------------------------------------
PLANT_LIFE_Y = 20
CONSTRUCTION_Y = 3
DISCOUNT_RATE = 0.10
TAX_RATE = 0.21                   # US federal corporate
OPERATING_HOURS = 8000            # per year
ELECTRICITY_USD_KWH = 0.07
LNG_PRICE_USD_MMBTU = 10.0
HELIUM_PRICE_USD_KG = 30.0        # liquid helium, bulk contract
COST_INDEX_YEAR = 2024            # CEPCI basis

# economics - added with the cost model                        TYPICAL
FEED_GAS_USD_MMBTU = 3.0          # Henry Hub-type feed gas price
CEPCI = {2001: 397.0, 2018: 603.1, 2024: 800.0}
LNG_USD_PER_TPA_2018 = {"greenfield": 1050.0, "brownfield": 700.0}
                                  # Songhurst, OIES NG137 (2018), US lean gas
MAINTENANCE_FRAC = 0.03           # of fixed capital per year
FIXED_OPEX_FRAC = 0.02            # labour, insurance, overheads
DEPRECIATION_Y = 10               # straight line
