"""
Overall economic optimization of the LNG section (C3MR), on MY ASPEN DESIGN.

Maximize annual profit  =  revenue - feed gas - electricity - other - annualized capital
With the feed fixed (600 t/h) and the LNG product fixed (-152 C leaving the main
exchanger), revenue, feed gas and "other" do not change, so maximizing profit is
the same as minimizing the total annual cost (TAC) of the section:

    TAC = electricity + CRF x (compressor capital + main-exchanger capital)

Decision variables (7)
    y_N2, y_CH4, y_C2H6   mixed-refrigerant composition (C3H8 = 1 - sum)
    P_lo                  MR compressor suction pressure, bar
    P_hi                  MR compressor discharge pressure, bar
    dT_warm               warm-end approach of the main exchanger, K
    T_pre                 propane precooling temperature, C
    (the MR flow is not free: the energy balance fixes it)

Constraints
    minimum temperature approach in the main exchanger  >= 2 K
    MR compressor suction all vapour (vapour fraction >= 0.999)
    each MR component >= 2 mol %, propane >= 2 mol %
    bounds on pressures and temperatures (below)
    LNG leaves at -152 C; feed fixed

Everything else is my Aspen design: feed and pre-treatment streams, LNG outlet
temperature, propane cycle condenser (12 bar) and upper levels, compressor
efficiencies (MR stage 1: 0.85 x 0.9, stages 2-3: 0.72 x 0.9, propane 0.72 x 0.9).
"""
import json, sys, time, numpy as np, CoolProp.CoolProp as CP
from props import PKG, C
from flash import flash_TP, flash_PS, flash_PH, mixture_H, mixture_S
from lngp import mche, precool_loads, propane_cycle, h_molar

K0 = 273.15; MW = PKG.MW
import os
S = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'aspen', 'streams.json')))
MAP = {'METHANE':'CH4','ETHANE':'C2H6','PROPANE':'C3H8','BUTANE-N':'nC4','BUTANE-I':'iC4','PENTAN-N':'nC5',
       'PENTAN-I':'iC5','HELIUM':'He','NITROGEN':'N2'}
def vecof(st):
    z = np.zeros(len(C))
    for k, v in S[st]['molefrac'].items():
        if k in MAP and v: z[C.index(MAP[k])] = v
    return z / z.sum()
z_ng = vecof('PURE-NG'); F_ng = S['PURE-NG']['tph']*1e6/3600/float(z_ng @ MW)
z_feed = vecof('NG'); F_feed = S['NG']['tph']*1e6/3600/float(z_feed @ MW)
mr_aspen = vecof('MR-O-MCX')
Z_NG0, F_NG0, Z_FEED0, F_FEED0 = z_ng.copy(), F_ng, z_feed.copy(), F_feed
HE_WT0, N2_WT0 = 0.1366643889268512, 7.63

def set_feed(flow=600.0, he_wt=HE_WT0, n2_wt=N2_WT0):
    """Change the plant feed: flow scales every stream; helium and nitrogen contents are
    scaled in the feed and in the treated gas entering the main exchanger (both pass the
    pre-treatment unchanged), the hydrocarbons make up the rest."""
    global z_ng, F_ng, z_feed, F_feed
    iHe, iN = C.index('He'), C.index('N2')
    def adj(z0):
        w = z0 * MW; w = w / w.sum() * 100
        w2 = w.copy(); w2[iHe] = w[iHe] * he_wt / HE_WT0; w2[iN] = w[iN] * n2_wt / N2_WT0
        rest = [k for k in range(len(C)) if k not in (iHe, iN)]
        w2[rest] *= (100 - w2[iHe] - w2[iN]) / w[rest].sum()
        z = w2 / MW; return z / z.sum()
    z_feed = adj(Z_FEED0); z_ng = adj(Z_NG0)
    F_feed = F_FEED0 * flow / 600 * float(Z_FEED0 @ MW) / float(z_feed @ MW)
    F_ng = F_NG0 * flow / 600 * float(Z_NG0 @ MW) / float(z_ng @ MW)

iN2, iC1, iC2, iC3 = (C.index(c) for c in ('N2', 'CH4', 'C2H6', 'C3H8'))

# ---------------- economic basis (assumptions stated, see sensitivity) ----------------
HOURS = 8322.0
ELEC = 0.0775                 # $/kWh, the paper's basis
CRF = 0.1 / (1 - 1.1 ** -20)  # 10 %, 20 years -> 0.1175 per year
COMP_COST = 1000.0            # $ per kW of shaft power, installed compressor + electric motor (assumption)
MCHE_COST_BASE = 100e6        # $ installed main exchanger at my Aspen design (assumption), linear in UA
P_C3 = lambda T: CP.PropsSI('P', 'T', T, 'Q', 0, 'Propane')
T_C3 = lambda P: CP.PropsSI('T', 'P', P, 'Q', 0, 'Propane') - K0
LEV_ASPEN = (T_C3(8e5), T_C3(4.2e5), T_C3(1.1e5)); TCOND = T_C3(12e5)
COLD_APPROACH = -34.0 - LEV_ASPEN[2]          # 5.8 K between the coldest propane level and the precooled gas

def ua_of(curves):
    q, Th, Tc = curves
    q = np.asarray(q); dT = np.maximum(np.asarray(Th) - np.asarray(Tc), 1e-3)
    dq = np.diff(q); a, b = dT[:-1], dT[1:]
    lm = np.where(np.abs(a - b) > 1e-6, (a - b) / np.log(a / b), a)
    return float(np.sum(np.abs(dq) / lm))      # W/K

def mr_compress(z, F, T_in, P_lo, P_hi):
    """3 stages, equal pressure ratios, intercooled to 40 C; Aspen efficiencies."""
    effs = (0.85 * 0.9, 0.72 * 0.9, 0.72 * 0.9); r = (P_hi / P_lo) ** (1 / 3)
    T, P, W = T_in, P_lo, 0.0
    for e in effs:
        b, x, y, _ = flash_TP(PKG, T, P, z); h1 = mixture_H(PKG, T, P, b, x, y); s1 = mixture_S(PKG, T, P, b, x, y)
        P2 = P * r; Ts, bs, xs, ys = flash_PS(PKG, P2, z, s1, lo=60, hi=800)
        W += F * (mixture_H(PKG, Ts, P2, bs, xs, ys) - h1) / e
        T, P = 40 + K0, P2
    return W

def evaluate(x, mche_cost=MCHE_COST_BASE, ua_base=None):
    yN2, yC1, yC2, P_lo, P_hi, dTw, T_pre = x
    yC3 = 1 - yN2 - yC1 - yC2
    out = dict(x=list(map(float, x)), yC3=float(yC3))
    if min(yN2, yC1, yC2, yC3) < 0.02:
        out.update(ok=False, why='composition'); return out
    mr = np.zeros(len(C)); mr[iN2], mr[iC1], mr[iC2], mr[iC3] = yN2, yC1, yC2, yC3
    try:
        r = mche(z_ng * F_ng, mr, P_ng=49.6e5, P_hi=P_hi * 1e5 - 0.31e5, P_lo=P_lo * 1e5, T_in=T_pre + K0,
                 T_ko=T_pre + K0 - 0.01, T_out=-152 + K0, dT_warm=dTw)
        if not r.get('F_mr'):
            out.update(ok=False, why=r.get('reason', 'mche')); return out
        F = r['F_mr']
        b, *_ = flash_TP(PKG, r['T_warm'], P_lo * 1e5, mr)
        W_mr = mr_compress(mr, F, r['T_warm'], P_lo * 1e5, P_hi * 1e5)
        lev = (LEV_ASPEN[0], LEV_ASPEN[1], T_pre - COLD_APPROACH)
        steps = (lev[0] + 4, lev[1] + 4, T_pre)
        loads = np.array(precool_loads(z_feed, F_feed, 50e5, 25 + K0, T_steps=steps)) + \
                np.array(precool_loads(mr, F, P_hi * 1e5, 40 + K0, T_steps=steps))
        if min(loads) < 0:
            out.update(ok=False, why='precool order'); return out
        W_c3 = propane_cycle(list(loads), T_evap_C=lev, T_cond_C=TCOND, eff=0.72 * 0.9)
        UA = ua_of(r['curves'])
    except Exception as e:
        out.update(ok=False, why=type(e).__name__); return out
    W = (W_mr + W_c3) / 1e6                                   # MW
    ub = ua_base or UA
    capex = COMP_COST * W * 1e3 + mche_cost * UA / ub         # $
    elec = W * 1e3 * HOURS * ELEC                             # $/yr
    tac = elec + CRF * capex
    out.update(ok=True, W=W, W_mr=W_mr / 1e6, W_c3=W_c3 / 1e6, dTmin=float(r['dTmin']), vf=float(b),
               F_tph=F * float(mr @ MW) * 3.6 / 1000, UA_MWK=UA / 1e6, elec=elec / 1e6, capex=capex / 1e6,
               tac=tac / 1e6, T_jt=r['T_jt'] - K0, T_warm=r['T_warm'] - K0)
    return out

# Aspen design point
X_ASPEN = [mr_aspen[iN2], mr_aspen[iC1], mr_aspen[iC2], 1.2, 42.0, 13.3, -34.0]
LB = np.array([0.02, 0.10, 0.10, 1.1, 30.0, 2.0, -40.0])
UB = np.array([0.25, 0.60, 0.65, 6.0, 60.0, 25.0, -25.0])

def penalized(x, mche_cost, ua_base, log=None):
    x = np.clip(x, LB, UB)
    r = evaluate(x, mche_cost, ua_base)
    if not r['ok']:
        f = 1e4
    else:
        pen = 0.0
        if r['dTmin'] < 2.0: pen += 200 * (2.0 - r['dTmin']) + 20
        if r['vf'] < 0.999: pen += 500 * (0.999 - r['vf']) + 20
        f = r['tac'] + pen
        r['feasible'] = pen == 0
    if log is not None: log.append(r)
    return f
