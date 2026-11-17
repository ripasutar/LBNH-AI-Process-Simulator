"""Process flow diagram in engineering style (fixed layout, stream tags) + stream tables, from streams_<case>.json.
    python draw_pfd.py   ->  LBNH_PFD.pdf (4 pages) and PFD_<case>.png"""
import json, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon, Circle, Rectangle
from matplotlib.backends.backend_pdf import PdfPages
HERE = os.path.dirname(os.path.abspath(__file__))
TITLE = {"my_design": "My design (Aspen settings)"}
INK, PROC, MRC, C3C, UTIL = "#26343C", "#1F4E79", "#9C2C77", "#B26A00", "#7F8C93"
SEC = "#12968B"

class PFD:
    def __init__(s, ax): s.ax = ax
    def box(s, x, y, w, h, title, sub="", fc="#EAF3FA", ec="#2F7FD0", r=0.12, bold=True, fs=9.5, ty=None):
        s.ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=1.4, zorder=3))
        s.ax.text(x, y + (ty if ty is not None else (0.12 if sub else 0)), title, ha="center", va="center", fontsize=fs, weight="bold" if bold else "normal", color=INK, zorder=4)
        if sub: s.ax.text(x, y - 0.22 - 0.1 * sub.count("\n"), sub, ha="center", va="center", fontsize=8, color="#35505E", zorder=4)
    def comp(s, x, y, w, h, title, power):
        s.ax.add_patch(Polygon([(x - w / 2, y + h / 2), (x + w / 2, y + h / 2), (x + w * 0.32, y - h / 2), (x - w * 0.32, y - h / 2)], closed=True, fc="#FFF1DC", ec="#B86E12", lw=1.4, zorder=3))
        s.ax.text(x, y + 0.13, title, ha="center", va="center", fontsize=8.5, color=INK, zorder=4)
        s.ax.text(x, y - 0.2, power, ha="center", va="center", fontsize=9.5, weight="bold", color="#8A4B00", zorder=4)
    def vessel(s, x, y, w, h, title, sub=""):
        s.ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle=f"round,pad=0,rounding_size={w / 2.2}", fc="#EEF7F4", ec=SEC, lw=1.4, zorder=3))
        s.ax.text(x, y + (0.1 if sub else 0), title, ha="center", va="center", fontsize=8.5, weight="bold", color=INK, zorder=4)
        if sub: s.ax.text(x, y - 0.2, sub, ha="center", va="center", fontsize=7.8, color="#35505E", zorder=4)
    def valve(s, x, y, label=""):
        d = 0.18
        s.ax.add_patch(Polygon([(x - d, y + d), (x + d, y - d), (x + d, y + d), (x - d, y - d)], closed=True, fc="white", ec=INK, lw=1.2, zorder=4))
        if label: s.ax.text(x, y + 0.32, label, ha="center", va="bottom", fontsize=7.8, color=INK, zorder=4)
    def line(s, pts, tag=None, at=None, col=PROC, dash=False, arrow=True, tagside=0):
        xs, ys = zip(*pts)
        s.ax.plot(xs, ys, color=col, lw=1.5, ls=(0, (4, 3)) if dash else "-", zorder=2, solid_capstyle="butt")
        if arrow:
            (x0, y0), (x1, y1) = pts[-2], pts[-1]
            s.ax.annotate("", xy=(x1, y1), xytext=(x1 - (x1 - x0) * 0.001 - (0.25 if x1 > x0 else -0.25 if x1 < x0 else 0), y1 - (0.25 if y1 > y0 else -0.25 if y1 < y0 else 0)),
                          arrowprops=dict(arrowstyle="-|>", color=col, lw=1.5, mutation_scale=13), zorder=5)
        if tag:
            k = at if at is not None else 0
            (x0, y0), (x1, y1) = pts[k], pts[k + 1]
            s.tag((x0 + x1) / 2, (y0 + y1) / 2, tag, col)
    def tag(s, x, y, t, col=PROC):
        s.ax.text(x, y, t, ha="center", va="center", fontsize=8.2, weight="bold", color=col, zorder=6,
                  bbox=dict(boxstyle="round,pad=0.22", fc="white", ec=col, lw=1.1))
    def text(s, x, y, t, **k):
        s.ax.text(x, y, t, ha=k.pop("ha", "center"), va=k.pop("va", "center"), fontsize=k.pop("fs", 9), color=k.pop("color", "#0F5E57"), zorder=4, **k)
    def section(s, x0, y0, x1, y1, title):
        s.ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle="round,pad=0,rounding_size=0.25", fc="none", ec=SEC, lw=1.1, ls=(0, (5, 4)), zorder=1))
        s.ax.text(x0 + 0.2, y1 - 0.12, title, ha="left", va="top", fontsize=11, weight="bold", color="#0F5E57", zorder=4)

def draw(case, pdf):
    d = json.load(open(os.path.join(HERE, f"streams_{case}.json"))); S, U, W = d["streams"], d["units"], d["power_MW"]
    fig, ax = plt.subplots(figsize=(24, 13.5)); ax.set_xlim(0, 24); ax.set_ylim(-0.3, 13.6); ax.axis("off")
    p = PFD(ax)
    fig.text(0.5, 0.975, "LBNH plant: process flow diagram generated from the AI simulator", ha="center", fontsize=18, weight="bold", color=INK)
    fig.text(0.5, 0.952, f"{TITLE[case]}  ·  feed 600 t/h, 0.137 wt % He  ·  total power {d['total_MW']:.1f} MW  ·  LNG {d['lng_tph']:.1f} t/h  ·  liquid helium {d['he_kgh']:.0f} kg/h  ·  stream conditions in the stream table",
             ha="center", fontsize=11.5, color="#35505E")
    # ---------- sections
    p.section(0.15, 7.9, 11.4, 13.3, "1 · LNG production (C3MR)")
    p.section(11.6, 7.9, 23.85, 13.3, "2 · LNG storage and boil-off gas (BOG)")
    p.section(0.15, -0.2, 19.6, 7.6, "3 · Nitrogen rejection and helium recovery (NRHR)")
    p.section(19.8, -0.2, 23.85, 7.6, "4 · Helium liquefaction (HLS)")
    # ---------- 1 LNG
    u = U["MCHE"]
    p.text(0.75, 10.15, "Natural gas\n600 t/h", fs=10, weight="bold")
    p.box(3.3, 10.45, 2.3, 1.25, "Propane pre-coolers", f"3 levels · {U['E-C3']['duty_MW']:.0f} MW\n+ heavy-HC removal")
    p.box(7.6, 10.2, 1.7, 3.4, "Main\ncryogenic\nheat\nexchanger", f"{u['duty_MW']:.0f} MW\nΔTmin {u['dTmin_K']:.1f} K\nUA {u['UA_MWK']:.1f} MW/K", fs=9.5, ty=0.6)
    p.comp(7.6, 12.6, 2.3, 0.9, "MR compressor, 3 stages", f"{W['MR']:.1f} MW")
    p.comp(2.6, 12.3, 2.1, 0.9, "Propane compressor", f"{W['C3']:.1f} MW")
    p.box(4.6, 12.3, 1.2, 0.55, "C3 condenser", "", fs=8.5)
    p.valve(9.9, 9.0, "JT valve")
    p.line([(1.25, 10.25), (2.15, 10.25)], "S1", 0)
    p.line([(4.45, 10.25), (6.75, 10.25)], "S2", 0)
    p.line([(3.3, 9.82), (3.3, 8.55), (4.6, 8.55)], "S3", 0); p.text(5.45, 8.55, "Heavy HC to\nfractionation", fs=9)
    p.line([(7.6, 11.9), (7.6, 12.15)], None, col=MRC); p.tag(8.1, 12.0, "S5", MRC)
    p.line([(6.7, 12.6), (6.0, 12.6), (6.0, 11.35), (4.8, 11.35), (4.8, 10.85), (4.45, 10.85)], "S6", 1, col=MRC)
    p.line([(4.45, 10.7), (6.75, 10.7)], "S7", 0, col=MRC)
    p.line([(8.45, 9.0), (9.72, 9.0)], "S8", 0, col=MRC)
    p.line([(10.08, 9.0), (10.5, 9.0), (10.5, 9.6), (8.45, 9.6)], "S9", 2, col=MRC)
    p.line([(3.65, 12.3), (4.0, 12.3)], None, col=C3C); p.line([(4.6, 12.02), (4.6, 11.5), (3.8, 11.5), (3.8, 11.08)], "C1", 1, col=C3C)
    p.line([(2.5, 11.08), (2.5, 11.85)], None, col=C3C); p.tag(1.95, 11.45, "C2–C4", C3C)
    p.line([(8.45, 11.3), (11.4, 11.3), (12.25, 11.3)], "S4", 0)
    # ---------- 2 storage
    p.valve(12.45, 11.3, "End flash")
    p.vessel(14.7, 10.9, 1.7, 1.7, "LNG tank", f"{S['11']['P_bar']:.1f} bar · {S['11']['T_C']:.0f} °C")
    p.ax.add_patch(Circle((17.4, 12.55), 0.32, fc="white", ec=INK, lw=1.2, zorder=4)); p.text(17.4, 12.55, "BOG\nheader", fs=7, color=INK)
    p.comp(19.9, 12.55, 2.0, 0.9, "BOG compressor, 3 st.", f"{W['BOG']:.1f} MW")
    p.valve(14.7, 8.6, "Rich-LNG valve")
    p.line([(12.63, 11.3), (13.85, 11.3)], "S11", 0)
    p.line([(12.45, 11.48), (12.45, 12.55), (17.08, 12.55)], "S10", 1)
    p.line([(14.7, 11.75), (14.7, 12.55)], None); p.tag(14.25, 12.15, "S12"); p.ax.add_patch(Circle((14.7, 12.55), 0.06, color=PROC, zorder=5))
    p.line([(15.55, 10.6), (16.6, 10.6), (16.6, 9.7), (18.6, 9.7)], "S15", 2); p.text(19.6, 9.7, "LNG product\nto ship", fs=10, weight="bold")
    p.line([(14.7, 8.78), (14.7, 10.05)], "S14", 0)
    p.line([(14.88, 8.6), (17.4, 8.6), (17.4, 12.23)], "S13", 1)
    p.ax.add_patch(Rectangle((17.3, 9.55), 0.2, 0.3, fc="white", ec="none", zorder=2.5))
    p.line([(17.72, 12.55), (18.9, 12.55)], "S16", 0)
    p.line([(20.95, 12.55), (22.4, 12.55), (22.4, 9.2)], "S17", 1); p.text(22.4, 8.85, "to NRHR\ncold box", fs=8.5, color=INK)
    p.line([(13.2, 8.6), (14.52, 8.6)], "S21", 0); p.text(12.55, 8.6, "from NRHR\ncold box", fs=8.5, color=INK)
    # ---------- 3 NRHR
    p.text(0.85, 5.15, "S17 from\nBOG\ncompressor", fs=8.5, color=INK)
    p.box(3.2, 4.0, 1.9, 3.0, "NRHR\ncold box", "plate-fin\nΔTmin 2 K", fc="#E6F3FB")
    t = U["T-NRU"]
    p.box(6.4, 4.0, 1.1, 4.2, "N₂\nrejection\ncolumn", f"20 stages\nboil-up {t['boilup']:.2f}\nreboiler\n{t['Q_reb_MW']:.1f} MW", fc="#EEF7F4", ec=SEC, fs=9)
    p.box(8.9, 6.6, 1.6, 0.8, "Overhead condenser", f"{t['condenser_C']:.1f} °C · {t['Q_cond_MW']:.1f} MW", fs=8.2)
    xs = (11.3, 13.4, 15.5)
    for k, x in enumerate(xs):
        p.vessel(x, 5.55, 0.95, 1.5, f"He sep. {k + 1}", f"{S[str(23 + 2 * k)]['T_C'] + 273.15:.0f} K")
    p.vessel(17.8, 5.55, 1.05, 1.5, "Adsorber", "He ≥ 99.999 %")
    p.valve(16.6, 2.6, "N₂ flash 1.3 bar")
    n2, he = U["K-N2"], U["K-He"]
    p.comp(2.25, 0.75, 2.6, 0.9, f"N₂ cycle {n2['P_hi_bar']:.0f}/{n2['P_lo_bar']:.1f} bar, 2 expanders", f"{W['N2']:.1f} MW")
    p.comp(5.0, 0.75, 2.4, 0.9, f"He cycle {he['P_hi_bar']:.1f}/{he['P_lo_bar']:.1f} bar", f"{W['He']:.2f} MW")
    p.line([(1.5, 5.0), (2.25, 5.0)], None)
    p.line([(4.15, 4.9), (5.85, 4.9)], "S18", 0)
    p.line([(6.4, 6.1), (6.4, 6.6), (8.1, 6.6)], "S19", 1)
    p.line([(6.4, 1.9), (6.4, 1.55), (3.7, 1.55), (3.7, 2.5)], "S20", 1)
    p.line([(2.25, 3.0), (1.3, 3.0)], "S21", 0); p.text(0.7, 3.0, "to rich-\nLNG valve", fs=8.5, color=INK)
    p.line([(2.25, 1.2), (2.25, 2.5)], None, col=UTIL, dash=True, arrow=False); p.line([(4.4, 1.2), (4.4, 1.35), (3.35, 1.35), (3.35, 2.5)], None, col=UTIL, dash=True, arrow=False)
    p.text(1.55, 1.85, "cold", fs=8, color=UTIL); p.line([(4.15, 4.2), (5.0, 4.2), (5.0, 7.1), (8.9, 7.1), (8.9, 7.0)], None, col=UTIL, dash=True, arrow=False); p.text(7.6, 7.25, "cold from cold box", fs=8, color=UTIL)
    p.line([(9.7, 6.6), (10.4, 6.6), (10.4, 5.95), (10.82, 5.95)], "S22", 0)
    for k in range(2):
        a, b = xs[k], xs[k + 1]
        p.line([(a, 6.3), (a, 6.9), (b - 0.75, 6.9), (b - 0.75, 5.95), (b - 0.47, 5.95)], f"S{23 + 2 * k}", 1)
    p.line([(xs[2], 6.3), (xs[2], 6.9), (17.0, 6.9), (17.0, 5.95), (17.27, 5.95)], "S27", 1)
    for k, x in enumerate(xs):
        p.line([(x, 4.8), (x, 2.6 if k < 2 else 2.6), (16.42, 2.6)] if k == 2 else [(x, 4.8), (x, 2.6)], f"S{24 + 2 * k}", 0, arrow=(k == 2))
    p.line([(xs[0], 2.6), (xs[2], 2.6)], None, arrow=False)
    p.line([(16.78, 2.6), (18.2, 2.6)], "S32", 0); p.text(19.0, 2.6, "Nitrogen\nproduct", fs=10, weight="bold")
    p.line([(17.8, 4.8), (17.8, 3.7)], "S30", 0); p.text(17.8, 3.35, "tail gas", fs=8.5, color=INK)
    # ---------- 4 HLS
    p.comp(21.6, 4.6, 2.6, 1.1, "Collins liquefier, 2 exp. + JT", f"{W['HLS']:.1f} MW")
    p.line([(18.33, 5.55), (21.6, 5.55), (21.6, 5.15)], "S29", 0)
    p.line([(21.6, 4.05), (21.6, 2.9)], "S31", 0); p.text(21.6, 2.35, "Liquid helium\n4.2 K", fs=10, weight="bold")
    # legend
    for k, (c, lab, ds) in enumerate([(PROC, "process stream", False), (MRC, "mixed refrigerant loop", False), (C3C, "propane loop", False), (UTIL, "refrigeration (cold) supply", True)]):
        ax.plot([12.0 + k * 2.9, 12.6 + k * 2.9], [0.15, 0.15], color=c, lw=2, ls=(0, (4, 3)) if ds else "-"); ax.text(12.7 + k * 2.9, 0.15, lab, va="center", fontsize=8.5, color=INK)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.94, bottom=0.01)
    fig.savefig(os.path.join(HERE, f"PFD_{case}.png"), dpi=150); pdf.savefig(fig); plt.close(fig)

ORDER = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "C1", "C2", "C3", "C4", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32"]
COMPS = ["He", "N2", "CH4", "C2H6", "C3H8", "iC4", "nC4", "iC5", "nC5"]
def sid(i): return i if i.startswith("C") else f"S{i}"
def table_pages(case, pdf):
    d = json.load(open(os.path.join(HERE, f"streams_{case}.json"))); S = d["streams"]
    ids = [i for i in ORDER if i in S]
    hdr = ["Stream", "Description", "From → To", "T, °C", "P, bar", "Vap. frac.", "Flow, t/h", "kmol/h", "He", "N₂", "CH₄", "C₂H₆", "C₃H₈", "C₄+C₅"]
    widths = [0.045, 0.19, 0.22, 0.05, 0.05, 0.05, 0.06, 0.06, 0.045, 0.045, 0.045, 0.045, 0.045, 0.045]
    for part in (ids[:18], ids[18:]):
        fig, ax = plt.subplots(figsize=(24, 13.5)); ax.axis("off")
        fig.text(0.5, 0.965, f"Stream table: {TITLE[case]}  (composition in mol %)", ha="center", fontsize=17, weight="bold", color=INK)
        cell = []
        for i in part:
            s = S[i]; x = s["x"]; heavy = sum(x[c] for c in ("iC4", "nC4", "iC5", "nC5"))
            f2 = lambda v: (f"{v * 100:.3f}" if v * 100 < 1 else f"{v * 100:.2f}") if v > 5e-6 else "–"
            cell.append([sid(i), s["name"], f"{s['frm']} → {s['to']}", f"{s['T_C']:.1f}", f"{s['P_bar']:.2f}", "–" if s["vf"] is None else f"{s['vf']:.3f}",
                         f"{s['tph']:.3f}" if s["tph"] < 10 else f"{s['tph']:,.1f}", f"{s['kmolh']:,.1f}"] + [f2(x[c]) for c in ("He", "N2", "CH4", "C2H6", "C3H8")] + [f2(heavy)])
        tb = ax.table(cellText=cell, colLabels=hdr, colWidths=widths, loc="upper center", cellLoc="center", bbox=[0, 0.93 - 0.048 * (len(part) + 1), 1, 0.048 * (len(part) + 1)])
        tb.auto_set_font_size(False); tb.set_fontsize(10)
        for (r, c), ce in tb.get_celld().items():
            ce.set_edgecolor("#C9D6DC"); ce.set_linewidth(0.6)
            if r == 0: ce.set_facecolor("#0F5E57"); ce.get_text().set_color("white"); ce.get_text().set_weight("bold")
            else:
                ce.set_facecolor("#F4F7F8" if r % 2 else "white")
                if c == 0: ce.get_text().set_weight("bold")
                if c in (1, 2): ce.get_text().set_ha("left"); ce.PAD = 0.02
        fig.text(0.01, 0.01, "Generated by the AI simulator (Peng–Robinson; CoolProp for pure N₂/He). LNG section on my Aspen settings; helium part in the AI simulator layout. Mass balances of every section close within 0.01 t/h.", fontsize=10, color="#6B8592")
        pdf.savefig(fig); plt.close(fig)

if __name__ == "__main__":
    with PdfPages(os.path.join(HERE, "LBNH_PFD.pdf")) as pdf:
        draw("my_design", pdf); table_pages("my_design", pdf)
    print("ok")
