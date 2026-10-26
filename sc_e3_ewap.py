from pathlib import Path
import csv, gzip, json, math, collections, statistics as st
import numpy as np
from scipy import stats
import openpyxl
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

W   = "anaconda3/envs/prolif-env"
OUT = "anaconda3/envs/prolif-env/Scripts/-f/results"
Path(OUT).mkdir(parents=True, exist_ok=True)

GTEX = f"{W}/gtex.gct.gz"
HPAB = f"{W}/rna_tissue_consensus.tsv"
HPAS = f"{W}/rna_single_cell_type.tsv"
MMC1 = f"{W}/mmc1.xlsx"
CXG  = f"{W}/cxg_all"

# pseudocontagem por escala: pequena em relacao a unidade de cada fonte.
# TPM/nTPM/nCPM vao de 0 a milhares -> 0.1 e um piso de ruido razoavel.
# O valor reconstruido do CELLxGENE fica entre 1e-4 e 0.3 -> 0.1 afogaria o sinal.
EPS_TPM = 0.1
EPS_CXG = 1e-6


def carregar_gtex(p):
    g = collections.defaultdict(dict)
    with gzip.open(p, "rt") as fh:
        fh.readline(); fh.readline()
        cab = fh.readline().rstrip("\n").split("\t"); cols = cab[2:]
        for ln in fh:
            r = ln.rstrip("\n").split("\t")
            if len(r) < len(cab): continue
            d = g[r[1]]
            for i, c in enumerate(cols):
                try: d[c] = float(r[i+2])
                except ValueError: pass
    return g

def carregar_hpa(p):
    g = collections.defaultdict(dict)
    with open(p) as fh:
        fh.readline()
        for ln in fh:
            r = ln.rstrip("\n").split("\t")
            if len(r) >= 4:
                try: g[r[1]][r[2]] = float(r[3])
                except ValueError: pass
    return g

def carregar_cellxgene(pasta):
    g = collections.defaultdict(dict)
    for f in sorted(Path(pasta).glob("*.csv")):
        with open(f) as fh:
            linhas = [l for l in fh if not l.startswith("#")]
        for r in csv.DictReader(linhas):
            if not r["Expression"]: continue
            s = r["Gene Symbol"].split("_ENSG")[0]
            frac = float(r["Tissue Composition"].rstrip("%")) / 100
            g[s][r["Tissue"]] = frac * math.expm1(float(r["Expression"]))
    return g

def carregar_e3ome(p):
    wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
    linhas = list(wb["1A_E3-ome"].iter_rows(values_only=True)); cab = list(linhas[0])
    iG, iC = cab.index("Gene"), cab.index("Protein class")
    e3 = {}
    for r in linhas[1:]:
        if r[iG] and str(r[iG]).strip():
            s = str(r[iG]).strip()
            e3[s] = str(r[iC]).strip() if r[iC] else "Unclassified"
    return e3

def carregar_tabula_sapiens(p, tipo="hepatocyte"):
    wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
    out = {}
    for r in wb["1T_Single cell"].iter_rows(min_row=2, values_only=True):
        if r[0] and str(r[0]).strip() == tipo:
            out[str(r[2]).strip()] = r[4]
    return out

def contraste(valores, foco, eps, excluir=None):
    out = {}
    for g, d in valores.items():
        if foco not in d: continue
        o = [v for t, v in d.items() if t != foco and not (excluir and excluir(t))]
        if o: out[g] = math.log2((d[foco] + eps) / (st.mean(o) + eps))
    return out


gtex_val = carregar_gtex(GTEX)
hpa_val  = carregar_hpa(HPAB)
hsc_val  = carregar_hpa(HPAS)
cxg_val  = carregar_cellxgene(CXG)
e3ome    = carregar_e3ome(MMC1)
ts_hep   = carregar_tabula_sapiens(MMC1)

E3 = set(e3ome)

gtex_b = contraste(gtex_val, "Liver", EPS_TPM,
                   excluir=lambda t: t.split("_")[0] == "Liver")
hpa_b  = contraste(hpa_val,  "liver",       EPS_TPM)
hpa_sc = contraste(hsc_val,  "hepatocytes", EPS_TPM)
cxg_c  = contraste(cxg_val,  "liver",       EPS_CXG)
hpa_liv = {g: d["liver"] for g, d in hpa_val.items() if "liver" in d}

UNI = sorted(set(gtex_b) & set(hpa_b) & set(hpa_sc) & set(cxg_c) & E3)
len(UNI)


rho = lambda a, b, gs: stats.spearmanr([a[g] for g in gs], [b[g] for g in gs])[0]

FONTES = {"HPA tecido": hpa_b, "GTEx tecido": gtex_b,
          "CELLxGENE": cxg_c, "HPA célula única": hpa_sc}
MODAL  = {"HPA tecido": "b", "GTEx tecido": "b",
          "CELLxGENE": "s", "HPA célula única": "s"}

mat, nomes = [], list(FONTES)
for i in range(len(nomes)):
    for j in range(i+1, len(nomes)):
        a, b = nomes[i], nomes[j]
        k = "".join(sorted(MODAL[a] + MODAL[b]))
        mat.append({"a": a, "b": b, "rho": rho(FONTES[a], FONTES[b], UNI),
                    "tipo": {"bb": "BB", "bs": "BS", "ss": "SS"}[k], "alt": False})
u_ts = sorted(set(ts_hep) & set(hpa_sc))
rep_ts = rho(ts_hep, hpa_sc, u_ts)
mat.append({"a": "HPA célula única", "b": "Tabula Sapiens",
            "rho": rep_ts, "tipo": "SS", "alt": True})

sweep = []
for t in (0, 1, 2, 5, 10, 20, 30, 50):
    gs = [g for g in UNI if hpa_liv.get(g, 0) > t]
    if len(gs) < 25: continue
    sweep.append({"thr": t, "n": len(gs),
                  "bulk": rho(hpa_b, gtex_b, gs), "sc": rho(cxg_c, hpa_sc, gs)})


MAPA = [("fígado","liver","liver","Liver"),("rim","kidney","kidney","Kidney_Cortex"),
    ("pulmão","lung","lung","Lung"),("baço","spleen","spleen","Spleen"),
    ("pâncreas","pancreas","pancreas","Pancreas"),("testículo","testis","testis","Testis"),
    ("próstata","prostate gland","prostate","Prostate"),("ovário","ovary","ovary","Ovary"),
    ("estômago","stomach","stomach","Stomach"),("cólon","colon","colon","Colon_Transverse"),
    ("esôfago","esophagus","esophagus","Esophagus_Mucosa"),
    ("intestino delgado","small intestine","small intestine","Small_Intestine_Terminal_Ileum"),
    ("pele","skin of body","skin","Skin_Sun_Exposed_Lower_leg"),
    ("mama","breast","breast","Breast_Mammary_Tissue"),
    ("adrenal","adrenal gland","adrenal gland","Adrenal_Gland"),
    ("tec. adiposo","adipose tissue","adipose tissue","Adipose_Subcutaneous"),
    ("coração","heart","heart muscle","Heart_Left_Ventricle"),
    ("tuba uterina","fallopian tube","fallopian tube","Fallopian_Tube"),
    ("bexiga","urinary bladder","urinary bladder","Bladder"),
    ("medula espinhal","spinal cord","spinal cord","Brain_Spinal_cord_cervical_c-1"),
    ("córtex cerebral","brain","cerebral cortex","Brain_Cortex")]

TTR = {"fígado":"liver","rim":"kidney","pulmão":"lung","baço":"spleen","pâncreas":"pancreas",
  "testículo":"testis","próstata":"prostate","ovário":"ovary","estômago":"stomach",
  "cólon":"colon","esôfago":"esophagus","intestino delgado":"small intestine","pele":"skin",
  "mama":"breast","adrenal":"adrenal gland","tec. adiposo":"adipose tissue","coração":"heart",
  "tuba uterina":"fallopian tube","bexiga":"urinary bladder","medula espinhal":"spinal cord",
  "córtex cerebral":"cerebral cortex"}

tecidos = []
for nome, t_cx, t_hp, t_gt in MAPA:
    pre = t_gt.split("_")[0]
    Cc = contraste(cxg_val,  t_cx, EPS_CXG)
    H  = contraste(hpa_val,  t_hp, EPS_TPM)
    G  = contraste(gtex_val, t_gt, EPS_TPM, excluir=lambda t: t.split("_")[0] == pre)
    gs = [g for g in E3 if g in Cc and g in H and g in G]
    if len(gs) < 100: continue
    tecidos.append({"tecido": nome, "n": len(gs), "bb": rho(H, G, gs),
                    "bs_h": rho(H, Cc, gs), "bs_g": rho(G, Cc, gs)})
tecidos.sort(key=lambda r: r["bb"])


U4 = sorted(set(gtex_b) & set(hpa_sc) & set(hpa_liv))
alvo = [g for g in U4 if g in E3]
niveis = np.array([hpa_liv[g] for g in alvo])
pool = [g for g in U4 if g not in E3]
pv = np.array([hpa_liv[g] for g in pool]); ordem = pv.argsort(kind="stable")
pool_s = [pool[i] for i in ordem]; pv_s = pv[ordem]

rng = np.random.default_rng(42)
nulos = np.array([rho(gtex_b, hpa_sc,
                      [pool_s[rng.integers(max(0, i-25), min(len(pool_s), i+25))]
                       for i in np.searchsorted(pv_s, niveis)]) for _ in range(2000)])
r_e3 = rho(gtex_b, hpa_sc, alvo)
p_emp = float((np.sum(nulos <= r_e3) + 1) / (len(nulos) + 1))

ordenado = sorted([g for g in U4 if hpa_liv[g] > 0], key=lambda g: hpa_liv[g])
faixas = []
for b in np.array_split(np.array(ordenado), 8):
    gs = list(b)
    faixas.append({"lo": hpa_liv[gs[0]], "hi": hpa_liv[gs[-1]], "n": len(gs),
                   "bb": rho(hpa_b, gtex_b, gs), "bs": rho(gtex_b, hpa_sc, gs)})
r_e3, p_emp


FULL = 7.087
SURF = "#ffffff"; INK = "#000000"; INK2 = "#1a1a1a"; INK3 = "#4d4d4d"; GRID = "#c8c8c4"
BLUE = "#1a68c8"; RED = "#a32121"; GOLD = "#b8911f"; GREY = "#9a9a94"; ALT = "#3d3c38"
AX = 8.5; TK = 7.5; LG = 7.5; AN = 7.5; PL = 11.0
LW = 1.5; MS = 4.6; BARLW = .55

plt.rcParams.update({"font.size": TK, "axes.labelsize": AX, "xtick.labelsize": TK,
                     "ytick.labelsize": TK, "legend.fontsize": LG,
                     "font.family": "DejaVu Sans", "savefig.facecolor": SURF,
                     "figure.facecolor": SURF, "svg.fonttype": "path",
                     "axes.linewidth": BARLW, "xtick.major.width": BARLW,
                     "ytick.major.width": BARLW, "lines.solid_capstyle": "round"})

def tidy(ax, grid="both"):
    ax.set_facecolor(SURF)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID); ax.spines[s].set_linewidth(BARLW)
    ax.tick_params(colors=INK2, labelsize=TK, width=BARLW, length=2.6, pad=2.2)
    if grid: ax.grid(True, axis=grid, color=GRID, lw=.45, alpha=.9); ax.set_axisbelow(True)

def panel(ax, letra, x=-.085, y=1.03):
    ax.text(x, y, letra, transform=ax.transAxes, fontsize=PL, fontweight="bold",
            color=INK, va="bottom", ha="left")

def save(fig, nome):
    fig.savefig(f"{OUT}/{nome}.svg", facecolor=SURF, bbox_inches="tight", pad_inches=.03)
    print(f"{OUT}/{nome}.svg")

TR = {"HPA tecido":"HPA bulk","GTEx tecido":"GTEx bulk",
      "HPA célula única":"HPA single-cell","CELLxGENE":"CELLxGENE single-cell",
      "Tabula Sapiens":"Tabula Sapiens"}


fig = plt.figure(figsize=(FULL, 3.55), facecolor=SURF)
gs = gridspec.GridSpec(1, 2, width_ratios=[1.02, 1], wspace=.30,
                       top=.945, bottom=.275, left=.085, right=.988)

ax = fig.add_subplot(gs[0]); tidy(ax, grid="y")
x = np.arange(len(sweep)); xb = [s["bulk"] for s in sweep]; xs = [s["sc"] for s in sweep]
ax.plot(x, xb, "-o", color=BLUE, lw=LW, ms=MS, zorder=4, mec=SURF, mew=.7)
ax.plot(x, xs, "-o", color=RED,  lw=LW, ms=MS, zorder=4, mec=SURF, mew=.7)
ax.fill_between(x, xb, xs, color=INK3, alpha=.07, zorder=2)
for i, (a, b) in enumerate(zip(xb, xs)):
    ax.text(i, a+.040, f"{a:.2f}", ha="center", fontsize=AN-.5, color=BLUE, fontweight="bold")
    ax.text(i, b-.072, f"{b:.2f}", ha="center", fontsize=AN-.5, color=RED,  fontweight="bold")
ax.set_xticks(x)
ax.set_xticklabels([("all" if s["thr"] == 0 else f">{s['thr']}") + f"\n{s['n']}"
                    for s in sweep], fontsize=TK-.5, color=INK2, linespacing=1.25)
ax.set_xlabel("E3 ligases retained, by liver expression (nTPM)",
              fontsize=AX, color=INK2, labelpad=3)
ax.set_ylabel("agreement between two\nindependent resources · Spearman ρ",
              fontsize=AX, color=INK2, labelpad=3)
ax.set_ylim(0, 1.03); ax.set_xlim(-.5, len(sweep)-.5)
ax.text(.26, .70, "two BULK atlases\nhold steady", fontsize=AN, color=BLUE,
        fontweight="bold", va="center", linespacing=1.35)
ax.text(.26, .50, "two SINGLE-CELL atlases\nfall apart", fontsize=AN, color=RED,
        fontweight="bold", va="center", linespacing=1.35)
panel(ax, "A", x=-.135)

ax = fig.add_subplot(gs[1]); tidy(ax, grid="x")
ordem = sorted(mat, key=lambda m: m["rho"]); col = {"BB": BLUE, "BS": GOLD, "SS": RED}
for i, m in enumerate(ordem):
    c = ALT if m["alt"] else col[m["tipo"]]
    ax.barh(i, m["rho"], height=.60, color=c, zorder=3,
            hatch="////" if m["alt"] else None,
            edgecolor=SURF if m["alt"] else "none", linewidth=.4)
    ax.text(m["rho"]+.018, i, f"{m['rho']:.3f}", va="center", fontsize=AN,
            color=c, fontweight="bold")
ax.set_yticks(np.arange(len(ordem)))
ax.set_yticklabels([f"{TR[m['a']]}\n× {TR[m['b']]}" for m in ordem],
                   fontsize=TK-1, color=INK, linespacing=1.25)
ax.set_xlim(0, 1.0); ax.set_xlabel("Spearman ρ", fontsize=AX, color=INK2, labelpad=3)
ax.spines["left"].set_visible(False)
ax.legend(handles=[Line2D([], [], marker='s', ls='', ms=6, c=BLUE, label="bulk × bulk"),
                   Line2D([], [], marker='s', ls='', ms=6, c=GOLD, label="bulk × single-cell"),
                   Line2D([], [], marker='s', ls='', ms=6, c=RED,
                          label="single-cell × single-cell"),
                   Patch(facecolor=ALT, hatch="////", edgecolor=SURF,
                         label="single-cell × single-cell (independent dataset)")],
          loc="upper center", ncol=2, frameon=False, fontsize=LG-.5, labelcolor=INK2,
          handletextpad=.5, columnspacing=1.6, handlelength=1.4, bbox_to_anchor=(.42, -.165))
panel(ax, "B", x=-.40)
save(fig, "figure_08_modality")


n = len(tecidos); y = np.arange(n)
fig, ax = plt.subplots(figsize=(FULL*.72, 4.5), facecolor=SURF); tidy(ax, grid="x")
fig.subplots_adjust(top=.965, bottom=.155, left=.235, right=.975)
for i, r in enumerate(tecidos):
    lo, hi = min(r["bs_h"], r["bs_g"]), max(r["bs_h"], r["bs_g"])
    ax.plot([lo, r["bb"]], [i, i], color=GRID, lw=.8, zorder=1)
    ax.plot([lo, hi], [i, i], color=RED, lw=2.0, zorder=3, alpha=.32)
    ax.plot([lo, hi], [i, i], "o", color=RED, ms=3.6, zorder=4, mec=SURF, mew=.55)
    ax.plot(r["bb"], i, "o", color=BLUE, ms=5.0, zorder=5, mec=SURF, mew=.6)
mb = st.median([r["bb"] for r in tecidos])
ms_ = st.median([(r["bs_h"]+r["bs_g"])/2 for r in tecidos])
for v, c in ((mb, BLUE), (ms_, RED)):
    ax.axvline(v, color=c, lw=.75, ls=(0, (4, 2.5)), alpha=.7, zorder=2)
ax.text(mb,  n-.05, f"median {mb:.2f}",  color=BLUE, fontsize=AN-.5, ha="center", fontweight="bold")
ax.text(ms_, n-.05, f"median {ms_:.2f}", color=RED,  fontsize=AN-.5, ha="center", fontweight="bold")
ax.set_yticks(y); ax.set_yticklabels([TTR[r["tecido"]] for r in tecidos], fontsize=TK, color=INK)
ax.set_ylim(-.8, n+.5); ax.set_xlim(0, 1.0)
ax.set_xlabel("agreement between two independent resources · Spearman ρ",
              fontsize=AX, color=INK2, labelpad=3)
ax.spines["left"].set_visible(False)
ax.legend(handles=[Line2D([], [], marker='o', ls='', ms=5, c=BLUE,
                          label="bulk × bulk  (HPA × GTEx)"),
                   Line2D([], [], marker='o', ls='', ms=3.8, c=RED,
                          label="bulk × single-cell  (each atlas × CELLxGENE)")],
          loc="upper center", ncol=1, frameon=False, fontsize=LG-.5, labelcolor=INK2,
          handletextpad=.5, bbox_to_anchor=(.5, -.075))
save(fig, "figure_09_tissues")


fig = plt.figure(figsize=(FULL, 3.25), facecolor=SURF)
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.03], wspace=.27,
                       top=.945, bottom=.235, left=.078, right=.988)

ax = fig.add_subplot(gs[0]); tidy(ax, grid="y")
ax.hist(nulos, bins=46, color=GREY, edgecolor=SURF, lw=.25, zorder=3)
ax.axvline(r_e3, color=RED, lw=1.5, zorder=5)
ax.axvline(np.median(nulos), color=INK2, lw=.8, ls=(0, (4, 2.5)), zorder=4)
ax.set_xlabel("ρ of each matched control set · 2,000 draws", fontsize=AX, color=INK2, labelpad=3)
ax.set_ylabel("number of control sets", fontsize=AX, color=INK2, labelpad=3)
yl = ax.get_ylim()[1]
ax.text(r_e3+.003, yl*.98, f" E3-ome\n ρ = {r_e3:.3f}", color=RED, fontsize=AN,
        fontweight="bold", va="top", linespacing=1.3)
ax.text(np.median(nulos)-.042, yl*.50, f"control median\n{np.median(nulos):.3f}",
        color=INK2, fontsize=AN-1, va="top", ha="right", linespacing=1.3)
ax.text(.5, -.265, f"p = {p_emp:.2f} — E3 ligases fall within chance",
        transform=ax.transAxes, ha="center", fontsize=AX, color=INK, fontweight="bold")
panel(ax, "A", x=-.125)

ax = fig.add_subplot(gs[1]); tidy(ax, grid="y")
x = np.arange(len(faixas)); bb = [f["bb"] for f in faixas]; bs = [f["bs"] for f in faixas]
ax.plot(x, bb, "-o", color=BLUE, lw=LW, ms=MS-.6, zorder=4, mec=SURF, mew=.7)
ax.plot(x, bs, "-o", color=RED,  lw=LW, ms=MS-.6, zorder=4, mec=SURF, mew=.7)
ax.fill_between(x, bb, bs, color=INK3, alpha=.07, zorder=2)
for i, (a, b) in enumerate(zip(bb, bs)):
    ax.text(i, a+.038, f"{a:.2f}", ha="center", fontsize=AN-1.2, color=BLUE, fontweight="bold")
    ax.text(i, b-.070, f"{b:.2f}", ha="center", fontsize=AN-1.2, color=RED,  fontweight="bold")
fm = lambda v: "200k" if v >= 1000 else (f"{v:.0f}" if v >= 10 else
                                         (f"{v:.1f}" if v >= 1 else f"{v:.2f}"))
ax.set_xticks(x)
ax.set_xticklabels([f"{fm(f['lo'])}–\n{fm(f['hi'])}" for f in faixas],
                   fontsize=TK-1.5, color=INK2, linespacing=1.2)
ax.set_xlabel("liver expression bin (nTPM) · eight bins", fontsize=AX, color=INK2, labelpad=3)
ax.set_ylabel("agreement · Spearman ρ", fontsize=AX, color=INK2, labelpad=3)
ax.set_ylim(0, 1.03); ax.set_xlim(-.5, len(faixas)-.5)
ax.text(.15, .44, "bulk × bulk", color=BLUE, fontsize=AN, fontweight="bold")
ax.text(.15, .33, "bulk × single-cell", color=RED, fontsize=AN, fontweight="bold")
panel(ax, "B", x=-.125)
save(fig, "figure_10_controls")
