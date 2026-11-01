"""
bi3_plots.py
============
Figure builders for the Bi3+ Luminescence-Structure Calculator.

Pure matplotlib (no Streamlit): the same functions draw the figures on screen and
write the PNG / SVG / PDF files of the export, so both are identical.
"""
import io

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from bi3_format import formula_math

HC_EV_NM = 1239.84193  # h*c in eV nm

# Okabe-Ito colours (distinguishable with colour-vision deficiency)
C_A, C_C, C_MMCT, C_IVCT, C_EXP = "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#D55E00"

COL_EA = "E_A (eV) [Eq. 18]"
COL_EC = "E_C (eV) [Eq. 19]"
COL_MMCT = "E_MMCT (eV) [Eq. 4]"
COL_IVCT = "E_IVCT Pair (eV) [Eq. 4]"
COL_HE_BI = "he(Bi3+) [Calc]"

ST_A = r"$^{3}\mathrm{P}_{1}$"
ST_C = r"$^{1}\mathrm{P}_{1}$"
ST_G = r"$^{1}\mathrm{S}_{0}$"


def _viable(summary_df):
    return summary_df[summary_df["Viable Bi3+ Site"]].reset_index(drop=True)


def _spread(ys, gap):
    """Push label heights apart (bottom-up) so that neighbours are at least `gap` apart."""
    ys = np.asarray(ys, float)
    order = np.argsort(ys)
    out = ys.copy()
    for k in range(1, len(order)):
        i, j = order[k - 1], order[k]
        if out[j] - out[i] < gap:
            out[j] = out[i] + gap
    return out


def _clean(v):
    return None if (v is None or (isinstance(v, float) and np.isnan(v)) or pd.isna(v)) else float(v)


def diagram_table(summary_df, exp_ex_nm=None, exp_em_nm=None):
    """The numbers behind both diagrams, one row per viable site."""
    v = _viable(summary_df)
    t = pd.DataFrame({
        "Site": v["Doping Site"],
        "he(Bi3+)": v[COL_HE_BI],
        "E_A (eV)": v[COL_EA],
        "E_C (eV)": v[COL_EC],
        "E_MMCT (eV)": v[COL_MMCT],
        "E_IVCT pair (eV)": v[COL_IVCT],
    })
    t["Exp. excitation maximum (eV)"] = HC_EV_NM / exp_ex_nm if exp_ex_nm else np.nan
    t["Exp. emission maximum (eV)"] = HC_EV_NM / exp_em_nm if exp_em_nm else np.nan
    return t


def make_energy_diagram(summary_df, formula, exp_ex_nm=None, exp_em_nm=None):
    """Energy-level diagram of Bi3+ above the 1S0 ground state, one column per viable site."""
    v = _viable(summary_df)
    n = max(len(v), 1)
    use_exp = bool(exp_ex_nm and exp_em_nm)
    e_ex = HC_EV_NM / exp_ex_nm if use_exp else None
    e_em = HC_EV_NM / exp_em_nm if use_exp else None

    top = 8.0
    for _, r in v.iterrows():
        for c in (COL_EA, COL_EC, COL_MMCT, COL_IVCT):
            val = _clean(r[c])
            if val is not None:
                top = max(top, val + 0.9)
    if use_exp:
        top = max(top, e_ex + 0.9)

    fig = Figure(figsize=(2.55 * n + 3.3, 5.2), dpi=130, layout="constrained")
    ax = fig.subplots()
    gap = 0.040 * (top + 0.3)
    hw = 0.22  # half width of a level bar

    for i, r in v.iterrows():
        x0 = float(i)
        ea, ec = _clean(r[COL_EA]), _clean(r[COL_EC])
        ect, eiv = _clean(r[COL_MMCT]), _clean(r[COL_IVCT])
        ax.hlines(0.0, x0 - hw, x0 + hw, colors="black", linewidth=2.4, zorder=3)
        items = []  # (y, text, colour, bold)
        if ea is not None:
            ax.hlines(ea, x0 - hw, x0 + hw, colors=C_A, linewidth=2.6, zorder=3)
            items.append((ea, "%s  %.2f eV" % (ST_A, ea), C_A, True))
        if ec is not None:
            ax.hlines(ec, x0 - hw, x0 + hw, colors=C_C, linewidth=2.2, zorder=3)
            items.append((ec, "%s  %.2f eV" % (ST_C, ec), C_C, False))
        if ect is not None:
            ax.hlines(ect, x0 - hw, x0 + hw, colors=C_MMCT, linewidth=2.2, linestyles="--", zorder=3)
            items.append((ect, "MMCT  %.2f eV" % ect, C_MMCT, False))
        if eiv is not None:
            ax.hlines(eiv, x0 - hw, x0 + hw, colors=C_IVCT, linewidth=2.0, linestyles=":", zorder=3)
            items.append((eiv, "IVCT  %.2f eV" % eiv, C_IVCT, False))
        if use_exp:
            ax.plot([x0], [e_ex], marker="s", color="black", markersize=6, linestyle="none", zorder=4)
            ax.hlines(e_em, x0 - hw * 0.8, x0 + hw * 0.8, colors=C_EXP, linewidth=2.0, linestyles="--", zorder=3)
        if items:
            ys = [it[0] for it in items]
            ys_new = _spread(ys, gap)
            for (y, txt, col, bold), yn in zip(items, ys_new):
                if abs(yn - y) > 0.03:
                    ax.plot([x0 + hw, x0 + hw + 0.05, x0 + hw + 0.09], [y, yn, yn], color=col, linewidth=0.7, zorder=2)
                ax.text(x0 + hw + 0.10, yn, txt, va="center", ha="left", fontsize=8.5, color=col,
                        fontweight="bold" if bold else "normal")

    ax.set_xticks(np.arange(len(v)))
    ax.set_xticklabels(list(v["Doping Site"]), fontsize=10)
    ax.set_xlim(-0.55, (len(v) - 1) + 1.0)
    ax.set_ylim(-0.3, top)
    ax.set_ylabel("Energy above the %s ground state (eV)" % ST_G, fontsize=10)
    ax.set_title("Bi$^{3+}$ energy levels in %s" % formula_math(formula), fontsize=11, loc="left")
    ax.grid(axis="y", linestyle=":", alpha=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], color="black", lw=2.4, label="Ground state %s" % ST_G),
        Line2D([], [], color=C_A, lw=2.6, label="A band %s" % ST_A),
        Line2D([], [], color=C_C, lw=2.2, label="C band %s" % ST_C),
        Line2D([], [], color=C_MMCT, lw=2.2, ls="--", label="MMCT state"),
        Line2D([], [], color=C_IVCT, lw=2.0, ls=":", label="Bi$^{3+}$-Bi$^{3+}$ IVCT pair"),
    ]
    if use_exp:
        handles += [
            Line2D([], [], color="black", marker="s", ls="none", markersize=6, label="Experiment: excitation maximum"),
            Line2D([], [], color=C_EXP, lw=2.0, ls="--", label="Experiment: emission maximum"),
        ]
    fig.legend(handles=handles, loc="outside right upper", frameon=False, fontsize=8.5)
    return fig


def make_power_law_plot(summary_df, exp_ex_nm=None):
    """E_A against he(Bi3+): the universal power-law curve with the sites."""
    v = _viable(summary_df)
    fig = Figure(figsize=(7.4, 5.2), dpi=130, layout="constrained")
    ax = fig.subplots()
    he = np.linspace(0.01, 2.45, 300)
    ea = 9.41 * (1.0 - 0.555 * he ** 0.240)
    ax.plot(he, ea, color=C_A, linewidth=2.2, label=r"$E_A = 9.41\,(1 - 0.555\,h_e^{0.240})$", zorder=3)
    ax.fill_between(he, ea - 0.35, ea + 0.35, color=C_A, alpha=0.14, linewidth=0, label=r"$\pm 0.35$ eV (accuracy of the curve)")
    palette = ["#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#000000"]
    for i, r in v.iterrows():
        x, y = float(r[COL_HE_BI]), float(r[COL_EA])
        ax.scatter([x], [y], s=70, zorder=5, color=palette[i % len(palette)], edgecolor="black", linewidth=0.6,
                   label=r"%s: $h_e$ = %.2f, $E_A$ = %.2f eV" % (r["Doping Site"], x, y))
        if exp_ex_nm:
            ax.scatter([x], [HC_EV_NM / exp_ex_nm], marker="x", s=60, color=C_EXP, zorder=6, linewidth=1.8,
                       label="Experiment: excitation maximum" if i == 0 else None)
    ax.set_xlabel(r"Nephelauxetic function $h_e$(Bi$^{3+}$)", fontsize=10)
    ax.set_ylabel(r"A-band energy $E_A$(Bi$^{3+}$) (eV)", fontsize=10)
    ax.set_title(r"$E_A$ against $h_e$", fontsize=11, loc="left")
    ax.set_xlim(0.0, 2.45)
    ax.set_ylim(2.2, 9.6)
    ax.grid(True, linestyle=":", alpha=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.legend(*ax.get_legend_handles_labels(), loc="outside right upper", frameon=False, fontsize=8.5)
    return fig


def fig_bytes(fig, fmt="png", dpi=600):
    """Serialise a figure to PNG / SVG / PDF bytes (fonts embedded in PDF, text kept as outlines in SVG)."""
    buf = io.BytesIO()
    with matplotlib.rc_context({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "path"}):
        fig.savefig(buf, format=fmt, dpi=dpi, bbox_inches="tight")
    return buf.getvalue()
