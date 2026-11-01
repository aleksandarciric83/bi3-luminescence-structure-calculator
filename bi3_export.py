"""
bi3_export.py
=============
Files written by the export buttons of the Bi3+ calculator (no Streamlit): the PNG / SVG / PDF
figures, the tables, a short report, and the ZIP that holds them.
"""
import io
import re
import zipfile
from datetime import date

import bi3_format as fmt
import bi3_plots as plots

VIABLE_RULE = (
    "Viable site rule of the program: coordination number 5 or more, and the cation is not a network former "
    "(P, B, Si, S, V, Mo, W, As). The other cations are listed under 'Show All Cation Sites'. "
    "A viable site is a candidate, not a prediction that Bi³⁺ occupies it."
)

# Which band the experimental excitation maximum stored with a built-in benchmark belongs to, and whether the
# stored number is one value or taken from a published range. The numbers themselves are those of the
# engine; this table only describes them.
BENCHMARK_EXP_INFO = {
    "CaYAlO4 (ICDD 01-081-0742, I4/mmm)": ("A", "single value (authors' data)"),
    "CaYAlO4 (MP mp-1227044 JSON, I4mm)": ("A", "single value (authors' data)"),
    "YPO4 (ICSD-79754, Boutinaud Table 5 & 6)": ("C", "range 5.46-5.51 eV, S 0.34-0.43 eV (A band not observed)"),
    "Y3Al5O12 (YAG, ICSD-280104, Boutinaud Table 5)": ("A", "range 4.53-4.64 eV, S 0.49-0.58 eV"),
    "Y2O3 (ICSD-192862, S6 & C2 Sites, Boutinaud Table 5)": ("A", "range 3.31-3.32 eV, S 0.27-0.31 eV"),
    "Cs2NaYCl6 (ICSD-245353, Boutinaud Table 5)": ("A", "single value"),
    "CaF2 (ICSD-656448, Boutinaud Table 4 & 5)": ("A", "single value"),
    "SrF2 (ICSD-404414, Boutinaud Table 5)": ("A", "single value"),
    "BaF2 (ICSD-9009004, Boutinaud Table 5)": ("A", "single value"),
    "CaO (ICSD-7200686, Boutinaud Table 4 & 5)": ("A", "range 3.47-3.50 eV, S 0.34-0.40 eV"),
    "CaS (ICSD-603165, Boutinaud Table 4 & 5)": ("A", "range 3.01-3.02 eV, S 0.27 eV"),
    "CaSe (ICSD-41957, Boutinaud Table 4 & 5)": ("A", "single value"),
    "CaTiO3 (ICSD-165801, Boutinaud Table 3, 5 & 9)": ("CT", "range 3.65-3.87 eV, S 1.46-1.64 eV (A-type value: 3.35 eV, S 1.13 eV)"),
}


def ct_parameter_note(cn):
    """Note on the charge-transfer parameter row used at a site with coordination number `cn` ('' if none)."""
    try:
        cn = int(cn)
    except (TypeError, ValueError):
        return ""
    if cn == 11:
        return "CN = 11: interpolated parameters (the parameter set covers CN 5-10 and 12 only)"
    if cn > 12:
        return "CN > 12: the CN = 12 parameters are used"
    if cn < 5:
        return "CN < 5: the CN = 5 parameters are used"
    return ""


def safe_name(text):
    """File-name stem from a formula: letters, digits, dot, dash and underscore only."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(text)).strip("_") or "structure"


def figure_files(summary_pub, formula, exp_ex_nm, exp_em_nm):
    """PNG (600 dpi), SVG and PDF of both diagrams, from the same functions that draw them on screen."""
    figs = (
        ("energy_levels", plots.make_energy_diagram(summary_pub, formula, exp_ex_nm, exp_em_nm)),
        ("EA_vs_he", plots.make_power_law_plot(summary_pub, exp_ex_nm)),
    )
    out = {}
    for stem, fig in figs:
        for fmt in ("png", "svg", "pdf"):
            out["%s.%s" % (stem, fmt)] = plots.fig_bytes(fig, fmt)
    return out


def build_report(formula, space_group, cell_text, note, exp_ex_nm, exp_em_nm, version, ct_notes="", citation=""):
    base = safe_name(formula)
    citation = citation or fmt.citation_for_report()
    if exp_ex_nm and exp_em_nm:
        exp_text = "lambda_ex = %.1f nm (%.3f eV), lambda_em = %.1f nm (%.3f eV)" % (
            exp_ex_nm, plots.HC_EV_NM / exp_ex_nm, exp_em_nm, plots.HC_EV_NM / exp_em_nm)
    else:
        exp_text = "none"
    return "\n".join([
        "Bi3+ Luminescence-Structure Calculator, calculation report",
        "=" * 58,
        "Program version : %s" % version,
        "Date            : %s" % date.today().isoformat(),
        "Structure       : %s, space group %s" % (fmt.formula_unicode(formula), fmt.space_group_unicode(space_group)),
        "Unit cell       : %s" % cell_text,
        "Source note     : %s" % (note or "none"),
        "Experimental    : %s (maximum of the A-band excitation and of the emission)" % exp_text,
        "Site rule       : %s" % VIABLE_RULE,
        "CT parameters   : %s" % (ct_notes or "all sites use the tabulated charge-transfer parameters of their CN"),
        "",
        "Files",
        "-----",
        "%s_optical_summary.csv    optical transitions of all cation sites" % base,
        "%s_matrix_I_BVS.csv       coordination, K coefficients, bond valence sums" % base,
        "%s_subshells_he.csv       sub-shell dielectric and nephelauxetic table" % base,
        "%s_diagram_data.csv       numbers plotted in the two diagrams (viable sites)" % base,
        "%s_energy_levels.*        energy-level diagram (png 600 dpi, svg, pdf)" % base,
        "%s_EA_vs_he.*             E_A against h_e (png 600 dpi, svg, pdf)" % base,
        "",
        "How to cite",
        "-----------",
        "If results of this program are used, please cite (reference loaded from %s when this file was written):" % fmt.CITATION_URL,
        citation,
        "",
    ])


def build_zip(summary_pub, matrix_I_df, subshells_df, formula, space_group, cell_text, note,
              exp_ex_nm, exp_em_nm, version, citation=""):
    """One ZIP with the report, four CSV tables and the six figure files."""
    base = safe_name(formula)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        notes = "; ".join("%s: %s" % (r["Doping Site"], r["CT parameters"]) for _, r in summary_pub.iterrows()
                          if "CT parameters" in summary_pub and r["CT parameters"])
        z.writestr("%s_report.txt" % base,
                   build_report(formula, space_group, cell_text, note, exp_ex_nm, exp_em_nm, version, notes, citation))
        z.writestr("%s_optical_summary.csv" % base, fmt.plain_columns(summary_pub).to_csv(index=False))
        z.writestr("%s_matrix_I_BVS.csv" % base, fmt.plain_columns(matrix_I_df).to_csv(index=False))
        z.writestr("%s_subshells_he.csv" % base, fmt.plain_columns(subshells_df).to_csv(index=False))
        z.writestr("%s_diagram_data.csv" % base,
                   plots.diagram_table(summary_pub, exp_ex_nm, exp_em_nm).to_csv(index=False))
        for name, data in figure_files(summary_pub, formula, exp_ex_nm, exp_em_nm).items():
            z.writestr("%s_%s" % (base, name), data)
    return buf.getvalue()
