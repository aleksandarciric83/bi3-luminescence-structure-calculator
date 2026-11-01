"""
bi3_app.py
==========
Interactive Browser & Local Desktop GUI (Streamlit) for Calculating Bi3+ Optical Properties
from Host Crystal Structure Data (CIF, JSON, or the built-in benchmarks).

Start with:
    Bi3_Calculator.exe          (Windows program with everything included, double-click)
    Start_Calculator.bat        (Windows, from the source files)
    ./start_calculator.sh       (macOS / Linux)
    streamlit run bi3_app.py    (any system with the requirements installed)

The calculation in bi3_engine.py is unchanged since version 1.0; later versions change the
presentation layer (figures, export, typesetting, help) only.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

import bi3_engine as eng
import bi3_export as export
import bi3_format as fmt
import bi3_plots as plots
import bi3_validate as validate

__version__ = "2.0.0"
HERE = Path(__file__).resolve().parent
VIABLE_RULE = export.VIABLE_RULE


def _page_icon():
    try:
        from PIL import Image
        return Image.open(HERE / "bi3_icon.png")
    except Exception:
        return "🔬"


st.set_page_config(
    page_title="Bi3+ Luminescence-Structure Calculator",
    page_icon=_page_icon(),
    layout="wide",
)
st.markdown(fmt.TABLE_CSS, unsafe_allow_html=True)


# ==============================================================================
# HELPERS (presentation only)
# ==============================================================================
def _show_df(df, compound_col=None):
    """Show a table over the full width with typeset headers (subscripts, italics, Greek letters)."""
    st.markdown(fmt.table_html(df, compound_col), unsafe_allow_html=True)


def _citation():
    """The citation record from the web (read again after a minute at most), or None when offline."""
    return fmt.fetch_citation()


_safe_name = export.safe_name


@st.cache_data(show_spinner=False)
def benchmark_suite_table():
    """Run every built-in benchmark structure once and tabulate the viable sites (cached)."""
    rows = []
    for b_key, b_obj in eng.get_boutinaud_2023_benchmarks().items():
        b_res = eng.run_boutinaud_2023_calculation(b_obj)
        sdf = b_res["summary_df"]
        for _, r in sdf[sdf["Viable Bi3+ Site"]].iterrows():
            rows.append({
                "Benchmark Compound": b_key,
                "Site [CN]": r["Doping Site"],
                "Q(M) [BVS]": r["Q(M) [BVS]"],
                "Calc he(M)": r["he(M) [Calc]"],
                "Pub he(M)": r["he(M) [Pub]"],
                "Calc he(Bi3+)": r["he(Bi3+) [Calc]"],
                "Pub he(Bi3+)": r["he(Bi3+) [Pub]"],
                "Calc E_A (eV)": r["E_A (eV) [Eq. 18]"],
                "Calc E_C (eV)": r["E_C (eV) [Eq. 19]"],
                "Calc E_MMCT (eV)": r["E_MMCT (eV) [Eq. 4]"],
                "Exp E_exc (eV)": r["Exp. E_exc (eV)"],
                "Exp. band type": export.BENCHMARK_EXP_INFO.get(b_key, ("", ""))[0],
                "Exp. value from": export.BENCHMARK_EXP_INFO.get(b_key, ("", ""))[1],
                "Calc Delta_E (eV)": r["Delta_E(Bi3+) (eV) [Eq. 5]"],
                "Emission Origin": r["Predicted Emission Origin"],
            })
    return pd.DataFrame(rows)


@st.cache_data(show_spinner=False)
def figure_files(summary_pub, formula, exp_ex_nm, exp_em_nm):
    """PNG (600 dpi), SVG and PDF of both diagrams (cached)."""
    return export.figure_files(summary_pub, formula, exp_ex_nm, exp_em_nm)


@st.cache_data(show_spinner=False)
def export_package(summary_pub, matrix_I_df, subshells_df, formula, space_group, cell_text,
                   note, exp_ex_nm, exp_em_nm, citation):
    """ZIP with figures, tables, diagram data and a report (cached; `citation` is part of the cache key)."""
    return export.build_zip(summary_pub, matrix_I_df, subshells_df, formula, space_group, cell_text,
                            note, exp_ex_nm, exp_em_nm, __version__, citation)


def _fig_downloads(files, stem, key):
    """Three download buttons (PNG / SVG / PDF) for one diagram."""
    cols = st.columns(3)
    specs = (("png", "image/png", "PNG, 600 dpi"), ("svg", "image/svg+xml", "SVG"), ("pdf", "application/pdf", "PDF"))
    for col, (fmt, mime, label) in zip(cols, specs):
        col.download_button(label, data=files["%s.%s" % (stem, fmt)], file_name="%s.%s" % (_safe_name(stem), fmt),
                            mime=mime, key="dl_%s_%s" % (key, fmt))


# ==============================================================================
# HEADER & INTRODUCTION
# ==============================================================================
st.title("Bi³⁺ Luminescence–Structure Calculator")
st.markdown(
    """
    **Automated calculation of nephelauxetic functions $h_e(M)$ & $h_e(\\text{Bi}^{3+})$, intra-ionic excitation bands ($E_A, E_C$),
    metal-to-metal charge transfer ($E_{\\text{MMCT}}$), $\\text{Bi}^{3+}\\text{–Bi}^{3+}$ dimer intervalence charge transfer ($E_{\\text{IVCT}}$),
    and emission-origin criteria ($\\Delta E(\\text{Bi}^{3+})$) directly from crystal structure data.**
    """
)
_cit, _cit_err = _citation()
if _cit is not None:
    st.caption("Please cite: " + fmt.citation_text(_cit))
else:
    st.caption("The citation could not be loaded (no internet connection?). The current reference is at "
               + fmt.CITATION_URL)

# ==============================================================================
# SIDEBAR: INPUT MODE SELECTION & EXPERIMENTAL SPECTROSCOPY INPUTS
# ==============================================================================
st.sidebar.header("1. Structure Input")
input_mode = st.sidebar.radio(
    "Choose Input Source:",
    options=[
        "Upload .cif File (Default)",
        "Upload .json File (MP / Custom)",
        "Built-In Benchmarks",
    ],
    index=0,
)

benchmarks_dict = eng.get_boutinaud_2023_benchmarks()
active_crystal = None

if input_mode == "Upload .cif File (Default)":
    uploaded_cif = st.sidebar.file_uploader("Upload Crystal Structure (.cif)", type=["cif"])
    if uploaded_cif is not None:
        cif_text = uploaded_cif.read().decode("utf-8", errors="ignore")
        cif_errors, cif_warnings = validate.check_cif_text(cif_text)
        try:
            if cif_errors:
                raise ValueError("; ".join(cif_errors))
            active_crystal = eng.parse_cif_string(cif_text, filename=uploaded_cif.name)
            built_errors = validate.check_crystal(active_crystal)
            if built_errors:
                active_crystal = None
                raise ValueError("; ".join(built_errors))
            st.sidebar.success(f"Parsed CIF: {fmt.formula_unicode(active_crystal.formula)} "
                               f"({fmt.space_group_unicode(active_crystal.space_group)})")
            for w in cif_warnings:
                st.sidebar.warning(w)
        except Exception as e:
            active_crystal = None
            st.sidebar.error(
                f"The CIF file cannot be used: {e}. Check that the file is a complete CIF with cell parameters, "
                f"symmetry operations and atom sites, or enter the structure through the .json route."
            )
    else:
        st.sidebar.info("No .cif file uploaded yet. Showing the default CaYAlO₄ (ICDD 01-081-0742, I4/mmm) structure. Upload any .cif file above to analyze a new crystal.")
        active_crystal = benchmarks_dict["CaYAlO4 (ICDD 01-081-0742, I4/mmm)"]

elif input_mode == "Upload .json File (MP / Custom)":
    uploaded_json = st.sidebar.file_uploader("Upload Structure (.json)", type=["json"])
    if uploaded_json is not None:
        try:
            raw_json = json.loads(uploaded_json.read().decode("utf-8", errors="ignore"))
            json_errors = validate.check_json(raw_json)
            if json_errors:
                raise ValueError("; ".join(json_errors))
            inspection = eng.inspect_json_for_missing_parameters(raw_json)
            if not inspection["is_complete"]:
                st.warning(
                    f"⚠️ **Missing Unit Cell Parameters Detected in `{uploaded_json.name}`:** "
                    f"The uploaded JSON is missing: `{', '.join(inspection['missing_fields'])}`. "
                    f"We have automatically estimated the unit cell volume ($V \\approx {inspection['estimated_V']}\\text{{ \\AA}}^3$) "
                    f"and lattice parameters from the JSON density ($\\rho = {inspection['density_g_cm3']}\\text{{ g/cm}}^3$) and atom count ({inspection['number_of_atoms']} atoms). "
                    f"**Please confirm or adjust the parameters below:**"
                )
            with st.expander("📐 Confirm / Edit Unit Cell Parameters for JSON Structure", expanded=not inspection["is_complete"]):
                c1, c2, c3, c4 = st.columns(4)
                u_a = c1.number_input("Lattice a (Å)", value=float(inspection.get("estimated_a", 3.666)), format="%.4f")
                u_b = c2.number_input("Lattice b (Å)", value=float(inspection.get("estimated_b", 3.666)), format="%.4f")
                u_c = c3.number_input("Lattice c (Å)", value=float(inspection.get("estimated_c", 12.000)), format="%.4f")
                u_Z = c4.number_input("Formula Units Z", value=int(inspection.get("estimated_Z", 2)), min_value=1, step=1)
                c5, c6, c7 = st.columns(3)
                u_al = c5.number_input("Angle α (°)", value=float(inspection.get("alpha", 90.0)), format="%.2f")
                u_be = c6.number_input("Angle β (°)", value=float(inspection.get("beta", 90.0)), format="%.2f")
                u_ga = c7.number_input("Angle γ (°)", value=float(inspection.get("gamma", 90.0)), format="%.2f")

            active_crystal = eng.build_crystal_from_custom_json(
                raw_json, a=u_a, b=u_b, c=u_c, alpha=u_al, beta=u_be, gamma=u_ga, Z=int(u_Z), name=uploaded_json.name
            )
            built_errors = validate.check_crystal(active_crystal)
            if built_errors:
                active_crystal = None
                raise ValueError("; ".join(built_errors))
            st.sidebar.success(f"Built structure from JSON: V = {active_crystal.V:.2f} Å³, Z = {active_crystal.Z}")
        except Exception as e:
            active_crystal = None
            st.sidebar.error(f"The JSON file cannot be used: {e}")
    else:
        st.sidebar.info("Upload a Materials Project or custom .json file (e.g. mp-1227044.json). Showing mp-1227044.json preset by default.")
        active_crystal = benchmarks_dict["CaYAlO4 (MP mp-1227044 JSON, I4mm)"]

else:
    selected_bench = st.sidebar.selectbox(
        "Select Benchmark Compound:",
        options=list(benchmarks_dict.keys()),
        index=0,
        format_func=fmt.compound_unicode,
    )
    active_crystal = benchmarks_dict[selected_bench]

st.sidebar.markdown("---")
st.sidebar.header("2. Experimental Spectroscopy (Optional)")
st.sidebar.caption("Provide experimental PLE & PL maxima to compute the Stokes shift $S$ and the emission-origin "
                   "criterion $\\Delta E(\\mathrm{Bi}^{3+})$.")

default_ex_nm = round(1239.84193 / active_crystal.exp_EA_eV, 1) if (active_crystal and active_crystal.exp_EA_eV) else 280.0
default_em_nm = (
    round(1239.84193 / (active_crystal.exp_EA_eV - active_crystal.exp_S_eV), 1)
    if (active_crystal and active_crystal.exp_EA_eV and active_crystal.exp_S_eV)
    else 398.0
)

use_exp = st.sidebar.checkbox("Include Experimental PLE / PL Peaks", value=True)
if use_exp:
    exp_ex_nm = st.sidebar.number_input(
        "Experimental excitation maximum $\\lambda_\\mathrm{ex}$ (nm)", value=float(default_ex_nm), min_value=120.0, max_value=600.0, step=1.0,
        help="Enter the maximum of the lowest-energy (A-band) excitation band. The band assignment is not checked by the program.")
    exp_em_nm = st.sidebar.number_input(
        "Experimental emission maximum $\\lambda_\\mathrm{em}$ (nm)", value=float(default_em_nm), min_value=200.0, max_value=900.0, step=1.0)
else:
    exp_ex_nm, exp_em_nm = None, None
st.sidebar.markdown("---")
st.sidebar.caption(f"Version {__version__}")

# ==============================================================================
# EXECUTE CALCULATION
# ==============================================================================
if active_crystal is None:
    st.info("Upload a structure file in the sidebar, or choose a built-in benchmark compound.")

if active_crystal is not None:
    results = eng.run_boutinaud_2023_calculation(
        active_crystal,
        custom_exp_ex_nm=exp_ex_nm,
        custom_exp_em_nm=exp_em_nm,
    )
    summary_df = results["summary_df"]
    ct_flags = summary_df["CN"].map(export.ct_parameter_note).where(summary_df["Viable Bi3+ Site"], "")
    if (ct_flags != "").any():  # flag only; no number changes
        summary_df["CT parameters"] = ct_flags
    matrix_I_df = results["matrix_I_df"]
    subshells_df = results["subshells_df"]
    summary_pub = summary_df.drop(columns=["Raw Site"], errors="ignore")

    # Top Metric Banner
    m1, m2, m3, m4 = st.columns(4)
    viable_count = int(summary_df["Viable Bi3+ Site"].sum())
    m1.markdown(fmt.card_html("Compound formula", fmt.formula_html(active_crystal.formula)), unsafe_allow_html=True)
    m2.markdown(fmt.card_html("Space group", fmt.space_group_html(active_crystal.space_group)), unsafe_allow_html=True)
    m3.markdown(fmt.card_html("Unit cell volume <i>V</i>", f"{active_crystal.V:.2f} Å<sup>3</sup> (<i>Z</i> = {active_crystal.Z})"),
                unsafe_allow_html=True)
    m4.markdown(fmt.card_html("Viable Bi<sup>3+</sup> doping sites", fmt.site_count_text(viable_count)), unsafe_allow_html=True)
    st.write("")

    if active_crystal.reference_note:
        st.info(f"📖 **Reference / Structure Note:** {active_crystal.reference_note}")

    # ==========================================================================
    # TABS FOR ALL DELIVERABLES
    # ==========================================================================
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "📊 1. Optical Transitions & ΔE(Bi³⁺)",
        "🧱 2. Matrices I & II (Bond Valence)",
        "⚛️ 3. Dielectric & Nephelauxetic (he)",
        "📈 4. Energy Diagrams & Export",
        "✅ 5. Benchmark Suite",
        "❓ 6. Help & How to cite",
    ])

    with tab1:
        st.subheader("Optical Transition Energy Map & Emission Origin Diagnosis")
        st.markdown(
            """
            * **Intra-ionic A band** ($^1S_0 \\to {}^3P_1$): $E_A(\\text{Bi}^{3+}) = 9.41\\,(1 - 0.555\\,h_e(\\text{Bi}^{3+})^{0.240})$ (accuracy $\\pm 0.35\\text{ eV}$)
            * **Intra-ionic C band** ($^1S_0 \\to {}^1P_1$): $E_C(\\text{Bi}^{3+}) = 3.236 + 2.29\\,(E_A(\\text{Bi}^{3+}) - 2.972)^{0.856}$
            * **Outer-sphere charge transfer** ($E_{\\text{MMCT}}$ and $E_{\\text{IVCT}}$): $E_{\\text{CT}}(\\text{Bi}^{3+}) = k_{\\text{CN}'}^{\\text{CN}}\\left[\\chi_{\\text{CN}}(\\text{Bi}^{3+}) - \\alpha_{\\text{CN}'}^{\\text{CN}}\\,\\chi_{\\text{CN}'}(\\text{acceptor})/d_{\\text{Bi–acceptor}}\\right]$, where $k$ and $\\alpha$ are tabulated for each pair of coordination numbers CN (Bi³⁺) and CN′ (acceptor)
            * **Emission-origin criterion:** $\\Delta E(\\text{Bi}^{3+}) = |E_{\\text{CT}}(\\text{Bi}^{3+}) - E_A(\\text{Bi}^{3+}) - S|$ (threshold: $\\Delta E < 1.39\\text{ eV} \\implies$ **A-band emission**; $\\Delta E > 1.41\\text{ eV} \\implies$ **CT emission**).
            """
        )
        viable_df = summary_df[summary_df["Viable Bi3+ Site"]].copy()
        _show_df(viable_df.drop(columns=["Raw Site", "Viable Bi3+ Site"]))
        st.caption(VIABLE_RULE)
        if "CT parameters" in viable_df and (viable_df["CT parameters"] != "").any():
            st.warning("Charge-transfer energies of the marked sites (column 'CT parameters') use parameters that are not "
                       "tabulated for their coordination number: CN = 11 is interpolated between CN = 10 and CN = 12, "
                       "CN > 12 uses the CN = 12 parameters.")

        with st.expander("Show All Cation Sites (Including Tetrahedral / Non-Doping Network Formers)"):
            _show_df(summary_df.drop(columns=["Raw Site"]))

    with tab2:
        st.subheader("Matrix I (Coordination Connectivity & Multiplicity) & Bond Valence Sums Q(M)")
        st.markdown(
            """
            * **Bond valence sum:** $Q(M) = \\sum_{i=1}^{N(M)} \\exp\\left(-\\frac{d(M\\text{–}X_i) - \\text{BVP}(M)}{0.37}\\right)$
            * **Effective anion charge:** $Q(X_i) = \\frac{K(M\\text{–}X_i)}{K(X_i\\text{–}M)} Q(M)$
            """
        )
        _show_df(matrix_I_df)

    with tab3:
        st.subheader("Step-by-Step Sub-Shell Dielectric & Nephelauxetic Calculation")
        st.markdown(
            r"""
            For every bond type M–X$_i$ (cation site M, anion site X$_i$) of a cation site, with $Z$ formula units per cell,
            $w$ the multiplicity of a site, $N_\text{M}$ and $N_{\text{X}_i}$ the total coordination numbers, $N(\text{M–X}_i)$
            and $N(\text{X}_i\text{–M})$ the numbers of M–X$_i$ bonds seen from each side, $d$ the mean bond length and
            $V$ the cell volume:

            * $K(\text{M–X}_i) = \dfrac{w_\text{M}}{Z}\,\dfrac{N(\text{M–X}_i)}{N_\text{M}}$ and
              $K(\text{X}_i\text{–M}) = \dfrac{w_{\text{X}_i}}{Z}\,\dfrac{N(\text{X}_i\text{–M})}{N_{\text{X}_i}}$
            * $Q(\text{X}_i) = \dfrac{K(\text{M–X}_i)}{K(\text{X}_i\text{–M})}\,Q(\text{M})$
            * $N^*(\text{M–X}_i) = \dfrac{Q(\text{M})}{N_\text{M}} + \dfrac{V_\text{N}\,K(\text{M–X}_i)\,Q(\text{M})}{V_\text{S}\,K(\text{X}_i\text{–M})\,N_{\text{X}_i}}$,
              with $V_\text{S}$ (valence state) and $V_\text{N}$ (valence number) of the anion element
            * $v = \dfrac{d^{3}}{D}$ with the bond-volume denominator $D = \sum_{\text{C}}\sum_{\text{X}} \dfrac{w_\text{C}\,N(\text{C–X})}{V}\,\bar d_{\text{CX}}^{\,3}$
              (the sum runs over all cation sites C and anion sites X; for the Bi³⁺ column $d$ is replaced by $d_\text{Bi} = d_\text{host} + \Delta d$ at the doped site)
            * $\Gamma_1 = 2.56\,\dfrac{1}{d}\left[\dfrac{K(\text{M–X}_i)\,N_\text{M} + K(\text{X}_i\text{–M})\,N_{\text{X}_i}}{K(\text{M–X}_i) + K(\text{X}_i\text{–M})}\right]^{2}$
            * $\Gamma_2 = Q(\text{M})\left(1 - \dfrac{V_\text{N}}{V_\text{S}}\right)$
            * $\Gamma_3 = \exp\!\left[-1.364\,d\,\left(\dfrac{N^*}{v}\right)^{1/6}\right]$
            * $f_\text{c} = \dfrac{1}{1 + \dfrac{(\Gamma_1\Gamma_2\Gamma_3)^{2}}{1579\,d^{-4.96}}}$
            * $\alpha = \dfrac{0.208\,f_\text{c}\,d^{4.96}\,N^*}{3 + 0.874\,f_\text{c}\,d^{4.96}\,N^*/v}$
            * $h_e = \sqrt{\sum_i N(\text{M–X}_i)\,f_\text{c}\,\alpha\,Q(\text{X}_i)^{2}}$, evaluated with the host distances for $h_e(\text{M})$
              and with $d_\text{Bi}$ for $h_e(\text{Bi}^{3+})$
            """
        )
        st.markdown(f"Bond-volume denominator $D$ of the unperturbed host: `{results['denom_host']:.5f}`")
        _show_df(subshells_df)

    with tab4:
        st.subheader("Energy Diagrams & Export")
        st.caption(
            "The experimental excitation and emission maxima apply to the compound as a whole and are drawn at every "
            "site for comparison. They are not assigned to a particular site."
        )
        fig1 = plots.make_energy_diagram(summary_pub, active_crystal.formula, exp_ex_nm, exp_em_nm)
        fig2 = plots.make_power_law_plot(summary_pub, exp_ex_nm)
        files = figure_files(summary_pub, active_crystal.formula, exp_ex_nm, exp_em_nm)

        col_plot1, col_plot2 = st.columns(2)
        with col_plot1:
            st.pyplot(fig1)
            _fig_downloads(files, "energy_levels", "energy")
        with col_plot2:
            st.pyplot(fig2)
            _fig_downloads(files, "EA_vs_he", "power")

        st.markdown("### 📥 Export")
        base = _safe_name(active_crystal.formula)
        cell_text = f"V = {active_crystal.V:.2f} A^3, Z = {active_crystal.Z}"
        zip_bytes = export_package(
            summary_pub, matrix_I_df, subshells_df, active_crystal.formula, active_crystal.space_group,
            cell_text, active_crystal.reference_note, exp_ex_nm, exp_em_nm,
            fmt.citation_text(_cit) if _cit is not None else fmt.citation_for_report())
        st.download_button(
            "Download everything (.zip): figures, tables, diagram data, report",
            data=zip_bytes, file_name=f"{base}_Bi3_results.zip", mime="application/zip", key="dl_zip")
        c_dl1, c_dl2, c_dl3, c_dl4 = st.columns(4)
        c_dl1.download_button(
            "Optical Summary (.csv)",
            data=fmt.plain_columns(summary_pub).to_csv(index=False).encode("utf-8"),
            file_name=f"{base}_Bi3_Optical_Summary.csv",
            mime="text/csv",
        )
        c_dl2.download_button(
            "Matrix I & BVS (.csv)",
            data=fmt.plain_columns(matrix_I_df).to_csv(index=False).encode("utf-8"),
            file_name=f"{base}_Matrix_I_BVS.csv",
            mime="text/csv",
        )
        c_dl3.download_button(
            "Sub-Shell Dielectric Table (.csv)",
            data=fmt.plain_columns(subshells_df).to_csv(index=False).encode("utf-8"),
            file_name=f"{base}_Subshells_he.csv",
            mime="text/csv",
        )
        c_dl4.download_button(
            "Diagram data (.csv)",
            data=plots.diagram_table(summary_pub, exp_ex_nm, exp_em_nm).to_csv(index=False).encode("utf-8"),
            file_name=f"{base}_diagram_data.csv",
            mime="text/csv",
        )

    with tab5:
        st.subheader("Live Validation Against Published Values")
        st.markdown(
            "This table runs all built-in benchmark crystal structures live through `bi3_engine.py` and compares "
            "the calculated $h_e(M)$, $h_e(\\text{Bi}^{3+})$, $E_A$, $E_{\\text{MMCT}}$, and $\\Delta E(\\text{Bi}^{3+})$ against the published values."
        )
        bench_df = benchmark_suite_table()
        _show_df(bench_df, compound_col="Benchmark Compound")
        st.caption(
            "Exp. $E_\\text{exc}$ is the experimental excitation maximum stored with each benchmark. Exp. band type says which band it "
            "belongs to: A = intra-ionic A band, C = intra-ionic C band (YPO₄, the A band is not observed), CT = charge-transfer "
            "band (CaTiO₃; the A-type value is 3.35 eV, $S$ = 1.13 eV). 'Exp. value from' marks numbers that were taken from a range "
            "in the literature; the stored number is one value inside that range."
        )
        st.download_button("Benchmark table (.csv)", data=fmt.plain_columns(bench_df).to_csv(index=False).encode("utf-8"),
                           file_name="Bi3_benchmark_suite.csv", mime="text/csv", key="dl_bench_csv")

    with tab6:
        st.subheader("Help & How to cite")
        help_path = HERE / "HELP.md"
        if help_path.exists():
            st.markdown(help_path.read_text(encoding="utf-8"))
        else:
            st.info("The help file HELP.md was not found next to the program.")
        st.markdown("### How to cite")
        st.markdown("If results of this program are used, please cite:")
        _data, _err = fmt.fetch_citation()  # read from the web every time this tab is drawn (reused for a minute)
        if _data is not None:
            try:
                st.code(fmt.citation_text(_data), language=None, wrap_lines=True)
            except TypeError:  # older Streamlit versions have no wrap_lines
                st.code(fmt.citation_text(_data), language=None)
            st.caption("Loaded live from " + fmt.CITATION_URL)
        else:
            st.warning("The citation could not be loaded because there is no internet connection or the server did not answer. "
                       "Please open " + fmt.CITATION_URL + " in a browser to read the current reference.")
