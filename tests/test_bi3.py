"""Regression and robustness tests for the Bi3+ Luminescence-Structure Calculator v2.0.

Run:  python -m pytest tests          (needs: pip install pytest)
"""
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path

import matplotlib as mpl
import pandas as pd
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg

import bi3_engine as eng
import bi3_export
import bi3_format
import bi3_plots
import bi3_validate

ROOT = Path(__file__).resolve().parent.parent
GOLDEN = Path(__file__).resolve().parent / "golden"
ENGINE_SHA256_V10 = "a6d7f29d07ae0bc796265e0ff0210dab80d7a988afa019d6c4ecc0c4f5e6d1a8"  # SHA-256 of bi3_engine.py of version 1.0
ENGINE_SHA256 = "fc65f35d46fd10999b14fc2883a70eff281dd254cd03bdeefbc75a83c29d0568"  # bi3_engine.py of version 1.2.0 (notes reworded, no literature names in the user-visible notes): version 1.0 plus one corrected data value and notes (no physics change)
CIF_DEMO = ROOT / "examples" / "CaYAlO4_ICDD_01_081_0742.cif"
JSON_DEMO = ROOT / "examples" / "mp-1227044_CaYAlO4.json"
HC = 1239.84193


def _structures():
    s = dict(eng.get_boutinaud_2023_benchmarks())
    s["UPLOAD demo CIF"] = eng.parse_cif_string(CIF_DEMO.read_text(encoding="utf-8"), filename="demo.cif")
    raw = json.loads(JSON_DEMO.read_text(encoding="utf-8"))
    ins = eng.inspect_json_for_missing_parameters(raw)
    s["UPLOAD demo JSON (estimated cell)"] = eng.build_crystal_from_custom_json(
        raw, a=ins["estimated_a"], b=ins["estimated_b"], c=ins["estimated_c"],
        alpha=90.0, beta=90.0, gamma=90.0, Z=int(ins["estimated_Z"]), name="demo.json")
    return s


STRUCTS = _structures()
NAMES = list(STRUCTS)


def _exp(c):
    ex = HC / c.exp_EA_eV if c.exp_EA_eV else None
    em = HC / (c.exp_EA_eV - c.exp_S_eV) if (c.exp_EA_eV and c.exp_S_eV) else None
    return ex, em


def _run(name):
    c = STRUCTS[name]
    ex, em = _exp(c)
    return eng.run_boutinaud_2023_calculation(c, custom_exp_ex_nm=ex, custom_exp_em_nm=em)


def _roundtrip(df):
    return pd.read_csv(io.StringIO(df.to_csv(index=False)))


def test_engine_is_unchanged():
    """Version 1.1.1 differs from 1.0 only by one stored published value (Y2O3 C2 he(Bi3+), 1.72 -> 1.75),
    reference notes and comments. Any other edit changes the hash and must be reviewed."""
    assert hashlib.sha256((ROOT / "bi3_engine.py").read_bytes()).hexdigest() == ENGINE_SHA256


@pytest.mark.parametrize("table,key,drop", [
    ("baseline_summary.csv", "summary_df", ["Raw Site"]),
    ("baseline_matrix_I.csv", "matrix_I_df", []),
    ("baseline_subshells.csv", "subshells_df", []),
])
def test_numbers_equal_version_1_0(table, key, drop):
    gold = pd.read_csv(GOLDEN / table)
    parts = []
    for name in NAMES:
        d = _run(name)[key].drop(columns=drop, errors="ignore").copy()
        d.insert(0, "Structure", name)
        parts.append(d)
    now = _roundtrip(pd.concat(parts, ignore_index=True))
    pd.testing.assert_frame_equal(now, gold, check_dtype=False, rtol=1e-9, atol=1e-9)


def test_published_he_values():
    """Calculated against stored published values (Boutinaud 2023, Tables 5 and 9), computed live."""
    worst = {"he(M)": (0.0, ""), "he(Bi3+)": (0.0, "")}
    n_M = n_Bi = 0
    for name, c in eng.get_boutinaud_2023_benchmarks().items():
        sdf = eng.run_boutinaud_2023_calculation(c)["summary_df"]
        for _, r in sdf.iterrows():
            for key, calc, pub in (("he(M)", "he(M) [Calc]", "he(M) [Pub]"), ("he(Bi3+)", "he(Bi3+) [Calc]", "he(Bi3+) [Pub]")):
                if pd.notna(r[pub]):
                    n_M += key == "he(M)"
                    n_Bi += key == "he(Bi3+)"
                    d = abs(r[calc] - r[pub])
                    if d > worst[key][0]:
                        worst[key] = (d, "%s %s: calc %.3f, published %.3f" % (name, r["Doping Site"], r[calc], r[pub]))
    assert (n_M, n_Bi) == (13, 12)
    assert worst["he(M)"][0] <= 0.0105, "worst he(M): %s" % (worst["he(M)"],)
    assert worst["he(Bi3+)"][0] <= 0.012, "worst he(Bi3+): %s" % (worst["he(Bi3+)"],)


def test_y2o3_c2_published_value_is_table_5():
    sdf = eng.run_boutinaud_2023_calculation(
        [v for k, v in eng.get_boutinaud_2023_benchmarks().items() if k.startswith("Y2O3")][0])["summary_df"]
    c2 = sdf[sdf["Doping Site"].str.startswith("Y2_C2")].iloc[0]
    assert c2["he(Bi3+) [Pub]"] == 1.75 and c2["he(M) [Pub]"] == 1.61


def test_every_benchmark_has_exp_band_info():
    assert set(bi3_export.BENCHMARK_EXP_INFO) == set(eng.get_boutinaud_2023_benchmarks())
    info = bi3_export.BENCHMARK_EXP_INFO
    assert info["YPO4 (ICSD-79754, Boutinaud Table 5 & 6)"][0] == "C"
    assert info["CaTiO3 (ICSD-165801, Boutinaud Table 3, 5 & 9)"][0] == "CT"


def test_ct_parameter_note_flags_cn_11():
    assert "interpolated" in bi3_export.ct_parameter_note(11)
    assert "CN = 12 parameters" in bi3_export.ct_parameter_note(14)
    assert all(bi3_export.ct_parameter_note(n) == "" for n in (5, 6, 7, 8, 9, 10, 12))
    t = eng.TABLE2_CT_PARAMS
    for k in t[11]:  # the CN = 11 row is the midpoint of CN = 10 and CN = 12 (within rounding)
        assert abs(t[11][k] - (t[10][k] + t[12][k]) / 2) <= 0.006


def _unbalanced(s):
    bad = [o + c for o, c in ("()", "[]", "{}") if s.count(o) != s.count(c)]
    return bad + (["$"] if s.count("$") % 2 else [])


@pytest.mark.parametrize("name", NAMES)
def test_labels_are_balanced(name):
    c = STRUCTS[name]
    res = _run(name)
    texts = [name, c.formula, c.space_group, c.reference_note or ""]
    for key in ("summary_df", "matrix_I_df", "subshells_df"):
        df = res[key]
        texts += [str(x) for x in df.columns]
        texts += [v for col in df.columns for v in df[col] if isinstance(v, str)]
    assert [t for t in texts if _unbalanced(t)] == []


def _text_overlaps(fig):
    cv = FigureCanvasAgg(fig)
    cv.draw()
    r = cv.get_renderer()
    ticks = set(sum([ax.get_xticklabels() + ax.get_yticklabels() for ax in fig.axes], []))
    labs = {ax.xaxis.label for ax in fig.axes} | {ax.yaxis.label for ax in fig.axes}
    boxes = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible() and t not in ticks and t not in labs]
    over = [(a.get_text(), b.get_text()) for i, (a, ba) in enumerate(boxes) for b, bb in boxes[i + 1:] if ba.overlaps(bb)]
    legend_on_axes = any(l.get_window_extent(r).overlaps(ax.get_window_extent(r)) for l in fig.legends for ax in fig.axes)
    return over, legend_on_axes


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("with_exp", [True, False])
def test_figures_have_no_overlapping_text(name, with_exp):
    c = STRUCTS[name]
    ex, em = _exp(c) if with_exp else (None, None)
    s = eng.run_boutinaud_2023_calculation(c, custom_exp_ex_nm=ex, custom_exp_em_nm=em)["summary_df"].drop(columns=["Raw Site"])
    if not s["Viable Bi3+ Site"].any():
        pytest.skip("no viable site")
    for fig in (bi3_plots.make_energy_diagram(s, c.formula, ex, em), bi3_plots.make_power_law_plot(s, ex)):
        assert _text_overlaps(fig) == ([], False)


def test_energy_axis_shows_every_level():
    """Version 1.0 cut the C-band of CaF2 (8.31 eV) off at 8 eV."""
    c = STRUCTS["CaF2 (ICSD-656448, Boutinaud Table 4 & 5)"]
    s = eng.run_boutinaud_2023_calculation(c)["summary_df"].drop(columns=["Raw Site"])
    fig = bi3_plots.make_energy_diagram(s, c.formula)
    assert fig.axes[0].get_ylim()[1] > float(s["E_C (eV) [Eq. 19]"].max())


@pytest.mark.parametrize("name", ["CaYAlO4 (ICDD 01-081-0742, I4/mmm)", "CaF2 (ICSD-656448, Boutinaud Table 4 & 5)",
                                  "CaTiO3 (ICSD-165801, Boutinaud Table 3, 5 & 9)"])
def test_export_zip_is_complete(name):
    c = STRUCTS[name]
    res = eng.run_boutinaud_2023_calculation(c, custom_exp_ex_nm=280.0, custom_exp_em_nm=398.0)
    sp = res["summary_df"].drop(columns=["Raw Site"])
    data = bi3_export.build_zip(sp, res["matrix_I_df"], res["subshells_df"], c.formula, c.space_group,
                                "V = %.2f A^3, Z = %d" % (c.V, c.Z), c.reference_note, 280.0, 398.0, "2.0.0",
                                "A. Author, Title, Journal 1, 1-2 (2026)")
    z = zipfile.ZipFile(io.BytesIO(data))
    base = bi3_export.safe_name(c.formula)
    expected = {base + s for s in ("_report.txt", "_optical_summary.csv", "_matrix_I_BVS.csv", "_subshells_he.csv",
                                   "_diagram_data.csv", "_energy_levels.png", "_energy_levels.svg", "_energy_levels.pdf",
                                   "_EA_vs_he.png", "_EA_vs_he.svg", "_EA_vs_he.pdf")}
    assert set(z.namelist()) == expected
    assert all(z.getinfo(n).file_size > 0 for n in z.namelist())
    assert z.read(base + "_energy_levels.png")[:4] == b"\x89PNG"
    assert z.read(base + "_energy_levels.pdf")[:4] == b"%PDF"
    assert pd.read_csv(io.BytesIO(z.read(base + "_optical_summary.csv"))).shape == sp.shape
    report = z.read(base + "_report.txt").decode("utf-8")
    assert "A. Author, Title, Journal 1, 1-2 (2026)" in report and "Boutinaud" not in report
    for csv in ("_optical_summary.csv", "_matrix_I_BVS.csv", "_subshells_he.csv"):
        assert "[Eq." not in z.read(base + csv).decode("utf-8").splitlines()[0]
    from PIL import Image
    dpi = Image.open(io.BytesIO(z.read(base + "_energy_levels.png"))).info["dpi"]
    assert round(dpi[0]) == 600 and round(dpi[1]) == 600


def test_json_with_missing_cell_parameters_is_detected():
    ins = eng.inspect_json_for_missing_parameters(json.loads(JSON_DEMO.read_text(encoding="utf-8")))
    assert ins["is_complete"] is False and len(ins["missing_fields"]) == 4
    assert bi3_validate.check_json(json.loads(JSON_DEMO.read_text(encoding="utf-8"))) == []


CIF_OK = CIF_DEMO.read_text(encoding="utf-8")


@pytest.mark.parametrize("text", ["", "this is not a cif file", "data_x\n_cell_length_a 3.0\n",
                                  CIF_OK.replace("_cell_length_a", "_xx_a").replace("_cell_length_b", "_xx_b").replace("_cell_length_c", "_xx_c"),
                                  CIF_OK.replace("_atom_site_fract_x", "_atom_site_xx")])
def test_bad_cif_text_is_rejected(text):
    errors, _ = bi3_validate.check_cif_text(text)
    assert errors


def test_truncated_cif_is_rejected_after_parsing():
    text = CIF_OK[:600]
    errors, _ = bi3_validate.check_cif_text(text)
    if not errors:
        assert bi3_validate.check_crystal(eng.parse_cif_string(text, filename="t.cif"))


def test_good_cif_is_accepted_and_missing_z_warns():
    errors, warns = bi3_validate.check_cif_text(CIF_OK)
    assert errors == [] and warns == []
    assert bi3_validate.check_crystal(eng.parse_cif_string(CIF_OK, filename="ok.cif")) == []
    errors, warns = bi3_validate.check_cif_text(CIF_OK.replace("_cell_formula_units_Z", "_xx_Z"))
    assert errors == [] and any("Z = 1" in w for w in warns)


@pytest.mark.parametrize("raw", [[], 5, {}, {"a": 1}, {"crystal_structure": {}}])
def test_bad_json_is_rejected(raw):
    assert bi3_validate.check_json(raw)


def test_every_module_is_listed_in_pyproject():
    listed = set(re.findall(r'"(bi3_\w+)"', (ROOT / "pyproject.toml").read_text(encoding="utf-8")))
    used = set()
    for f in ROOT.glob("bi3_*.py"):
        used |= set(re.findall(r"^\s*import (bi3_\w+)", f.read_text(encoding="utf-8"), re.M))
        used.add(f.stem)
    assert used <= listed, used - listed


def test_requirements_do_not_contain_pymatgen():
    assert "pymatgen" not in (ROOT / "requirements.txt").read_text(encoding="utf-8").lower()


# ---------------------------------------------------------------------------- the application itself
@pytest.fixture(scope="module")
def app_modules():
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    return AppTest


FAKE_CITATION = {"authors": "A. Author", "title": "A Title", "journal": "J. Test", "volume": "7", "pages": "1-9",
                 "year": "2026", "doi": "10.1000/xyz", "note": ""}


@pytest.fixture
def fake_citation(monkeypatch):
    import bi3_format
    monkeypatch.setattr(bi3_format, "fetch_citation", lambda *a, **k: (dict(FAKE_CITATION), ""))


def _all_text(at):
    parts = []
    for coll in (at.markdown, at.caption, at.subheader, at.title, at.header, at.warning, at.info, at.success, at.code,
                 at.sidebar.markdown, at.sidebar.caption, at.sidebar.header, at.sidebar.info, at.sidebar.success):
        parts += [str(e.value) for e in coll]
    parts += [str(t.label) for t in at.tabs]
    parts += [str(o) for s in at.sidebar.selectbox for o in s.options]
    return "\n".join(parts)


def test_app_runs_in_every_input_mode(app_modules, fake_citation):
    AppTest = app_modules
    at = AppTest.from_file(str(ROOT / "bi3_app.py"), default_timeout=300).run()
    assert not at.exception                                              # default: CIF mode with the demo structure
    at.sidebar.radio[0].set_value("Upload .json File (MP / Custom)").run()
    assert not at.exception                                              # JSON mode with the demo structure
    at.sidebar.radio[0].set_value("Built-In Benchmarks").run()
    import bi3_format
    for key in eng.get_boutinaud_2023_benchmarks():
        at.sidebar.selectbox[0].set_value(key).run()
        assert not at.exception, key
        c = STRUCTS[key]
        ex = round(HC / c.exp_EA_eV, 1) if c.exp_EA_eV else 280.0
        em = round(HC / (c.exp_EA_eV - c.exp_S_eV), 1) if (c.exp_EA_eV and c.exp_S_eV) else 398.0
        ref = eng.run_boutinaud_2023_calculation(c, custom_exp_ex_nm=ex, custom_exp_em_nm=em)["summary_df"]
        ref = ref[ref["Viable Bi3+ Site"]].drop(columns=["Raw Site", "Viable Bi3+ Site"]).reset_index(drop=True)
        shown = [m.value for m in at.markdown]
        assert bi3_format.table_html(ref) in shown, key                  # the numbers of the first table
        text = _all_text(at)
        assert bi3_format.formula_html(c.formula) in text, key           # formula with real subscripts
        assert bi3_format.space_group_html(c.space_group) in text, key
        assert bi3_format.site_count_text(len(ref)) in text, key
        assert "Boutinaud" not in text and "Eq. " not in text and "Eqs." not in text, key
    at.sidebar.checkbox[0].set_value(False).run()                        # experimental peaks off
    assert not at.exception


def test_help_and_citation_come_from_the_live_record(app_modules, fake_citation):
    at = app_modules.from_file(str(ROOT / "bi3_app.py"), default_timeout=300).run()
    text = _all_text(at)
    assert bi3_format_citation() in text
    assert "Boutinaud" not in text


def bi3_format_citation():
    import bi3_format
    return bi3_format.citation_text(FAKE_CITATION)


def test_offline_citation_is_reported_not_invented(app_modules, monkeypatch):
    import bi3_format
    monkeypatch.setattr(bi3_format, "fetch_citation", lambda *a, **k: (None, "offline"))
    at = app_modules.from_file(str(ROOT / "bi3_app.py"), default_timeout=300).run()
    assert not at.exception
    text = _all_text(at)
    assert "could not be loaded" in text and bi3_format.CITATION_URL in text


def test_citation_is_downloaded_from_the_web_and_never_from_an_old_copy(monkeypatch):
    import bi3_format
    calls = []

    class Resp:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self): return json.dumps(FAKE_CITATION).encode("utf-8")

    def fake_open(req, timeout=0):
        calls.append(req.full_url)
        return Resp()

    monkeypatch.setattr(bi3_format.urllib.request, "urlopen", fake_open)
    bi3_format._cache.update(t=0.0, data=None, fail_t=0.0, err="")
    data, err = bi3_format.fetch_citation(force=True)
    assert calls == ["https://omasgroup.org/citation2.json"] and data["title"] == "A Title" and err == ""
    assert bi3_format.citation_text(data) == "A. Author, A Title, J. Test 7, 1-9 (2026), https://doi.org/10.1000/xyz"

    def broken(req, timeout=0):
        raise OSError("no network")

    monkeypatch.setattr(bi3_format.urllib.request, "urlopen", broken)
    data, err = bi3_format.fetch_citation(force=True)
    assert data is None and "no network" in err                           # no stale copy after a failure
    bi3_format._cache.update(t=0.0, data=None, fail_t=0.0, err="")


def test_formulas_are_typeset_with_subscripts():
    f = bi3_format
    assert f.formula_html("Y2O3") == "Y<sub>2</sub>O<sub>3</sub>"
    assert f.formula_unicode("Cs2NaYCl6") == "Cs₂NaYCl₆"
    assert f.formula_unicode("CaYAlO4") == "CaYAlO₄"
    assert f.formula_unicode("Ca1 Y1 Al1 O4") == "CaYAlO₄"               # sum formula: counts of 1 are dropped
    assert f.formula_unicode("Ca(OH)2") == "Ca(OH)₂"
    assert f.formula_unicode("CuSO4·5H2O") == "CuSO₄·5H₂O"                # coefficient stays at normal height
    assert f.formula_math("Y3Al5O12") == "Y$_{3}$Al$_{5}$O$_{12}$"
    assert f.compound_unicode("Y2O3 (ICSD-192862, S6 & C2 Sites, Boutinaud Table 5)") == "Y₂O₃ (ICSD-192862, S6 & C2 Sites)"
    assert f.compound_unicode("CaTiO3 (ICSD-165801, Boutinaud Table 3, 5 & 9)") == "CaTiO₃ (ICSD-165801)"


def test_space_groups_are_typeset():
    f = bi3_format
    assert f.space_group_unicode("Fm-3m (No. 225)") == "Fm3̄m (No. 225)"
    assert f.space_group_unicode("I41/amd (No. 141)") == "I4₁/amd (No. 141)"
    assert f.space_group_unicode("P 21/c") == "P2₁/c"
    assert f.space_group_unicode("P 4_2/m n m") == "P4₂/mnm"
    assert f.space_group_unicode("P63/mmc") == "P6₃/mmc"
    assert f.space_group_unicode("I4/mmm (No. 139)") == "I4/mmm (No. 139)"
    assert f.space_group_html("Ia-3d (No. 230)") == '<i>Ia</i><span style="text-decoration:overline">3</span><i>d</i> (No. 230)'


def test_site_counts_are_written_in_words():
    assert [bi3_format.site_count_text(n) for n in (0, 1, 2, 3, 4)] == ["No sites", "1 site", "2 sites", "3 sites", "4 sites"]


def test_table_headers_are_typeset_and_have_no_equation_numbers():
    res = eng.run_boutinaud_2023_calculation(STRUCTS["CaF2 (ICSD-656448, Boutinaud Table 4 & 5)"])
    for key in ("summary_df", "matrix_I_df", "subshells_df"):
        html = bi3_format.table_html(res[key])
        assert "[Eq." not in html and "Eq. " not in html and "[Table" not in html
        assert not any("[Eq." in c for c in bi3_format.plain_columns(res[key]).columns)
    html = bi3_format.table_html(res["summary_df"])
    assert "<i>E</i><sub>A</sub> (eV)" in html and "Δ<i>E</i>(Bi<sup>3+</sup>) (eV)" in html


def test_zip_report_names_interpolated_cn11_parameters():
    """A site with CN = 11 (not in Table 2 of the paper) is flagged in the table and in the ZIP report."""
    import copy
    c = copy.deepcopy(eng.get_boutinaud_2023_benchmarks()["CaF2 (ICSD-656448, Boutinaud Table 4 & 5)"])
    c.sites_cations["Ca"]["N"] = 11
    sdf = eng.run_boutinaud_2023_calculation(c)["summary_df"]
    assert int(sdf["CN"].iloc[0]) == 11
    sdf["CT parameters"] = sdf["CN"].map(bi3_export.ct_parameter_note)
    res = eng.run_boutinaud_2023_calculation(c)
    z = zipfile.ZipFile(io.BytesIO(bi3_export.build_zip(
        sdf.drop(columns=["Raw Site"]), res["matrix_I_df"], res["subshells_df"], "CaF2", "Fm-3m", "V = 1", "", 280.0, 398.0, "2.0.0", "A. Author, Title")))
    report = z.read([n for n in z.namelist() if n.endswith("_report.txt")][0]).decode("utf-8")
    assert "interpolated" in report


def test_display_mode_by_platform(monkeypatch):
    """Own window on Windows (with WebView2) and macOS; web browser on Linux; the options override."""
    import bi3_desktop as d
    monkeypatch.setattr(d, "webview2_available", lambda: True)
    assert d.use_browser_mode([], "win32") is False and d.use_browser_mode([], "darwin") is False
    assert d.use_browser_mode([], "linux") is True
    assert d.use_browser_mode(["--window"], "linux") is False and d.use_browser_mode(["--browser"], "darwin") is True
    monkeypatch.setattr(d, "webview2_available", lambda: False)
    assert d.use_browser_mode([], "win32") is True


def test_start_files_for_every_system_exist_and_use_unix_line_endings():
    for name in ("start_calculator.sh", "Start_Calculator.command", "build_app.sh"):
        raw = (ROOT / name).read_bytes()
        assert raw.startswith(b"#!/usr/bin/env bash") and b"\r" not in raw, name
    assert (ROOT / "Start_Calculator.bat").exists() and (ROOT / "build_app.py").exists()
    assert (ROOT / "bi3_icon.icns").exists()
