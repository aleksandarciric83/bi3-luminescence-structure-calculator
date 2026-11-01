# Changelog

## 2.0.0 (2026-11-01)

Version 2: the 1.2.0 program (live citation, 600 dpi PNG, typeset formulas, full equations) now also runs on Linux and macOS. The calculation is unchanged; the engine file is identical to 1.2.0.

- Linux: start_calculator.sh works on Debian / Ubuntu even when the Python `venv` module or pip is missing (it fetches pip); tested from scratch on Ubuntu (Python 3.14) with all 83 tests passing. The packaged program opens the calculator in the default web browser (or prints the address when there is no display) and stops its server on Ctrl+C or `kill`.
- macOS: Start_Calculator.command (double-click); the packaged program Bi3_Calculator.app shows the calculator in its own window (WKWebView). Built by the GitHub workflow; not run on a Mac by the authors of this change.
- build_app.py (Windows, macOS, Linux) and build_app.sh build the self-contained program; build_exe.bat still builds the Windows one. New icon file bi3_icon.icns.
- `Bi3_Calculator --selftest [file]`: calculates one structure and writes every export file (PNG, SVG, PDF, ZIP); a packaged build that lacks a library fails here (this is how the missing Matplotlib writers of the first 1.1.1 build would have been caught).
- GitHub Actions workflow (.github/workflows/build.yml): runs the tests on Windows, macOS and Linux, builds the three programs, runs the self-test on each, and attaches them to a release when a tag `v*` is pushed.
- .gitattributes keeps LF line endings (the engine hash test reads the file bytes).
- 83 automatic tests.

## 1.2.0 (2026-11-01)

Presentation release. No formula, coefficient, parameter or number of the calculation is changed; the tests still compare every number with the frozen version 1.0 results.

### Interface
- Citation: the paper by the method's author is no longer cited in the program. "How to cite" (tab 6), the line under the title and the report in the ZIP show the reference loaded from https://omasgroup.org/citation2.json every time (a successful download is reused for one minute; offline, the program says so and gives the address; an old copy is never shown).
- Names: "Built-In Benchmarks", "5. Benchmark Suite", "1. Structure Input".
- Equations: no equation numbers anywhere (tables, CSV files, figures, help). The formulas are printed in full above the tables of tabs 1 to 3 and in the help; column headers that carried "[Eq. n]" lost the tag.
- Typesetting: chemical formulas with subscripts (Y2O3 -> Y₂O₃) in the cards, the compound list, the benchmark table, the figure titles and the report; space-group symbols in Hermann-Mauguin form (italic letters, bar over -3, subscript screw axes, e.g. I4₁/amd, Fm3̄m); table headers with italic symbols, subscripts, superscripts and Greek letters (E_A, h_e(Bi3+), Δd, Γ1, α ...); tables are now HTML tables (not sortable).
- "Viable Bi3+ doping sites" is written "No sites", "1 site", "2 sites", ...
- PNG export at 600 dpi (was 300).
- Reference notes of the built-in compounds no longer name a paper (engine notes only; the hash of bi3_engine.py changed, the numbers did not).

### Program
- New icon (bi3_icon.ico / bi3_icon.png) for the program, the window and the browser tab.
- The standalone program is rebuilt with every library inside (also the Matplotlib SVG/PDF/PNG writers, which were missing in the first 1.1.1 build and made the diagram downloads fail), and it opens the calculator in the web browser when the Edge WebView2 component is missing.
- New module bi3_format.py.

### Tests
81 automatic tests.

## 1.1.1 (2026-11-01)

Release after a check of the built-in benchmarks against Boutinaud, PCCP 2023, 25, 11027. No formula, coefficient or parameter of the calculation is changed.

### Data correction (no physics change)
- bi3_engine.py: the stored published he(Bi3+) of Y2O3, site Y2 (C2), was 1.72; Table 5 of the paper gives 1.75. The S6 site (1.60 / 1.72) was correct. Largest |calc - published| is now 0.008 for he(M) (CaF2, Y2O3 S6) and 0.011 for he(Bi3+) (Cs2NaYCl6), mean about 0.004; Y2O3 C2 is -0.006.
- Notes only: the CaTiO3 note gave a single Delta_E of 2.19 eV, which is not in the paper; it now quotes the range 2.10-2.28 eV (Table 5, CT type) and 2.28 eV (Table 9). The CaYAlO4 note no longer names an internal file. Comments: electronegativities (not optical ones) and the CN = 11 row of Table 2.
- All other stored published numbers (Tables 4, 5, 6, 9 and Table 2) were rechecked and agree.

### Interface
- Benchmark tab and a new CSV download: column "Exp. band type" (A, C or CT) and "Exp. value from" (single value or range of the paper). YPO4 5.48 eV is the C band, CaTiO3 3.76 eV the CT band; numbers unchanged. Help step 3 explains that S is the excitation-emission separation and that the absorbing state need not be the emitting state.
- Charge-transfer parameters for CN = 11 are interpolated (Table 2 lists CN 5-10 and 12 only). A viable site that uses them is marked in column "CT parameters", with a warning under the table and a line in the ZIP report. Documented: CN > 12 uses the CN = 12 row, an acceptor with CN' other than 4 uses the CN' >= 6 columns.
- Help and README: 11 compounds (13 sites) plus CaYAlO4 in two structure models; electronegativities wording.

### Installation and release
- New: bi3_desktop.py and build_exe.bat, a Windows program in its own window (no browser, no Python needed to run it).
- The launcher listens on localhost only. make_screenshots.py ends cleanly (stops the server tree), waits for the page, and takes full-height pictures (use --width 1800).
- start_calculator.sh (macOS and Linux) is not tested.

### Tests
74 tests (python -m pytest tests). The published-value test is computed live (he(M) 0.0105, he(Bi3+) 0.012, the message names the worst site). The engine hash test now pins version 1.1.1 (version 1.0 plus the one data value and notes); one cell of the frozen golden table is corrected and every other number still equals version 1.0 to 1e-9.

## 1.1.0 (2026-10-31)

### Calculation
Unchanged. bi3_engine.py is byte-identical to version 1.0 and the tests compare every number of 15 structures with the frozen version 1.0 results.

### Installation and start
- One double-click on Windows (Start_Calculator.bat) and one command on macOS and Linux (start_calculator.sh). The first start sets up a private environment and runs a self-test; later starts open at once. If Python is missing, the Windows file offers to install it.
- The unused pymatgen requirement is removed; only streamlit, numpy, pandas and matplotlib are installed, with upper version bounds.
- pip install . provides the command bi3-calculator.

### Interface
- Both diagrams are redrawn: real subscripts in the formula, no overlapping labels (leader lines where levels lie close), colour-blind-safe colours, the legend outside the plot, the energy axis taken from the data.
- Experimental maxima are labelled as compound properties; the viable-site rule is stated under the table; a Help tab and a how-to-cite box are added.
- Streamlit calls that are being removed (use_container_width) are replaced by a version-safe form. The benchmark table is calculated once and cached.

### Export
- PNG (300 dpi), SVG and PDF for both diagrams, a CSV with the numbers behind the diagrams, and one ZIP with figures, tables and a short report (program version, inputs, how to cite).

### Fixed
- Levels above 8 eV (C-band of CaF2 8.31 eV, SrF2 8.18 eV, BaF2 8.14 eV) are no longer cut off.
- A CIF without cell lengths no longer becomes a 5 Angstrom cube, a CIF without Z gives a warning, and a file without atom positions or sites is refused with a message instead of silent defaults or a crash.

### Tests
70 automatic tests (python -m pytest tests), documented in the test report.

## Open points for the authors
See HANDOVER.md.
