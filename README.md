# Bi3+ Luminescence-Structure Calculator (version 2.0.0)

A program that calculates the optical properties of Bi3+ in inorganic solids directly from the crystal structure: nephelauxetic functions, the A- and C-band energies, metal-to-metal and intervalence charge-transfer energies, and the emission-origin criterion.

Input is a CIF file, a Materials Project style JSON file, or one of the built-in benchmarks (11 compounds with 13 sites with published values, plus CaYAlO4 in two structure models). Output is a set of tables, two diagrams (PNG at 600 dpi, SVG, PDF) and a short report. The reference to cite is loaded live from https://omasgroup.org/citation2.json each time it is shown.

## Start

Windows (recommended): unzip `Bi3_Calculator_v2.0.0_windows.zip` and double-click **Bi3_Calculator.exe** in the folder `Bi3_Calculator`. Everything is inside the folder (Python, all libraries, the web components); nothing has to be installed or downloaded, and no internet is needed except to show the citation. Keep the whole folder together.

Windows, from the source files (developers): double-click **Start_Calculator.bat**. The first start installs what is needed (internet required, 2 to 5 minutes, once). If Python is missing, the file offers to install it. Keep the folder path short (for example `C:\Bi3`); very long paths make Windows fail to load libraries.

macOS and Linux (programs): download `Bi3_Calculator_v2.0.0_macos.zip` or `Bi3_Calculator_v2.0.0_linux.zip` from the Releases page, unzip, and start `Bi3_Calculator.app` (macOS; first time right-click, Open; own window) or `./Bi3_Calculator` (Linux; opens in the default web browser, stop with Ctrl+C).

macOS and Linux (from the source files): `./start_calculator.sh` (macOS also: double-click `Start_Calculator.command`). Needs Python 3.9 or newer; the first start installs four libraries (internet, 2 to 5 minutes, once). Tested from scratch on Ubuntu with Python 3.14 (also when the `venv` module is missing); not yet run on a Mac.

Building a standalone program yourself: `python build_app.py` (Windows, macOS, Linux; run it on the system the program is for) or `build_exe.bat` (Windows). The result is `dist/Bi3_Calculator` with all libraries inside; `Bi3_Calculator --selftest` checks that it can write all export files. If the Microsoft Edge WebView2 component is missing on a very old Windows, the program opens the calculator in the web browser instead and shows a small window with a Stop button. The GitHub workflow in `.github/workflows/build.yml` builds the three programs and runs the tests on all three systems.

With pip: `pip install .` and then `bi3-calculator`.

Python 3.9 or newer is required. The libraries are streamlit, numpy, pandas and matplotlib (see requirements.txt). The program runs only on your own computer and sends nothing anywhere.

## Contents

| File | Purpose |
|---|---|
| bi3_app.py | the web interface |
| bi3_engine.py | the calculation (formulas and parameters unchanged from version 1.0; one stored published value corrected in 1.1.1) |
| bi3_plots.py | the two diagrams |
| bi3_format.py | typesetting (formulas with subscripts, space-group symbols, table headers) and the live citation |
| bi3_icon.ico, bi3_icon.png | the program icon |
| bi3_export.py | figure files, tables, report and ZIP |
| bi3_validate.py | checks of uploaded files |
| bi3_launcher.py | starts the program in the web browser; `--check` runs a self-test |
| bi3_desktop.py, build_exe.bat | the same program in its own window; build of the .exe |
| HELP.md | the help text shown in the program |
| examples/ | CaYAlO4 as CIF and as JSON |
| tests/ | 83 automatic tests (`python -m pytest tests`) |
| make_screenshots.py | takes screenshots of every tab (needs `pip install playwright`) |

## Changes from version 1.0

See CHANGELOG.md. Version 1.2.0 changed the presentation only (live citation, 600 dpi PNG, typeset formulas, space groups and table headers, full equations instead of equation numbers, program icon, self-contained program). The calculation is identical; the tests compare every number with the frozen version 1.0 results. Version 1.1.1 corrects one stored published value (Y2O3, C2 site) and labels the experimental values of the benchmarks. Charge-transfer parameters for CN = 11 are interpolated (the parameter set has no CN = 11 entry); the program marks such sites.

## How to cite

If results of this program are used, please cite the reference published at https://omasgroup.org/citation2.json (the program loads it from there).

## Licence

GNU General Public License, version 3 (GPL-3.0-only). See the file LICENSE. You may use, study, change and share the program under the terms of that licence; changed versions that you distribute must be offered under the same licence.
