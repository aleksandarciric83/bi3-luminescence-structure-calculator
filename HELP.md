### Quick guide

**1. Start the program.** Windows: open the folder `Bi3_Calculator` and double-click `Bi3_Calculator.exe`. The program has everything it needs inside the folder: Python, all libraries and the web components are included, and nothing has to be downloaded or installed. It works without an internet connection; only the citation at the top and in tab 6 is loaded live from the web (see "How to cite"). Keep the whole folder together (copy or zip the folder, not the `.exe` alone). To close the program, close its window. Windows may show a "protected your PC" notice the first time because the program is not digitally signed: choose *More info*, then *Run anyway*.

If your Windows is old and lacks the Microsoft Edge WebView2 component, the program opens the calculator in your web browser instead and shows a small window with a *Stop* button.

**macOS and Linux.** Download the program for your system (`Bi3_Calculator_v2.0.0_macos.zip` or `Bi3_Calculator_v2.0.0_linux.zip` from the GitHub release page), unzip it and start `Bi3_Calculator.app` (macOS: the first time right-click, then *Open*) or the file `Bi3_Calculator` (Linux: in a terminal, `./Bi3_Calculator`). On macOS the calculator opens in its own window; on Linux it opens in your default web browser, and you stop it with Ctrl+C in the terminal.

The source version (for developers) starts with `Start_Calculator.bat` (Windows), `Start_Calculator.command` (macOS, double-click) or `./start_calculator.sh` (macOS and Linux). It installs four Python libraries on the first start and needs internet for that one time; on Linux it needs Python 3.9 or newer (if the `venv` module is missing, the script fetches pip itself).

**2. Choose the structure** in the sidebar.
- *Upload .cif*: any CIF with cell lengths, symmetry operations and atom sites. A short message appears if something essential is missing, for example the cell lengths. If `_cell_formula_units_Z` is missing, Z = 1 is assumed and a warning is shown; check it, because Z changes the stoichiometric coefficients.
- *Upload .json*: a Materials Project style file. If the unit cell is not in the file, the program estimates it from the density and asks you to confirm or edit a, b, c, the angles and Z.
- *Built-in benchmarks*: 11 compounds (13 sites) with published values, plus CaYAlO₄ in two structure models, useful for a first test.

**3. Experimental peaks (optional).** Enter the maximum of the lowest-energy (A-band) excitation band and of the emission. They give the Stokes shift $S$ that is needed for the emission-origin criterion. $S$ is the separation of the excitation and emission maxima; the absorbing state need not be the emitting state. The program does not check which band the number belongs to.

**4. Read the tabs.**
1. *Optical transitions*: for every viable site, the nephelauxetic functions $h_e(\mathrm{M})$ and $h_e(\mathrm{Bi}^{3+})$, the A-band energy $E_A$ and C-band energy $E_C$, the charge-transfer energies $E_\mathrm{MMCT}$ and $E_\mathrm{IVCT}$, and the criterion $\Delta E(\mathrm{Bi}^{3+})$. $\Delta E$ below 1.39 eV points to A-band emission, above 1.41 eV to charge-transfer emission. The formulas are printed above the table.
2. *Matrices I and II*: coordination numbers, decomposition coefficients and bond valence sums.
3. *Dielectric and nephelauxetic*: the step-by-step sub-shell table, with all formulas printed above it.
4. *Energy diagrams and export*: the two diagrams and all downloads.
5. *Benchmark suite*: all built-in compounds, calculated against the published values.

**5. Export.** Under each diagram: PNG (600 dpi), SVG and PDF. Below them: one ZIP with both diagrams in all three formats, the tables as CSV, the numbers behind the diagrams and a short report that records the program version, the inputs and how to cite.

### The formulas
- A band: $E_A(\mathrm{Bi}^{3+}) = 9.41\,\left(1 - 0.555\,h_e(\mathrm{Bi}^{3+})^{0.240}\right)$, accuracy about 0.35 eV.
- C band: $E_C(\mathrm{Bi}^{3+}) = 3.236 + 2.29\,\left(E_A(\mathrm{Bi}^{3+}) - 2.972\right)^{0.856}$.
- Charge transfer: $E_\mathrm{CT}(\mathrm{Bi}^{3+}) = k\,\left[\chi(\mathrm{Bi}^{3+}) - \alpha\,\chi(\mathrm{acceptor})/d\right]$, with $k$ and $\alpha$ tabulated for each pair of coordination numbers.
- Emission origin: $\Delta E(\mathrm{Bi}^{3+}) = \left|E_\mathrm{CT}(\mathrm{Bi}^{3+}) - E_A(\mathrm{Bi}^{3+}) - S\right|$.

The formulas of the sub-shell calculation ($K$, $Q$, $N^*$, $v$, $\Gamma_1$, $\Gamma_2$, $\Gamma_3$, $f_\mathrm{c}$, $\alpha$, $h_e$) are printed in tab 3.

### What to keep in mind
- The method is empirical: it rests on relations fitted to 177 compounds (242 sites). The formula for $E_A$ has a stated accuracy of about 0.35 eV. The results indicate trends and candidate sites; they do not replace a measurement.
- A *viable* site is a candidate (coordination number 5 or more, and the cation is not a network former such as P, B, Si, S, V, Mo, W or As). It is not a prediction that Bi³⁺ occupies it. The number of viable sites is shown at the top as "No sites", "1 site", "2 sites" and so on.
- The experimental excitation maximum stored with a benchmark is not always the A band. The benchmark tab gives the band type of each value: for YPO₄ it is the C band (5.46–5.51 eV; the A band is not observed), and for CaTiO₃ the charge-transfer band (3.65–3.87 eV, $S$ = 1.46–1.64 eV; the A-type value is 3.35 eV, $S$ = 1.13 eV). Where a range is published, the stored number lies inside it and is marked.
- Charge-transfer energies use electronegativities (Li and Xue 2006, J. Phys. Chem. A 110, 11332) and tabulated parameters for the coordination numbers CN = 5, 6, 7, 8, 9, 10 and 12. The parameters for CN = 11 are a linear interpolation between CN = 10 and CN = 12 and are marked in the results; sites with CN > 12 use the CN = 12 parameters; an acceptor with CN′ other than 4 uses the columns for CN′ ≥ 6.
- The experimental peaks are drawn at every site of the energy diagram for comparison. They are properties of the compound, not of one site.
- Different structural models of the same compound give slightly different numbers, for example CaYAlO₄ described in I4/mmm or in I4mm.

### If something does not work
- *The program does not open*: make sure the whole `Bi3_Calculator` folder was unzipped (not started from inside the zip file), and that antivirus software did not block it.
- *A CIF is refused*: read the message in the sidebar; it names the missing item.
- *The page stays blank*: wait a few seconds; in the browser version reload the tab.
- *The citation is not shown*: there is no internet connection or the website did not answer; the rest of the program works normally.
- Source version only: *the installation fails* — check the internet connection (a company proxy may block pip), delete the folder `.venv` and start again. Keep the folder path short (for example `C:\Bi3`), because Windows cannot load some libraries from very long paths.
