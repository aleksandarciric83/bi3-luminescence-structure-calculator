"""Build the self-contained program (all libraries inside) with PyInstaller, on Windows, macOS or Linux.

    python build_app.py

Result:  dist/Bi3_Calculator/   (Windows: Bi3_Calculator.exe, Linux: Bi3_Calculator, macOS: also dist/Bi3_Calculator.app)
Run it with the Python that has the requirements installed (the .venv made by Start_Calculator.bat or
start_calculator.sh, or any environment with `pip install -r requirements.txt`). pywebview and PyInstaller
are installed here if they are missing. The build has to run on the system the program is meant for:
a Windows build does not run on macOS or Linux and the other way round.
"""
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEP = os.pathsep  # ';' on Windows, ':' elsewhere (the separator PyInstaller expects in --add-data)


def run(*cmd):
    print(">", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True, cwd=str(HERE))


def main():
    py = sys.executable
    run(py, "-m", "pip", "install", "--disable-pip-version-check", "pywebview", "pyinstaller")
    icon = "bi3_icon.ico" if sys.platform == "win32" else ("bi3_icon.icns" if sys.platform == "darwin" else "bi3_icon.png")
    cmd = [py, "-m", "PyInstaller", "--noconfirm", "--clean", "--windowed", "--name", "Bi3_Calculator", "--icon", icon,
           "--collect-all", "streamlit", "--copy-metadata", "streamlit", "--collect-submodules", "streamlit",
           "--collect-submodules", "matplotlib.backends",
           "--hidden-import", "matplotlib.backends.backend_svg", "--hidden-import", "matplotlib.backends.backend_pdf",
           "--hidden-import", "matplotlib.backends.backend_agg", "--hidden-import", "tkinter",
           "--add-data", "bi3_app.py" + SEP + ".", "--add-data", "HELP.md" + SEP + ".",
           "--add-data", "examples" + SEP + "examples",
           "--add-data", "bi3_icon.ico" + SEP + ".", "--add-data", "bi3_icon.png" + SEP + "."]
    if sys.platform == "darwin":
        cmd += ["--osx-bundle-identifier", "org.omasgroup.bi3calculator"]
    cmd.append("bi3_desktop.py")
    run(*cmd)
    print("\nDone: dist%sBi3_Calculator. Give the users the whole folder%s." %
          (os.sep, " (or dist/Bi3_Calculator.app on macOS)" if sys.platform == "darwin" else ""))


if __name__ == "__main__":
    main()
