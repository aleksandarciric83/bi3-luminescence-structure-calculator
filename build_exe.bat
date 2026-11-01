@echo off
setlocal EnableExtensions
cd /d "%~dp0"
rem Builds dist\Bi3_Calculator\Bi3_Calculator.exe: a program with its own window and ALL libraries inside
rem (Python, Streamlit, NumPy, pandas, Matplotlib, pywebview). Users need no Python, no pip and no download.
rem Needs internet once, on the PC where it is built. Run Start_Calculator.bat once first so that .venv exists.
if not exist ".venv\Scripts\python.exe" (
    echo  Please run Start_Calculator.bat once first.
    pause
    exit /b 1
)
".venv\Scripts\python.exe" -m pip install --disable-pip-version-check pywebview pyinstaller
if errorlevel 1 goto fail
".venv\Scripts\python.exe" -m PyInstaller --noconfirm --clean --windowed --name Bi3_Calculator --icon bi3_icon.ico ^
  --collect-all streamlit --copy-metadata streamlit --collect-submodules streamlit ^
  --collect-submodules matplotlib.backends --hidden-import matplotlib.backends.backend_svg ^
  --hidden-import matplotlib.backends.backend_pdf --hidden-import matplotlib.backends.backend_agg ^
  --hidden-import tkinter ^
  --add-data "bi3_app.py;." --add-data "HELP.md;." --add-data "examples;examples" ^
  --add-data "bi3_icon.ico;." --add-data "bi3_icon.png;." ^
  bi3_desktop.py
if errorlevel 1 goto fail
echo.
echo  Done: dist\Bi3_Calculator\Bi3_Calculator.exe
echo  Zip the whole folder dist\Bi3_Calculator and give it to the users.
pause
exit /b 0
:fail
echo.
echo  The build failed. Please send the lines above to the authors.
pause
exit /b 1
