@echo off
setlocal EnableExtensions
cd /d "%~dp0"
title Bi3+ Luminescence-Structure Calculator

rem ---------------------------------------------------------------
rem  One double-click: the first start sets everything up (one time,
rem  needs internet); every later start opens the calculator at once.
rem  Start_Calculator.bat check     set up, run the self-test, stop
rem  set BI3_PYTHON=C:\path\python.exe   to choose a particular Python
rem ---------------------------------------------------------------

set "VENV=%~dp0.venv"
set "PYEXE=%VENV%\Scripts\python.exe"

if exist "%PYEXE%" goto have_venv

echo.
echo  First start: setting up the calculator. This needs internet and takes
echo  about 2 to 5 minutes. It is done only once.
echo.
call :find_python
if not defined BASEPY goto no_python
echo  Using Python: %BASEPY%
%BASEPY% -m venv "%VENV%"
if errorlevel 1 goto fail_venv

:have_venv
fc /b requirements.txt "%VENV%\requirements.installed.txt" >nul 2>&1
if not errorlevel 1 goto run
echo  Installing components (streamlit, numpy, pandas, matplotlib) ...
"%PYEXE%" -m pip install --disable-pip-version-check -r requirements.txt
if errorlevel 1 goto fail_pip
copy /y requirements.txt "%VENV%\requirements.installed.txt" >nul
echo.
"%PYEXE%" bi3_launcher.py --check
if errorlevel 1 goto fail_check
echo.

:run
if /i "%~1"=="check" goto selftest_only
"%PYEXE%" bi3_launcher.py
if errorlevel 1 goto fail_run
exit /b 0

:selftest_only
"%PYEXE%" bi3_launcher.py --check
exit /b %errorlevel%

:find_python
set "BASEPY="
if defined BI3_PYTHON "%BI3_PYTHON%" -c "import sys" >nul 2>&1 && set "BASEPY="%BI3_PYTHON%""
if defined BASEPY exit /b 0
for %%V in (3.12 3.13 3.11 3.10 3.9) do (
    if not defined BASEPY (
        py -%%V -c "import sys" >nul 2>&1 && set "BASEPY=py -%%V"
    )
)
if defined BASEPY exit /b 0
python -c "import sys; sys.exit(0 if sys.version_info[:2]>=(3,9) else 1)" >nul 2>&1 && set "BASEPY=python"
exit /b 0

:no_python
echo.
echo  Python 3.9 or newer was not found on this computer.
echo.
echo   Option 1 (automatic): install Python now with the Windows Package Manager.
echo   Option 2 (manual):    download Python from https://www.python.org/downloads/windows/
echo                         and tick "Add python.exe to PATH" on the first screen.
echo.
where winget >nul 2>&1
if errorlevel 1 goto no_winget
choice /c YN /m "Install Python 3.12 automatically now"
if errorlevel 2 goto end_fail
winget install --id Python.Python.3.12 -e --source winget --accept-package-agreements --accept-source-agreements
echo.
echo  Python has been installed. Please CLOSE this window and double-click
echo  Start_Calculator.bat again.
goto end_fail

:no_winget
echo  The Windows Package Manager is not available here. Please use option 2,
echo  then double-click Start_Calculator.bat again.
goto end_fail

:fail_venv
echo.
echo  The virtual environment could not be created in:
echo    %VENV%
echo  Check that the folder is writable (copy the whole folder to your Documents
echo  folder if it sits in a protected location) and try again.
goto end_fail

:fail_pip
echo.
echo  The components could not be installed. Please check the internet
echo  connection (a proxy or firewall may block pip) and double-click again.
echo  If the problem stays, delete the folder ".venv" and start again.
goto end_fail

:fail_check
echo.
echo  The self-test failed. Please send the lines above to the authors.
goto end_fail

:fail_run
echo.
echo  The calculator stopped with an error. Please send the lines above to the authors.
goto end_fail

:end_fail
echo.
pause
exit /b 1
