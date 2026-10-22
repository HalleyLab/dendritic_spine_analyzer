@echo off
setlocal
set "PYTHONDONTWRITEBYTECODE=1"
set "TEMP=%~dp0.notebook_runtime\tmp"
set "TMP=%TEMP%"
set "MPLCONFIGDIR=%~dp0.notebook_runtime\matplotlib"
set "PIP_CACHE_DIR=%~dp0.notebook_runtime\pip-cache"
if not exist "%TEMP%" mkdir "%TEMP%"
cd /d "%~dp0"
set "BUNDLED_PY=%USERPROFILE%\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"
if exist "%BUNDLED_PY%" (
    "%BUNDLED_PY%" -B app.py %*
    if errorlevel 1 pause
    exit /b
)
if not exist ".venv\Scripts\python.exe" (
    py -3.12 -m venv .venv
    if errorlevel 1 (
        echo Please install Python 3.12 64-bit with tkinter, then try again.
        pause
        exit /b 1
    )
)
".venv\Scripts\python.exe" -c "import numpy,scipy,matplotlib,skimage,tifffile,czifile,imagecodecs" >nul 2>&1
if errorlevel 1 (
    ".venv\Scripts\python.exe" -m pip install -r requirements.txt
    if errorlevel 1 (
        echo Dependency installation failed. See IMAGEJ_GUIDE.md for setup.
        pause
        exit /b 1
    )
)
".venv\Scripts\python.exe" -B app.py %*
if errorlevel 1 pause
