@echo off
REM Launch the interactive supplement on Windows (double-click or run from a command prompt).
REM Creates a local virtual environment on first run, installs the requirements, starts Streamlit.
cd /d "%~dp0"
where python >nul 2>nul || (echo Python 3.9 or newer is required: https://www.python.org/downloads/ & pause & exit /b 1)
if not exist .venv\Scripts\python.exe (
    echo Creating virtual environment ...
    python -m venv .venv
)
echo Installing / checking requirements ...
.venv\Scripts\python.exe -m pip install --quiet --upgrade pip
.venv\Scripts\python.exe -m pip install --quiet -r requirements.txt
echo Starting the dashboard at http://localhost:8501  (close this window to stop)
.venv\Scripts\python.exe -m streamlit run dashboard.py --server.port 8501 --server.address localhost
pause
