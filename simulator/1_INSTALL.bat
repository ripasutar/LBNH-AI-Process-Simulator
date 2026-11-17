@echo off
cd /d "%~dp0"
echo Installing the Python packages (needs internet, about 2 minutes)...
python -m pip install -r requirements.txt
pause
