@echo off
cd /d "%~dp0"
echo Starting the LBNH live simulator. Keep this window open; close it to stop.
python live_app.py
pause
