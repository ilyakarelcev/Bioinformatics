@echo off
setlocal

cd /d "%~dp0protein_alignment_streamlit"

echo Installing Protein Pairwise Alignment Lab requirements...
echo.

python -m pip install -r requirements.txt

echo.
echo Done. Press any key to close this window.
pause >nul
