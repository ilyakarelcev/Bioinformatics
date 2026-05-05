@echo off
setlocal

cd /d "%~dp0protein_alignment_streamlit"

echo Starting Protein Pairwise Alignment Lab...
echo.
echo If this is the first launch on this computer, run:
echo   python -m pip install -r requirements.txt
echo.
echo App URL:
echo   http://localhost:8502
echo.

python -m streamlit run app.py --server.port 8502 --browser.gatherUsageStats false

echo.
echo Streamlit stopped. Press any key to close this window.
pause >nul
