# How to run the Protein Alignment app

## Normal launch

Double-click:

```text
run_alignment_app.bat
```

Then open:

```text
http://localhost:8502
```

Keep the black terminal window open while you use the app. Closing it stops Streamlit.

The app starts in dark theme by default. Use the Theme switch in the left sidebar to change it to Light.

## First launch or after reinstalling Python

Double-click:

```text
install_alignment_app_requirements.bat
```

Then double-click:

```text
run_alignment_app.bat
```

## Manual launch from terminal

```powershell
cd C:\Users\53k\Desktop\BioSeq\protein_alignment_streamlit
python -m pip install -r requirements.txt
python -m streamlit run app.py --server.port 8502
```

The app file is:

```text
C:\Users\53k\Desktop\BioSeq\protein_alignment_streamlit\app.py
```
