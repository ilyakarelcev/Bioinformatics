# Protein Pairwise Alignment Lab

Streamlit app for comparing two amino-acid sequences with pairwise alignment.

## Features

- Smith-Waterman local alignment for fragments, peptides, domains, and partial proteins.
- Needleman-Wunsch global alignment for full-length protein comparison.
- BLOSUM62 scoring with editable gap open and gap extension penalties.
- Alignment metrics: score, identity, similarity, gaps, coverage, aligned region.
- Tile-based visual alignment viewer with match/similarity/gap highlighting.
- Optional residue coloring by biochemical class.
- Downloadable text report.
- Dark theme by default, with a Light/Dark switch in the sidebar.

## Run

Double-click the launcher in the parent folder:

```text
C:\Users\53k\Desktop\BioSeq\run_alignment_app.bat
```

Then open:

```text
http://localhost:8502
```

If dependencies are missing, double-click first:

```text
C:\Users\53k\Desktop\BioSeq\install_alignment_app_requirements.bat
```

Manual launch:

```powershell
cd C:\Users\53k\Desktop\BioSeq\protein_alignment_streamlit
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

If `python` is not in PATH, run the same commands from the Python environment you normally use.

The app can still run without Biopython because it includes a small built-in BLOSUM62 aligner. Installing
the requirements switches the calculation layer to Biopython's pairwise alignment implementation.
