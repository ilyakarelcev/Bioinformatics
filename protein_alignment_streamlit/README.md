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

## Run

```powershell
cd C:\Users\53k\Desktop\BioSeq\protein_alignment_streamlit
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

If `python` is not in PATH, run the same commands from the Python environment you normally use.

The app can still run without Biopython because it includes a small built-in BLOSUM62 aligner. Installing
the requirements switches the calculation layer to Biopython's pairwise alignment implementation.
