# UniProt Streamlit Lab

Small local Streamlit prototype for exploring UniProt REST API responses.

## Run

```powershell
cd C:\Users\53k\Desktop\BioSeq\uniprot_streamlit
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

## What it does

- Fetches UniProtKB records by accession, for example `P01308`.
- Shows JSON, FASTA, and compact TSV-style fields.
- Searches UniProtKB with query strings.
- Hydrates a list of similar protein hits returned by an embedding/graph search.
