# XtalHop-Predictor

A Streamlit web tool that takes a protein amino acid sequence, searches the RCSB PDB for structural homologs via MMseqs2, extracts their crystallization conditions, and aggregates them into a consensus starting formulation. It also runs a local biophysical risk panel (isoelectric point, instability index, GRAVY score, amino acid composition flags) powered by Biopython.

**Live app:** https://xtalhop.streamlit.app

---

## How it works

1. **Sequence search** — submits your sequence to the RCSB PDB Search API (MMseqs2 backend, ≥30% identity)
2. **Metadata fetch** — queries the RCSB GraphQL API for crystallization records (`exptl_crystal_grow`), classifies X-ray vs. non-X-ray vs. AlphaFold/CSM entries
3. **Condition normalization** — parses free-text `pdbx_details` fields to extract precipitants, buffers, concentrations (96% extraction rate on 200-entry benchmark)
4. **Consensus** — aggregates pH, temperature, precipitant/buffer frequencies across contributing X-ray structures with a confidence tier (HIGH / MEDIUM / LOW)
5. **Biophysics panel** — pI, MW, instability index, GRAVY, and flags for high Cys/Pro/Gly content

All data comes from public RCSB APIs — no API keys required.

---

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

Open http://localhost:8501 in your browser.

---

## Project structure

| File | Purpose |
|------|---------|
| `app.py` | Streamlit UI |
| `sequence_search.py` | RCSB MMseqs2 search with async polling |
| `metadata_fetcher.py` | GraphQL metadata fetch + AlphaFold filter |
| `condition_normalizer.py` | Free-text condition parser (precipitants, buffers, units) |
| `consensus.py` | Condition aggregation and confidence scoring |
| `biophysics.py` | Biopython biophysical risk panel |
| `parser_spec.md` | Field format specification from 60-entry RCSB audit |
| `test_normalizer.py` | 200-entry benchmark for the condition normalizer |
| `step0_audit.py` | Development audit script (GraphQL field survey) |
