"""
app.py — XtalHop-Predictor
Streamlit UI: protein sequence -> PDB homologs -> consensus crystallization conditions.
"""

import io
import time
import pandas as pd
import streamlit as st

from sequence_search import search as seq_search, SearchResult, SequenceHit
from metadata_fetcher import fetch as meta_fetch, FetchResult, EntryMetadata
from condition_normalizer import normalize, summarize
from consensus import build_consensus, ConsensusResult
from biophysics import analyze as bio_analyze, BiophysicsResult

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="XtalHop-Predictor",
    page_icon="🔬",
    layout="wide",
)

st.title("🔬 XtalHop-Predictor")
st.caption(
    "Protein sequence → PDB structural homologs → consensus crystallization conditions"
)

# ---------------------------------------------------------------------------
# Sidebar — search parameters
# ---------------------------------------------------------------------------

with st.sidebar:
    st.header("Search Parameters")
    n_hits = st.slider("Max homologs to retrieve", 10, 50, 20, step=5)
    identity_pct = st.number_input(
        "Min sequence identity (%)", min_value=10, max_value=100, value=30, step=5
    )
    st.caption(
        "Identity threshold is fixed at 30% in the search API for this release. "
        "The slider filters results after retrieval."
    )
    st.divider()
    st.header("About")
    st.markdown(
        "XtalHop queries the RCSB PDB for structural homologs of your sequence "
        "and aggregates their crystallization conditions into a consensus formulation.\n\n"
        "**Provenance is always shown** — you can see exactly how many structures "
        "contributed and the data coverage rate."
    )

# ---------------------------------------------------------------------------
# Sequence input
# ---------------------------------------------------------------------------

EXAMPLE_LYSOZYME = (
    "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRT"
    "PGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
)

if "sequence" not in st.session_state:
    st.session_state["sequence"] = ""

col_run, col_example, _ = st.columns([1, 1, 4])
with col_example:
    if st.button("Load lysozyme example", use_container_width=True):
        st.session_state["sequence"] = EXAMPLE_LYSOZYME
        st.rerun()

sequence_input = st.text_area(
    "Protein sequence (FASTA or plain amino acids)",
    placeholder="Paste your sequence here...",
    height=150,
    key="sequence",
)

with col_run:
    run_clicked = st.button("Search PDB", type="primary", use_container_width=True)

# ---------------------------------------------------------------------------
# Run pipeline
# ---------------------------------------------------------------------------

if run_clicked:
    if not sequence_input.strip():
        st.error("Please enter a protein sequence.")
        st.stop()

    # Strip FASTA header
    lines = sequence_input.strip().splitlines()
    seq = "".join(l.strip() for l in lines if not l.startswith(">"))

    # ---- Biophysics panel ------------------------------------------------
    with st.expander("Biophysical Risk Panel", expanded=True):
        try:
            bio: BiophysicsResult = bio_analyze(seq)
            bcol1, bcol2, bcol3, bcol4 = st.columns(4)
            bcol1.metric("pI", f"{bio.pi:.2f}")
            bcol2.metric("MW", f"{bio.molecular_weight_da/1000:.1f} kDa")
            bcol3.metric("Instability Index", f"{bio.instability_index:.1f}",
                         delta="stable" if bio.instability_index <= 40 else "unstable",
                         delta_color="normal" if bio.instability_index <= 40 else "inverse")
            bcol4.metric("GRAVY", f"{bio.gravy:+.3f}")

            if bio.risks:
                st.warning("**Crystallizability risks detected:**\n\n" +
                           "\n".join(f"- {r}" for r in bio.risks))
            else:
                st.success("No biophysical risk flags — good crystallization candidate.")
        except Exception as e:
            st.warning(f"Biophysics panel failed: {e}")

    # ---- Sequence search -------------------------------------------------
    st.subheader("Step 1 — Sequence search (MMseqs2)")
    search_placeholder = st.empty()

    with search_placeholder.container():
        with st.spinner("Searching RCSB PDB for structural homologs..."):
            t_search = time.time()
            sr: SearchResult = seq_search(seq, n=n_hits)
            elapsed_search = round(time.time() - t_search, 1)

    if sr.error:
        st.error(f"Sequence search failed: {sr.error}")
        st.stop()

    # Filter by user-selected identity cutoff (post-hoc; API floor is 30%)
    identity_floor = identity_pct / 100.0
    hits_filtered = [h for h in sr.hits if h.identity >= identity_floor]

    search_placeholder.success(
        f"Found **{sr.total_count:,}** total homologs in {elapsed_search}s. "
        f"Showing top {len(hits_filtered)} with identity ≥ {identity_pct}%."
    )

    if not hits_filtered:
        st.warning("No hits above the selected identity threshold.")
        st.stop()

    entry_ids = list(dict.fromkeys(h.entry_id for h in hits_filtered))

    # ---- Metadata fetch --------------------------------------------------
    st.subheader("Step 2 — Metadata & crystallization data")
    meta_placeholder = st.empty()

    with meta_placeholder.container():
        with st.spinner(f"Fetching metadata for {len(entry_ids)} entries..."):
            t_meta = time.time()
            fr: FetchResult = meta_fetch(entry_ids)
            elapsed_meta = round(time.time() - t_meta, 1)

    meta_placeholder.success(
        f"Metadata fetched in {elapsed_meta}s — "
        f"{len(fr.xray_entries)} X-ray, "
        f"{len(fr.non_xray_entries)} non-X-ray, "
        f"{len(fr.alphafold_entries)} AlphaFold/CSM. "
        f"Missing crystallization data: {fr.missing_data_pct:.1f}% of X-ray entries."
    )

    # Build identity lookup for the table
    identity_map: dict[str, float] = {}
    for h in hits_filtered:
        if h.entry_id not in identity_map:
            identity_map[h.entry_id] = h.identity

    # ---- Homolog table ---------------------------------------------------
    st.subheader("Homolog Table")

    def _flag(e: EntryMetadata) -> str:
        if e.is_alphafold:
            return "CSM — no xtal data"
        if e.is_non_crystallographic:
            return f"Non-X-ray ({e.experimental_method or '?'})"
        if not e.has_crystallization_data:
            return "X-ray — missing data"
        return ""

    def _condition_short(e: EntryMetadata) -> str:
        if not e.crystal_grow:
            return ""
        for rec in e.crystal_grow:
            if rec.pdbx_details:
                cond = normalize(rec.pdbx_details)
                s = summarize(cond)
                return s[:80] + "..." if len(s) > 80 else s
        return ""

    rows = []
    for e in fr.entries:
        ident = identity_map.get(e.entry_id, 0.0)
        condition_short = _condition_short(e) if e.is_xray else ""
        ph_val = ""
        if e.crystal_grow:
            phs = [r.ph for r in e.crystal_grow if r.ph is not None]
            if phs:
                ph_val = f"{phs[0]:.1f}"
        rows.append({
            "Entry": e.entry_id,
            "Method": e.experimental_method or "?",
            "Identity %": f"{ident*100:.1f}",
            "Resolution Å": f"{e.resolution_a:.2f}" if e.resolution_a else "N/A",
            "UniProt": ", ".join(e.uniprot_ids[:2]) or "—",
            "pH": ph_val or "—",
            "Condition (preview)": condition_short or (e.exclusion_reason or ""),
            "Flag": _flag(e),
        })

    df_table = pd.DataFrame(rows)

    def _row_style(row):
        flag = row["Flag"]
        if "CSM" in flag:
            return ["background-color: #fff3cd"] * len(row)  # amber
        if "Non-X-ray" in flag:
            return ["background-color: #e2e3e5"] * len(row)  # grey
        if "missing" in flag:
            return ["background-color: #f8d7da"] * len(row)  # light red
        return [""] * len(row)

    st.dataframe(
        df_table.style.apply(_row_style, axis=1),
        use_container_width=True,
        hide_index=True,
    )

    # CSV export — homolog table
    csv_table = df_table.to_csv(index=False).encode("utf-8")
    st.download_button(
        "Download homolog table (CSV)",
        csv_table,
        file_name="xtalhop_homologs.csv",
        mime="text/csv",
    )

    # ---- Consensus -------------------------------------------------------
    st.subheader("Consensus Crystallization Formulation")

    if not fr.xray_entries:
        st.warning("No X-ray entries found — cannot compute consensus.")
    else:
        with st.spinner("Aggregating conditions..."):
            cons: ConsensusResult = build_consensus(fr.xray_entries)

        # Provenance banner
        conf_colors = {
            "HIGH": "success",
            "MEDIUM": "info",
            "LOW": "warning",
            "INSUFFICIENT": "error",
        }
        conf_fn = getattr(st, conf_colors.get(cons.confidence, "info"))
        conf_fn(
            f"**Confidence: {cons.confidence}** — "
            f"{cons.n_contributing} of {cons.n_xray_entries} X-ray entries "
            f"contributed crystallization data ({cons.coverage_pct:.0f}% coverage)."
        )

        if cons.warnings:
            for w in cons.warnings:
                st.warning(w)

        # Core conditions
        ccol1, ccol2, ccol3, ccol4 = st.columns(4)
        ccol1.metric(
            "Median pH",
            f"{cons.median_ph:.1f}" if cons.median_ph else "N/A",
            help=f"Range: {cons.ph_range[0]}–{cons.ph_range[1]}" if cons.ph_range else None,
        )
        ccol2.metric(
            "Median Temp",
            f"{cons.median_temp_c:.0f} °C" if cons.median_temp_c else "N/A",
        )
        ccol3.metric(
            "Top Method",
            cons.top_method.title() if cons.top_method else "N/A",
        )
        ccol4.metric("X-ray Coverage", f"{cons.coverage_pct:.0f}%")

        # Precipitants
        if cons.precipitants:
            st.markdown("**Top Precipitants**")
            prec_rows = []
            for c in cons.precipitants[:10]:
                conc_str = ""
                if c.median_conc is not None and c.unit:
                    conc_str = f"{c.median_conc} {c.unit}"
                elif c.median_conc is not None:
                    conc_str = str(c.median_conc)
                prec_rows.append({
                    "Precipitant": c.name,
                    "Frequency": f"{c.frequency} entries ({c.frequency_pct:.0f}%)",
                    "Median Conc.": conc_str,
                })
            st.dataframe(pd.DataFrame(prec_rows), use_container_width=True, hide_index=True)

        # Buffers
        if cons.buffers:
            st.markdown("**Top Buffers**")
            buf_rows = []
            for c in cons.buffers[:8]:
                conc_str = ""
                if c.median_conc is not None and c.unit:
                    conc_str = f"{c.median_conc} {c.unit}"
                elif c.median_conc is not None:
                    conc_str = str(c.median_conc)
                buf_rows.append({
                    "Buffer": c.name,
                    "Frequency": f"{c.frequency} entries ({c.frequency_pct:.0f}%)",
                    "Median Conc.": conc_str,
                })
            st.dataframe(pd.DataFrame(buf_rows), use_container_width=True, hide_index=True)

        # Additives
        if cons.additives:
            with st.expander(f"Additives ({len(cons.additives)})"):
                add_rows = [
                    {"Additive": c.name, "Frequency": f"{c.frequency} entries"}
                    for c in cons.additives[:10]
                ]
                st.dataframe(pd.DataFrame(add_rows), use_container_width=True, hide_index=True)

        # Provenance links
        with st.expander("Contributing PDB entries"):
            links = " · ".join(
                f"[{eid}](https://www.rcsb.org/structure/{eid})"
                for eid in cons.contributing_entry_ids
            )
            st.markdown(links or "None")

        # CSV export — consensus
        if cons.precipitants or cons.buffers:
            cons_rows = []
            for c in cons.precipitants + cons.buffers + cons.additives:
                cons_rows.append({
                    "Chemical": c.name,
                    "Role": c.role,
                    "Frequency (entries)": c.frequency,
                    "Frequency (%)": c.frequency_pct,
                    "Median Concentration": c.median_conc,
                    "Unit": c.unit,
                })
            df_cons = pd.DataFrame(cons_rows)
            csv_cons = df_cons.to_csv(index=False).encode("utf-8")
            st.download_button(
                "Download consensus matrix (CSV)",
                csv_cons,
                file_name="xtalhop_consensus.csv",
                mime="text/csv",
            )

else:
    st.info("Enter a protein sequence and click **Search PDB** to begin.")
