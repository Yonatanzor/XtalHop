"""
sequence_search.py — Module 1
Submits a protein sequence to the RCSB PDB Search API (MMseqs2 backend),
polls for completion, and returns ranked homolog hits.

RCSB sequence search is asynchronous: the API returns a job URL on submission
and requires polling until the job is ready. This module handles the full
submit → poll → retrieve cycle transparently.
"""

import time
import requests
from dataclasses import dataclass, field

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"

# Council-mandated constants (never change without re-benchmarking)
IDENTITY_CUTOFF = 0.30   # >= 30% sequence identity
EVALUE_CUTOFF = 0.1
POLL_INTERVAL_S = 3      # seconds between polls
MAX_POLL_ATTEMPTS = 40   # 40 × 3 s = 2-minute max wait


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class SequenceHit:
    entry_id: str           # PDB entry ID, e.g. "1LYZ"
    entity_id: str          # polymer entity ID, e.g. "1LYZ_1"
    identity: float         # sequence identity 0.0–1.0
    evalue: float | None
    bitscore: float | None
    rank: int


@dataclass
class SearchResult:
    query_sequence: str
    hits: list[SequenceHit] = field(default_factory=list)
    total_count: int = 0
    elapsed_s: float = 0.0
    error: str | None = None


# ---------------------------------------------------------------------------
# Payload builder
# ---------------------------------------------------------------------------

def _build_payload(sequence: str, n: int) -> dict:
    """Build the RCSB v2 search API payload for MMseqs2 sequence search."""
    return {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": EVALUE_CUTOFF,
                "identity_cutoff": IDENTITY_CUTOFF,
                "sequence_type": "protein",
                "value": sequence.strip(),
            },
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": min(n, 50)},
            "sort": [
                {"sort_by": "score", "direction": "desc"}
            ],
            "scoring_strategy": "sequence",
        },
    }


# ---------------------------------------------------------------------------
# Submission + polling
# ---------------------------------------------------------------------------

def _submit(payload: dict, timeout: int = 20) -> requests.Response:
    """POST the search payload and return the raw response."""
    return requests.post(SEARCH_URL, json=payload, timeout=timeout)


def _parse_hits(data: dict) -> list[SequenceHit]:
    hits = []
    for rank, item in enumerate(data.get("result_set", []), start=1):
        identifier = item.get("identifier", "")
        # identifier format: "ENTRYID_ENTITYIDX", e.g. "1LYZ_1"
        parts = identifier.split("_")
        entry_id = parts[0] if parts else identifier
        services = item.get("services", [])
        identity = None
        evalue = None
        bitscore = None
        for svc in services:
            for node in svc.get("nodes", []):
                orig = node.get("original_score")
                if orig is not None:
                    identity = orig  # sequence identity 0–1
                stats = node.get("match_context", [{}])[0] if node.get("match_context") else {}
                evalue = stats.get("evalue", evalue)
                bitscore = stats.get("bitscore", bitscore)
        # Fallback: score field at top level
        if identity is None:
            identity = item.get("score", 0.0)
        hits.append(SequenceHit(
            entry_id=entry_id,
            entity_id=identifier,
            identity=round(float(identity), 4),
            evalue=evalue,
            bitscore=bitscore,
            rank=rank,
        ))
    return hits


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def search(
    sequence: str,
    n: int = 20,
    verbose: bool = False,
) -> SearchResult:
    """
    Search the RCSB PDB for protein homologs using MMseqs2.

    Args:
        sequence: Raw amino acid sequence (FASTA or plain text).
        n: Number of top hits to return (capped at 50 for v1).
        verbose: Print polling status to stdout.

    Returns:
        SearchResult with ranked SequenceHit list.
    """
    # Strip FASTA header if present
    lines = sequence.strip().splitlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))

    if not seq:
        return SearchResult(query_sequence=sequence, error="empty sequence after parsing")

    payload = _build_payload(seq, n)
    t0 = time.time()

    if verbose:
        print(f"[search] Submitting sequence ({len(seq)} aa) to RCSB MMseqs2...")

    # --- Submit ---
    try:
        resp = _submit(payload, timeout=30)
    except requests.exceptions.Timeout:
        return SearchResult(query_sequence=seq, error="submission timed out")
    except requests.exceptions.RequestException as e:
        return SearchResult(query_sequence=seq, error=f"submission failed: {e}")

    # --- Handle async (202 Accepted with Location header) ---
    if resp.status_code == 202:
        poll_url = resp.headers.get("Location")
        if not poll_url:
            return SearchResult(query_sequence=seq,
                                error="202 Accepted but no Location header for polling")
        if verbose:
            print(f"[search] Async job accepted. Polling: {poll_url}")

        for attempt in range(1, MAX_POLL_ATTEMPTS + 1):
            time.sleep(POLL_INTERVAL_S)
            if verbose:
                print(f"[search]   poll attempt {attempt}/{MAX_POLL_ATTEMPTS}...")
            try:
                poll_resp = requests.get(poll_url, timeout=20)
            except requests.exceptions.RequestException as e:
                return SearchResult(query_sequence=seq, error=f"poll failed: {e}")

            if poll_resp.status_code == 200:
                resp = poll_resp
                break
            elif poll_resp.status_code in (202, 204):
                continue  # still processing
            else:
                return SearchResult(
                    query_sequence=seq,
                    error=f"poll returned HTTP {poll_resp.status_code}: {poll_resp.text[:200]}",
                )
        else:
            return SearchResult(
                query_sequence=seq,
                error=f"timed out after {MAX_POLL_ATTEMPTS * POLL_INTERVAL_S}s polling",
            )

    # --- Handle synchronous success ---
    if resp.status_code == 200:
        try:
            data = resp.json()
        except Exception as e:
            return SearchResult(query_sequence=seq, error=f"JSON parse error: {e}")

        hits = _parse_hits(data)
        elapsed = round(time.time() - t0, 2)

        if verbose:
            print(f"[search] Done in {elapsed}s. Total hits: {data.get('total_count', '?')}")

        return SearchResult(
            query_sequence=seq,
            hits=hits,
            total_count=data.get("total_count", len(hits)),
            elapsed_s=elapsed,
        )

    # --- Handle 204 No Content (no hits) ---
    if resp.status_code == 204:
        return SearchResult(
            query_sequence=seq,
            hits=[],
            total_count=0,
            elapsed_s=round(time.time() - t0, 2),
        )

    return SearchResult(
        query_sequence=seq,
        error=f"HTTP {resp.status_code}: {resp.text[:300]}",
    )


# ---------------------------------------------------------------------------
# CLI test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Hen egg-white lysozyme — canonical test case, well-crystallized,
    # thousands of PDB homologs expected
    LYSOZYME = (
        "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRT"
        "PGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
    )

    print("=" * 60)
    print("sequence_search.py — live test with hen lysozyme")
    print("=" * 60)

    result = search(LYSOZYME, n=20, verbose=True)

    if result.error:
        print(f"\nERROR: {result.error}")
    else:
        print(f"\nTotal PDB homologs found: {result.total_count}")
        print(f"Elapsed: {result.elapsed_s}s")
        print(f"\nTop {len(result.hits)} hits (identity >= {IDENTITY_CUTOFF*100:.0f}%):\n")
        print(f"{'Rank':<5} {'Entry':<8} {'Entity ID':<12} {'Identity':>10}")
        print("-" * 40)
        for h in result.hits:
            print(f"{h.rank:<5} {h.entry_id:<8} {h.entity_id:<12} {h.identity*100:>9.1f}%")
