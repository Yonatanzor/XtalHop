"""
test_normalizer.py — Step 0b benchmark
Tests condition_normalizer against 200 live RCSB PDB entries.
Reports coverage, extraction rates, and flags cases needing attention.
"""

import json
import time
import random
import requests
from collections import Counter
from condition_normalizer import normalize, summarize

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"

GRAPHQL_QUERY = """
query GetCrystalGrow($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    exptl_crystal_grow {
      pH
      pdbx_pH_range
      temp
      method
      details
      pdbx_details
    }
    rcsb_entry_info {
      experimental_method
      resolution_combined
    }
  }
}
"""


def fetch_ids(n=220) -> list[str]:
    payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl_crystal_grow.pH",
                        "operator": "range",
                        "value": {"from": 0, "to": 14,
                                  "include_lower": True, "include_upper": True},
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 100, "rows": n},  # skip first 100 (used in audit)
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
    }
    r = requests.post(SEARCH_URL, json=payload, timeout=30)
    r.raise_for_status()
    ids = [h["identifier"] for h in r.json().get("result_set", [])]
    return ids


def fetch_graphql(ids: list[str]) -> list[dict]:
    r = requests.post(
        GRAPHQL_URL,
        json={"query": GRAPHQL_QUERY, "variables": {"ids": ids}},
        timeout=30,
    )
    r.raise_for_status()
    return r.json().get("data", {}).get("entries") or []


def main():
    print("=" * 70)
    print("CONDITION NORMALIZER — 200-ENTRY BENCHMARK")
    print("=" * 70)

    print("\n[1] Fetching 200 X-ray entry IDs (offset 100, skipping audit set)...")
    ids = fetch_ids(220)
    random.seed(42)
    random.shuffle(ids)
    ids = ids[:200]
    print(f"    Using {len(ids)} entries")

    print("\n[2] Fetching GraphQL data in batches of 50...")
    all_entries = []
    for i in range(0, len(ids), 50):
        batch = ids[i:i+50]
        entries = fetch_graphql(batch)
        all_entries.extend(entries)
        time.sleep(0.5)
        print(f"    Batch {i//50 + 1}: +{len(entries)} entries")

    print(f"\n    Total entries returned: {len(all_entries)}")

    # --- Metrics ---
    total = len(all_entries)
    n_with_pdbx = 0
    n_got_precipitant = 0
    n_got_buffer = 0
    n_got_either = 0
    n_warned = 0
    n_no_extract = 0
    precipitant_names: Counter = Counter()
    buffer_names: Counter = Counter()
    warning_types: Counter = Counter()
    failures = []       # entries where nothing extracted + has pdbx_details
    samples = []        # random sample for manual inspection

    for entry in all_entries:
        eid = entry["rcsb_id"]
        grow_records = entry.get("exptl_crystal_grow") or []

        for rec in grow_records:
            raw = rec.get("pdbx_details")
            if raw:
                n_with_pdbx += 1
            cond = normalize(raw)

            has_p = bool(cond.precipitants)
            has_b = bool(cond.buffers)

            if has_p:
                n_got_precipitant += 1
            if has_b:
                n_got_buffer += 1
            if has_p or has_b:
                n_got_either += 1
            if cond.warnings:
                n_warned += 1
                for w in cond.warnings:
                    warning_types[w] += 1
            if raw and not has_p and not has_b:
                n_no_extract += 1
                if len(failures) < 30:
                    failures.append((eid, raw, cond))

            for p in cond.precipitants:
                precipitant_names[p.name] += 1
            for b in cond.buffers:
                buffer_names[b.name] += 1

            if raw and len(samples) < 25:
                samples.append((eid, raw, cond))

    # --- Report ---
    print("\n" + "=" * 70)
    print("BENCHMARK RESULTS")
    print("=" * 70)

    print(f"\nEntries processed:           {total}")
    print(f"Records with pdbx_details:   {n_with_pdbx}")
    print(f"\nExtraction rates (of records with pdbx_details):")
    pct = lambda n: f"{n/n_with_pdbx*100:.1f}%" if n_with_pdbx else "N/A"
    print(f"  Got >=1 precipitant:        {n_got_precipitant:3d}  ({pct(n_got_precipitant)})")
    print(f"  Got >=1 buffer:             {n_got_buffer:3d}  ({pct(n_got_buffer)})")
    print(f"  Got >=1 of either:          {n_got_either:3d}  ({pct(n_got_either)})")
    print(f"  Extracted nothing:         {n_no_extract:3d}  ({pct(n_no_extract)})")
    print(f"  Any warning:               {n_warned:3d}  ({pct(n_warned)})")

    print(f"\nTop 20 precipitants found:")
    for name, count in precipitant_names.most_common(20):
        bar = "#" * min(30, count)
        print(f"  {name:<35} {count:3d}  {bar}")

    print(f"\nTop 15 buffers found:")
    for name, count in buffer_names.most_common(15):
        bar = "#" * min(30, count)
        print(f"  {name:<35} {count:3d}  {bar}")

    print(f"\nWarning types:")
    for w, count in warning_types.most_common():
        print(f"  {w:<40} {count:3d}")

    # --- Failures (has pdbx_details but nothing extracted) ---
    print(f"\n{'='*70}")
    print(f"FAILURE CASES (has pdbx_details, extracted nothing) — first 20")
    print(f"{'='*70}")
    for eid, raw, cond in failures[:20]:
        print(f"\n  [{eid}] {raw!r}")
        if cond.warnings:
            print(f"        warnings: {cond.warnings}")

    # --- Random sample for manual inspection ---
    print(f"\n{'='*70}")
    print(f"RANDOM SAMPLE — 15 entries for manual review")
    print(f"{'='*70}")
    for eid, raw, cond in samples[:15]:
        print(f"\n[{eid}]")
        print(f"  Raw:     {raw!r}")
        print(f"  Cleaned: {cond.raw_cleaned!r}")
        print(f"  {summarize(cond)}")

    # --- Benchmark verdict ---
    print(f"\n{'='*70}")
    print("VERDICT")
    print(f"{'='*70}")
    extraction_rate = n_got_either / n_with_pdbx * 100 if n_with_pdbx else 0
    precipitant_rate = n_got_precipitant / n_with_pdbx * 100 if n_with_pdbx else 0

    print(f"\n  Overall extraction rate: {extraction_rate:.1f}%")
    print(f"  Precipitant identification rate: {precipitant_rate:.1f}%")

    if precipitant_rate >= 85:
        print(f"\n  PASS — precipitant rate {precipitant_rate:.1f}% >= 85% target")
        print("  Proceed to building sequence_search.py and metadata_fetcher.py.")
    elif precipitant_rate >= 70:
        print(f"\n  MARGINAL — precipitant rate {precipitant_rate:.1f}% (target: 85%)")
        print("  Review failure cases above before proceeding.")
    else:
        print(f"\n  FAIL — precipitant rate {precipitant_rate:.1f}% < 70%")
        print("  Do not proceed. Fix the normalizer first.")


if __name__ == "__main__":
    main()
