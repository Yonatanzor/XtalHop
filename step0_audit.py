"""
Step 0: Manual RCSB GraphQL audit of exptl_crystal_grow fields.
Fetches 30 real X-ray PDB entries and dumps raw crystallization condition data
so we can document every format variation before writing the parser.
"""

import json
import requests
import time

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"


def get_xray_entry_ids(n=60):
    """Fetch PDB IDs for X-ray entries via RCSB Search API."""
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
                        "negation": False,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl_crystal_grow.pH",
                        "operator": "range",
                        "value": {
                            "from": 0, "to": 14,
                            "include_lower": True, "include_upper": True
                        },
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": n},
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
    }
    r = requests.post(SEARCH_URL, json=payload, timeout=30)
    r.raise_for_status()
    data = r.json()
    return [hit["identifier"] for hit in data.get("result_set", [])]


def fetch_crystal_grow_data(entry_ids):
    """Fetch exptl_crystal_grow fields via RCSB GraphQL for a list of entry IDs."""
    query = """
    query GetCrystalGrow($ids: [String!]!) {
      entries(entry_ids: $ids) {
        rcsb_id
        exptl {
          method
        }
        exptl_crystal_grow {
          pH
          pdbx_pH_range
          temp
          temp_details
          method
          details
          pdbx_details
        }
        rcsb_entry_info {
          resolution_combined
          deposited_atom_count
        }
      }
    }
    """
    r = requests.post(
        GRAPHQL_URL,
        json={"query": query, "variables": {"ids": entry_ids}},
        timeout=30,
    )
    r.raise_for_status()
    return r.json()


def main():
    print("=" * 70)
    print("RCSB exptl_crystal_grow RAW FIELD AUDIT")
    print("=" * 70)

    print("\n[1] Fetching X-ray entry IDs from RCSB Search API...")
    all_ids = get_xray_entry_ids(60)
    print(f"    Got {len(all_ids)} IDs: {all_ids[:10]}...")

    # Fetch in two batches of 30 (GraphQL has limits)
    batch1 = all_ids[:30]
    batch2 = all_ids[30:60]

    print("\n[2] Fetching crystallization data via GraphQL (batch 1)...")
    data1 = fetch_crystal_grow_data(batch1)
    time.sleep(1)
    print("[3] Fetching crystallization data via GraphQL (batch 2)...")
    data2 = fetch_crystal_grow_data(batch2)

    all_entries = (data1.get("data", {}).get("entries") or []) + (
        data2.get("data", {}).get("entries") or []
    )

    print(f"\n    Total entries returned: {len(all_entries)}")

    # Filter to those with at least one exptl_crystal_grow record
    with_data = [
        e for e in all_entries if e.get("exptl_crystal_grow")
    ]
    without_data = [
        e for e in all_entries if not e.get("exptl_crystal_grow")
    ]

    print(f"    With exptl_crystal_grow data:    {len(with_data)}")
    print(f"    Without exptl_crystal_grow data: {len(without_data)}")
    print(f"    Missing data %: {len(without_data)/len(all_entries)*100:.1f}%")

    # Print IDs missing data
    print(f"\n    Entries with NO crystallization data:")
    for e in without_data:
        method = (e.get("exptl") or [{}])[0].get("method", "?")
        print(f"      {e['rcsb_id']}  method={method}")

    print("\n" + "=" * 70)
    print("RAW FIELD DUMP — with_data entries")
    print("=" * 70)

    for entry in with_data:
        eid = entry["rcsb_id"]
        method = (entry.get("exptl") or [{}])[0].get("method", "?")
        res = (entry.get("rcsb_entry_info") or {}).get("resolution_combined")
        grow_records = entry.get("exptl_crystal_grow") or []

        print(f"\n{'-'*60}")
        print(f"PDB: {eid}  |  Method: {method}  |  Resolution: {res} Å")
        print(f"  exptl_crystal_grow records: {len(grow_records)}")

        for i, rec in enumerate(grow_records, 1):
            print(f"\n  [Record {i}]")
            print(f"    pH:           {rec.get('pH')!r}")
            print(f"    pdbx_ph_range:{rec.get('pdbx_ph_range')!r}")
            print(f"    temp:         {rec.get('temp')!r}")
            print(f"    temp_details: {rec.get('temp_details')!r}")
            print(f"    method:       {rec.get('method')!r}")
            print(f"    details:      {rec.get('details')!r}")
            print(f"    pdbx_details: {rec.get('pdbx_details')!r}")

    # Summary statistics
    print("\n" + "=" * 70)
    print("FIELD COVERAGE SUMMARY")
    print("=" * 70)

    field_counts = {
        "pH": 0, "pdbx_ph_range": 0, "temp": 0, "temp_details": 0,
        "method": 0, "details": 0, "pdbx_details": 0,
    }
    total_records = 0

    for entry in with_data:
        for rec in (entry.get("exptl_crystal_grow") or []):
            total_records += 1
            for field in field_counts:
                if rec.get(field) is not None:
                    field_counts[field] += 1

    print(f"\n  Total exptl_crystal_grow records across {len(with_data)} entries: {total_records}")
    print(f"\n  Field population rates:")
    for field, count in field_counts.items():
        pct = count / total_records * 100 if total_records else 0
        bar = "#" * int(pct / 5)
        print(f"    {field:<20} {count:3d}/{total_records}  ({pct:5.1f}%)  {bar}")

    # Raw pH values
    print("\n  All raw pH values found:")
    ph_vals = []
    for entry in with_data:
        for rec in (entry.get("exptl_crystal_grow") or []):
            if rec.get("pH") is not None:
                ph_vals.append((entry["rcsb_id"], rec["pH"]))
    for eid, ph in ph_vals[:40]:
        print(f"    {eid}: {ph!r}")

    # Raw details samples
    print("\n  Sample 'details' field values (first 20 non-null):")
    details_shown = 0
    for entry in with_data:
        for rec in (entry.get("exptl_crystal_grow") or []):
            if rec.get("details") and details_shown < 20:
                print(f"\n    [{entry['rcsb_id']}] details=")
                print(f"      {rec['details']!r}")
                details_shown += 1

    print("\n  Sample 'pdbx_details' field values (first 20 non-null):")
    pdbx_shown = 0
    for entry in with_data:
        for rec in (entry.get("exptl_crystal_grow") or []):
            if rec.get("pdbx_details") and pdbx_shown < 20:
                print(f"\n    [{entry['rcsb_id']}] pdbx_details=")
                print(f"      {rec['pdbx_details']!r}")
                pdbx_shown += 1

    print("\n\n[DONE] Save this output as parser_spec_raw.txt for reference.")


if __name__ == "__main__":
    main()
