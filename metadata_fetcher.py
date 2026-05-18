"""
metadata_fetcher.py — Module 2
Fetches structural metadata and crystallization conditions for a list of PDB entry IDs.
CRITICAL: AlphaFold (CSM) entries are explicitly filtered BEFORE the condition pipeline
          and flagged with a user-visible warning. Including them silently corrupts output.
"""

import time
import requests
from dataclasses import dataclass, field

GRAPHQL_URL = "https://data.rcsb.org/graphql"
BATCH_SIZE = 50   # RCSB GraphQL handles up to ~100 per call; 50 is safe

EXPERIMENTAL_METHOD_QUERY = """
query GetMetadata($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    rcsb_entry_info {
      experimental_method
      resolution_combined
    }
    exptl {
      method
    }
    exptl_crystal_grow {
      pH
      pdbx_pH_range
      temp
      temp_details
      method
      pdbx_details
    }
    polymer_entities {
      rcsb_id
      uniprots {
        rcsb_id
      }
    }
    symmetry {
      space_group_name_H_M
    }
    cell {
      length_a
      length_b
      length_c
      angle_alpha
      angle_beta
      angle_gamma
    }
    nonpolymer_entities {
      nonpolymer_comp {
        chem_comp {
          id
          name
        }
      }
    }
  }
}
"""

# Methods that produce structural models WITHOUT experimental crystallization data.
# These must NEVER be passed to the condition aggregation pipeline.
NON_CRYSTALLOGRAPHIC_METHODS = {
    "ELECTRON MICROSCOPY",       # Cryo-EM
    "SOLUTION NMR",
    "SOLID-STATE NMR",
    "NEUTRON DIFFRACTION",
    "ELECTRON CRYSTALLOGRAPHY",
    "FIBER DIFFRACTION",
    "POWDER DIFFRACTION",
    "SOLUTION SCATTERING",
    "EPR",
    "THEORETICAL MODEL",
    "AB INITIO MODEL",
    "ALPHAFOLD",                 # Computed Structure Models (CSMs)
    "ESMFOLD",
}

# RCSB flags AlphaFold/CSM entries with this prefix in their IDs
CSM_ENTRY_PREFIXES = ("AF-", "MA-")


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class UnitCell:
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float

    def format_compact(self) -> str:
        dims = f"{self.a:.1f}×{self.b:.1f}×{self.c:.1f}"
        non_ortho = not (self.alpha == 90.0 and self.beta == 90.0 and self.gamma == 90.0)
        if non_ortho:
            dims += f" ({self.alpha:.0f}/{self.beta:.0f}/{self.gamma:.0f}°)"
        return dims + " Å"


@dataclass
class CrystallizationRecord:
    ph: float | None
    temp_k: float | None
    temp_c: float | None    # converted from Kelvin
    method: str | None
    pdbx_details: str | None
    ph_range: str | None


@dataclass
class EntryMetadata:
    entry_id: str
    experimental_method: str | None
    resolution_a: float | None
    uniprot_ids: list[str]
    crystal_grow: list[CrystallizationRecord]
    # Classification flags
    is_xray: bool = False
    is_alphafold: bool = False
    is_non_crystallographic: bool = False
    has_crystallization_data: bool = False
    # Human-readable warning for the UI
    exclusion_reason: str | None = None
    # Crystal geometry
    space_group: str | None = None
    unit_cell: UnitCell | None = None
    # Ligands (non-polymer entities, water excluded)
    ligand_ids: list[str] = field(default_factory=list)


@dataclass
class FetchResult:
    entries: list[EntryMetadata] = field(default_factory=list)
    xray_entries: list[EntryMetadata] = field(default_factory=list)
    non_xray_entries: list[EntryMetadata] = field(default_factory=list)
    alphafold_entries: list[EntryMetadata] = field(default_factory=list)
    missing_data_pct: float = 0.0
    elapsed_s: float = 0.0


# ---------------------------------------------------------------------------
# GraphQL fetch
# ---------------------------------------------------------------------------

def _fetch_batch(entry_ids: list[str]) -> list[dict]:
    resp = requests.post(
        GRAPHQL_URL,
        json={"query": EXPERIMENTAL_METHOD_QUERY, "variables": {"ids": entry_ids}},
        timeout=30,
    )
    resp.raise_for_status()
    return resp.json().get("data", {}).get("entries") or []


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def _extract_uniprot(entry: dict) -> list[str]:
    ids = []
    for entity in (entry.get("polymer_entities") or []):
        for uni in (entity.get("uniprots") or []):
            acc = uni.get("rcsb_id")
            if acc:
                ids.append(acc)
    return list(dict.fromkeys(ids))  # deduplicate, preserve order


def _parse_crystal_grow(records: list[dict]) -> list[CrystallizationRecord]:
    result = []
    for rec in (records or []):
        temp_k = rec.get("temp")
        temp_c = round(temp_k - 273.15, 2) if temp_k is not None else None
        result.append(CrystallizationRecord(
            ph=rec.get("pH"),
            temp_k=temp_k,
            temp_c=temp_c,
            method=rec.get("method"),
            pdbx_details=rec.get("pdbx_details"),
            ph_range=rec.get("pdbx_pH_range"),
        ))
    return result


def _extract_unit_cell(raw: dict) -> UnitCell | None:
    cell = raw.get("cell")
    if not cell:
        return None
    try:
        return UnitCell(
            a=float(cell["length_a"]),
            b=float(cell["length_b"]),
            c=float(cell["length_c"]),
            alpha=float(cell["angle_alpha"]),
            beta=float(cell["angle_beta"]),
            gamma=float(cell["angle_gamma"]),
        )
    except (TypeError, KeyError, ValueError):
        return None


_SOLVENT_IDS = {"HOH", "DOD", "H2O", "WAT"}


def _extract_ligands(raw: dict) -> list[str]:
    """Return deduplicated ligand IDs (3-letter codes), excluding water."""
    ids = []
    for entity in (raw.get("nonpolymer_entities") or []):
        comp = (entity.get("nonpolymer_comp") or {})
        chem = (comp.get("chem_comp") or {})
        chem_id = chem.get("id")
        if chem_id and chem_id.upper() not in _SOLVENT_IDS:
            ids.append(chem_id)
    return list(dict.fromkeys(ids))


def _classify(entry_id: str, raw: dict) -> EntryMetadata:
    """Parse a raw GraphQL entry dict into a classified EntryMetadata."""
    info = raw.get("rcsb_entry_info") or {}
    exptl_list = raw.get("exptl") or [{}]
    exptl_method = (exptl_list[0].get("method") or "").upper().strip()
    method_combined = info.get("experimental_method", exptl_method) or exptl_method

    resolution = info.get("resolution_combined")
    if isinstance(resolution, list):
        resolution = resolution[0] if resolution else None

    uniprot_ids = _extract_uniprot(raw)
    crystal_records = _parse_crystal_grow(raw.get("exptl_crystal_grow"))

    # rcsb_entry_info.experimental_method uses short codes: "X-ray", "EM", "NMR", etc.
    # exptl[].method uses long form: "X-RAY DIFFRACTION", "ELECTRON MICROSCOPY", etc.
    method_upper = method_combined.upper()

    is_alphafold = (
        entry_id.startswith(CSM_ENTRY_PREFIXES)
        or "ALPHAFOLD" in method_upper
        or "THEORETICAL" in method_upper
        or "AB INITIO" in method_upper
    )

    is_xray = (
        not is_alphafold
        and ("X-RAY" in method_upper or method_upper == "X-RAY")
    )

    is_non_cryst = (
        not is_xray
        and not is_alphafold
        and method_upper not in ("", "UNKNOWN")
    )

    has_data = bool(crystal_records and any(r.pdbx_details or r.ph for r in crystal_records))

    # Build exclusion reason for UI display
    exclusion_reason = None
    if is_alphafold:
        exclusion_reason = "Computed structure model (AlphaFold/CSM) — no experimental crystallization data"
    elif is_non_cryst:
        exclusion_reason = f"Non-crystallographic method ({method_combined}) — no crystallization conditions available"

    return EntryMetadata(
        entry_id=entry_id,
        experimental_method=method_combined or None,
        resolution_a=resolution,
        uniprot_ids=uniprot_ids,
        crystal_grow=crystal_records,
        is_xray=is_xray,
        is_alphafold=is_alphafold,
        is_non_crystallographic=is_non_cryst,
        has_crystallization_data=has_data,
        exclusion_reason=exclusion_reason,
        space_group=(raw.get("symmetry") or {}).get("space_group_name_H_M"),
        unit_cell=_extract_unit_cell(raw),
        ligand_ids=_extract_ligands(raw),
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def fetch(entry_ids: list[str], verbose: bool = False) -> FetchResult:
    """
    Fetch and classify metadata for a list of PDB entry IDs.

    AlphaFold and non-crystallographic entries are tagged but NOT passed
    to the crystallization data pipeline — they appear in the UI table
    with an explicit flag.

    Args:
        entry_ids: List of PDB entry IDs (e.g. ["1LYZ", "6VXX"]).
        verbose: Print progress to stdout.

    Returns:
        FetchResult with classified entry lists and metadata statistics.
    """
    if not entry_ids:
        return FetchResult()

    t0 = time.time()
    all_entries: list[EntryMetadata] = []

    if verbose:
        print(f"[fetch] Fetching metadata for {len(entry_ids)} entries in batches of {BATCH_SIZE}...")

    for i in range(0, len(entry_ids), BATCH_SIZE):
        batch = entry_ids[i:i + BATCH_SIZE]
        raw_entries = _fetch_batch(batch)
        for raw in raw_entries:
            eid = raw.get("rcsb_id", "")
            entry = _classify(eid, raw)
            all_entries.append(entry)
        if verbose:
            print(f"[fetch]   batch {i//BATCH_SIZE + 1}: +{len(raw_entries)}")
        if i + BATCH_SIZE < len(entry_ids):
            time.sleep(0.3)

    xray = [e for e in all_entries if e.is_xray]
    alphafold = [e for e in all_entries if e.is_alphafold]
    non_xray = [e for e in all_entries if e.is_non_crystallographic]

    # Missing data % counts X-ray entries that lack crystallization records
    xray_missing = [e for e in xray if not e.has_crystallization_data]
    missing_pct = (
        len(xray_missing) / len(xray) * 100 if xray else 0.0
    )

    elapsed = round(time.time() - t0, 2)
    if verbose:
        print(f"[fetch] Done in {elapsed}s")
        print(f"  X-ray:       {len(xray)}")
        print(f"  AlphaFold:   {len(alphafold)}  (excluded from condition pipeline)")
        print(f"  Non-X-ray:   {len(non_xray)}  (excluded from condition pipeline)")
        print(f"  Missing data: {missing_pct:.1f}% of X-ray entries lack crystallization records")

    return FetchResult(
        entries=all_entries,
        xray_entries=xray,
        non_xray_entries=non_xray,
        alphafold_entries=alphafold,
        missing_data_pct=missing_pct,
        elapsed_s=elapsed,
    )


def branching_decision(result: FetchResult) -> str:
    """
    Apply the Phase 2 branching thresholds from XtalHop spec.
    Returns: 'condition_a', 'condition_b', or 'buffer_zone'
    """
    pct = result.missing_data_pct
    if pct < 10:
        return "condition_a"    # < 10% missing: proceed normally, show N/A
    elif pct >= 20:
        return "condition_b"    # >= 20% missing: activate Cluster-Hopping Engine
    else:
        return "buffer_zone"    # 10–19%: Condition A + ensure >= 5 valid X-ray entries


# ---------------------------------------------------------------------------
# CLI test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Test with a mix: X-ray lysozyme entries, a Cryo-EM entry, and an AlphaFold entry
    TEST_IDS = [
        "1LYZ",   # X-ray lysozyme (old, likely no pdbx_details)
        "193L",   # X-ray lysozyme with crystallization data
        "1AKI",   # X-ray lysozyme
        "6VXX",   # Cryo-EM SARS-CoV-2 spike — must be flagged non-xray
        "7N3C",   # Cryo-EM — non-crystallographic
        "4HHB",   # X-ray hemoglobin
        "10GW",   # X-ray with pdbx_details confirmed in audit
        "10CH",   # X-ray with pdbx_details confirmed in audit
    ]

    print("=" * 65)
    print("metadata_fetcher.py — live test")
    print("=" * 65)

    result = fetch(TEST_IDS, verbose=True)

    print("\n" + "=" * 65)
    print("ENTRY CLASSIFICATION TABLE")
    print("=" * 65)
    print(f"{'Entry':<12} {'Method':<30} {'Res':<6} {'UniProt':<12} {'XtalData':<9} {'Flag'}")
    print("-" * 90)
    for e in result.entries:
        method_short = (e.experimental_method or "?")[:28]
        res = f"{e.resolution_a:.2f}" if e.resolution_a else "N/A"
        uniprot = ",".join(e.uniprot_ids[:2]) or "—"
        has_data = "YES" if e.has_crystallization_data else "NO"
        flag = ""
        if e.is_alphafold:
            flag = "[CSM — excluded]"
        elif e.is_non_crystallographic:
            flag = "[Non-Xray — excluded]"
        elif not e.has_crystallization_data:
            flag = "[missing data]"
        print(f"{e.entry_id:<12} {method_short:<30} {res:<6} {uniprot:<12} {has_data:<9} {flag}")

    print(f"\nMissing data (X-ray only): {result.missing_data_pct:.1f}%")
    branching = branching_decision(result)
    print(f"Branching decision: {branching.upper()}")

    print("\nSample crystallization record (1LYZ):")
    lyz = next((e for e in result.entries if e.entry_id == "1LYZ"), None)
    if lyz and lyz.crystal_grow:
        r = lyz.crystal_grow[0]
        print(f"  pH={r.ph}  temp={r.temp_c}C  method={r.method}")
        print(f"  pdbx_details={r.pdbx_details!r}")
    elif lyz:
        print("  No crystallization records (1LYZ pre-dates structured deposition)")
