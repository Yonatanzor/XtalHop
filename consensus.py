"""
consensus.py
Aggregates normalized crystallization conditions from multiple PDB entries
into a consensus starting formulation.
"""

import statistics
from collections import defaultdict
from dataclasses import dataclass, field
from condition_normalizer import normalize, NormalizedCondition
from metadata_fetcher import EntryMetadata


@dataclass
class ConsensusChemical:
    name: str
    role: str                   # 'precipitant' | 'buffer' | 'additive'
    frequency: int              # how many entries have it
    frequency_pct: float        # fraction of contributing entries
    concentrations: list[float] = field(default_factory=list)
    median_conc: float | None = None
    unit: str | None = None


@dataclass
class ConsensusResult:
    # Provenance
    n_xray_entries: int         # X-ray entries queried
    n_contributing: int         # X-ray entries with usable condition data
    coverage_pct: float         # n_contributing / n_xray_entries * 100

    # Confidence tier
    confidence: str             # 'HIGH' | 'MEDIUM' | 'LOW' | 'INSUFFICIENT'

    # Aggregated values
    median_ph: float | None
    ph_range: tuple[float, float] | None
    median_temp_c: float | None
    top_method: str | None

    # Top chemicals (sorted by frequency desc)
    precipitants: list[ConsensusChemical]
    buffers: list[ConsensusChemical]
    additives: list[ConsensusChemical]

    # Source entry IDs for provenance links
    contributing_entry_ids: list[str]

    warnings: list[str] = field(default_factory=list)


def _confidence_tier(n_contributing: int, coverage_pct: float) -> str:
    if n_contributing < 3:
        return "INSUFFICIENT"
    if n_contributing >= 10 and coverage_pct >= 80:
        return "HIGH"
    if n_contributing >= 5 and coverage_pct >= 50:
        return "MEDIUM"
    return "LOW"


def build_consensus(xray_entries: list[EntryMetadata]) -> ConsensusResult:
    """
    Compute consensus crystallization conditions from X-ray EntryMetadata list.
    Only entries with has_crystallization_data=True contribute to chemical stats.
    Structured pH and temp fields are used directly; pdbx_details drives
    chemical identification.
    """
    n_xray = len(xray_entries)
    contributing_ids: list[str] = []
    ph_values: list[float] = []
    temp_values: list[float] = []
    method_counts: dict[str, int] = defaultdict(int)

    # Per-chemical accumulators: {canonical_name: {unit: [conc, ...]}}
    chem_role: dict[str, str] = {}
    chem_units: dict[str, str] = {}
    chem_concs: dict[str, list[float]] = defaultdict(list)
    chem_entries: dict[str, set[str]] = defaultdict(set)

    for entry in xray_entries:
        if not entry.has_crystallization_data:
            continue
        contributing_ids.append(entry.entry_id)

        for rec in entry.crystal_grow:
            # Structured pH
            if rec.ph is not None:
                ph_values.append(rec.ph)

            # Structured temperature
            if rec.temp_c is not None:
                temp_values.append(rec.temp_c)

            # Method
            if rec.method:
                method_counts[rec.method.upper()] += 1

            # Chemical extraction from free text
            if rec.pdbx_details:
                cond: NormalizedCondition = normalize(rec.pdbx_details)
                for chem_list, role in (
                    (cond.precipitants, "precipitant"),
                    (cond.buffers, "buffer"),
                    (cond.additives, "additive"),
                ):
                    for chem in chem_list:
                        name = chem.name
                        if name not in chem_role:
                            chem_role[name] = role
                        chem_entries[name].add(entry.entry_id)
                        if chem.concentration is not None:
                            chem_concs[name].append(chem.concentration)
                        if chem.unit and name not in chem_units:
                            chem_units[name] = chem.unit

    n_contributing = len(contributing_ids)
    coverage_pct = (n_contributing / n_xray * 100) if n_xray else 0.0

    # pH stats
    median_ph = round(statistics.median(ph_values), 2) if ph_values else None
    ph_range = (round(min(ph_values), 1), round(max(ph_values), 1)) if ph_values else None

    # Temperature stats
    median_temp_c = round(statistics.median(temp_values), 1) if temp_values else None

    # Top method
    top_method = max(method_counts, key=method_counts.get) if method_counts else None

    # Build chemical objects, sorted by frequency
    chemicals: list[ConsensusChemical] = []
    for name, entry_set in chem_entries.items():
        concs = chem_concs[name]
        median_conc = round(statistics.median(concs), 3) if concs else None
        chemicals.append(ConsensusChemical(
            name=name,
            role=chem_role.get(name, "precipitant"),
            frequency=len(entry_set),
            frequency_pct=round(len(entry_set) / n_contributing * 100, 1) if n_contributing else 0,
            concentrations=concs,
            median_conc=median_conc,
            unit=chem_units.get(name),
        ))

    chemicals.sort(key=lambda c: -c.frequency)

    precipitants = [c for c in chemicals if c.role == "precipitant"]
    buffers = [c for c in chemicals if c.role == "buffer"]
    additives = [c for c in chemicals if c.role == "additive"]

    warnings: list[str] = []
    confidence = _confidence_tier(n_contributing, coverage_pct)
    if confidence == "INSUFFICIENT":
        warnings.append(
            f"Only {n_contributing} entries with crystallization data — "
            f"consensus is unreliable. Consider expanding the search."
        )
    elif confidence == "LOW":
        warnings.append(
            f"Low confidence: {n_contributing} contributing entries, "
            f"{coverage_pct:.0f}% coverage. Treat as a starting hint only."
        )

    if not ph_values:
        warnings.append("No structured pH data available from contributing entries.")
    if not precipitants:
        warnings.append("No precipitants identified — pdbx_details may be missing or unparseable.")

    return ConsensusResult(
        n_xray_entries=n_xray,
        n_contributing=n_contributing,
        coverage_pct=round(coverage_pct, 1),
        confidence=confidence,
        median_ph=median_ph,
        ph_range=ph_range,
        median_temp_c=median_temp_c,
        top_method=top_method,
        precipitants=precipitants,
        buffers=buffers,
        additives=additives,
        contributing_entry_ids=contributing_ids,
        warnings=warnings,
    )
