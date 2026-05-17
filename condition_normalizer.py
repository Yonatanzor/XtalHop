"""
condition_normalizer.py — Step 0b
Parses and normalizes raw pdbx_details strings from RCSB exptl_crystal_grow records.
Built directly from the format patterns documented in parser_spec.md.
"""

import re
from dataclasses import dataclass, field


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ChemComponent:
    name: str          # canonical name
    raw_name: str      # as found in text
    concentration: float | None
    unit: str | None   # 'M', 'mM', '%', 'uM', None
    role: str          # 'precipitant', 'buffer', 'additive', 'unknown'


@dataclass
class NormalizedCondition:
    precipitants: list[ChemComponent] = field(default_factory=list)
    buffers: list[ChemComponent] = field(default_factory=list)
    additives: list[ChemComponent] = field(default_factory=list)
    method: str | None = None
    raw_cleaned: str = ""   # pdbx_details after pre-processing
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Canonical name lookup tables (from parser_spec.md)
# ---------------------------------------------------------------------------

PEG_PATTERNS = [
    # Ordered longest-match first to avoid sub-matches
    (re.compile(r'POLYETHYLENE GLYCOL MONOMETHYL ETHER\s+(\d[\d,]*)'), 'PEG MME'),
    (re.compile(r'POLYETHYLENE GLYCOL\s+(\d[\d,]*)'), 'PEG'),
    (re.compile(r'PEG\s+MME\s*(\d[\d,]*)'), 'PEG MME'),        # PEG MME 500
    (re.compile(r'PEG\s+MONOMETHYL\s+ETHER\s*(\d[\d,]*)'), 'PEG MME'),
    (re.compile(r'PEG\s*(\d[\d,]*)\s+MME\b'), 'PEG MME'),      # PEG 500 MME
    (re.compile(r'PEG\s*(\d[\d,]*)'), 'PEG'),
]

SALT_CANONICAL = {
    'AMMONIUM SULFATE': 'AMMONIUM SULFATE',
    'AMMONIUM SULPHATE': 'AMMONIUM SULFATE',
    '(NH4)2SO4': 'AMMONIUM SULFATE',
    'NH4SO4': 'AMMONIUM SULFATE',
    'SODIUM MALONATE': 'SODIUM MALONATE',
    'NA MALONATE': 'SODIUM MALONATE',
    'SODIUM FORMATE': 'SODIUM FORMATE',
    'NA FORMATE': 'SODIUM FORMATE',
    'MAGNESIUM SULFATE': 'MAGNESIUM SULFATE',
    'MAGNESIUM SULPHATE': 'MAGNESIUM SULFATE',
    'MGSO4': 'MAGNESIUM SULFATE',
    'LITHIUM SULFATE': 'LITHIUM SULFATE',
    'LI2SO4': 'LITHIUM SULFATE',
    'LISO4': 'LITHIUM SULFATE',
    'LITHIUM SULPHATE': 'LITHIUM SULFATE',
    'SODIUM POTASSIUM TARTRATE': 'Na/K TARTRATE',
    'NA/K TARTRATE': 'Na/K TARTRATE',
    'POTASSIUM SODIUM TARTRATE': 'Na/K TARTRATE',
    'CALCIUM CHLORIDE': 'CALCIUM CHLORIDE',
    'CACL2': 'CALCIUM CHLORIDE',
    'MAGNESIUM ACETATE': 'MAGNESIUM ACETATE',
    'MAGNESIUM CHLORIDE': 'MAGNESIUM CHLORIDE',
    'MGCL2': 'MAGNESIUM CHLORIDE',
    'AMMONIUM ACETATE': 'AMMONIUM ACETATE',
    'SODIUM CHLORIDE': 'SODIUM CHLORIDE',
    'NACL': 'SODIUM CHLORIDE',
    'POTASSIUM CHLORIDE': 'POTASSIUM CHLORIDE',
    'KCL': 'POTASSIUM CHLORIDE',
    'SODIUM ACETATE': 'SODIUM ACETATE',
    'NAOAC': 'SODIUM ACETATE',
    'NA ACETATE': 'SODIUM ACETATE',
    'SODIUM CITRATE': 'SODIUM CITRATE',
    'NA CITRATE': 'SODIUM CITRATE',
    'ZINC ACETATE': 'ZINC ACETATE',
    'ZINC CHLORIDE': 'ZINC CHLORIDE',
    'ZNCL2': 'ZINC CHLORIDE',
    'POTASSIUM THIOCYANATE': 'POTASSIUM THIOCYANATE',
    'K-THIOCYANATE': 'POTASSIUM THIOCYANATE',
    'KSCN': 'POTASSIUM THIOCYANATE',
    'SODIUM THIOCYANATE': 'SODIUM THIOCYANATE',
    'LITHIUM CHLORIDE': 'LITHIUM CHLORIDE',
    'LICL': 'LITHIUM CHLORIDE',
    'AMMONIUM CHLORIDE': 'AMMONIUM CHLORIDE',
    'NH4CL': 'AMMONIUM CHLORIDE',
    'SILVER NITRATE': 'SILVER NITRATE',
    'CADMIUM CHLORIDE': 'CADMIUM CHLORIDE',
}

BUFFER_CANONICAL = {
    'HEPES': 'HEPES',
    'BIS-TRIS PROPANE': 'BIS-TRIS PROPANE',
    'BTP': 'BIS-TRIS PROPANE',
    'BIS-TRIS': 'BIS-TRIS',
    'BIS TRIS': 'BIS-TRIS',
    'BISTRIS': 'BIS-TRIS',
    'MES': 'MES',
    'MOPS': 'MOPS',
    'TRIS-HCL': 'TRIS',
    'TRIS HCL': 'TRIS',
    'TRIS': 'TRIS',
    'SODIUM CACODYLATE': 'SODIUM CACODYLATE',
    'CACODYLATE': 'SODIUM CACODYLATE',
    'SODIUM CITRATE': 'SODIUM CITRATE',
    'CITRIC ACID': 'CITRIC ACID',
    'IMIDAZOLE': 'IMIDAZOLE',
    'SODIUM ACETATE': 'SODIUM ACETATE',
    'ACETATE': 'ACETATE BUFFER',
    'PHOSPHATE': 'PHOSPHATE',
    'SODIUM PHOSPHATE': 'SODIUM PHOSPHATE',
    'POTASSIUM PHOSPHATE': 'POTASSIUM PHOSPHATE',
    'PIPES': 'PIPES',
    'BICINE': 'BICINE',
    'TRICINE': 'TRICINE',
    'CHES': 'CHES',
    'CAPS': 'CAPS',
    'ACETATE/ACIDIC ACID': 'ACETATE BUFFER',
}

ADDITIVE_NAMES = {
    'EDTA', 'DTT', 'BME', 'BETA-MERCAPTOETHANOL', 'DITHIOTHREITOL',
    'GLYCEROL', 'ETHYLENE GLYCOL', 'MPD', '2-METHYL-2,4-PENTANEDIOL',
    '1,6-HEXANEDIOL', '1-BUTANOL', '1,2-PROPANEDIOL',
    '2-PROPANOL', '1,4-BUTANEDIOL', '1,3-PROPANEDIOL', '6-AMINOHEXANOIC ACID',
    'TRIMETHYLAMINE N-OXIDE', 'TMAO', 'SPERMIDINE', 'SPERMINE',
    'COENZYME A', 'COA',
}
# PEG MW values typically used as cryo/additive (not precipitant)
CRYO_PEG_MW = {'200', '300', '400'}  # PEG 400 can be precipitant too — context-dependent

METHOD_CANONICAL = {
    'VAPOR DIFFUSION, HANGING DROP': 'hanging drop',
    'VAPOR DIFFUSION, SITTING DROP': 'sitting drop',
    'VAPOR DIFFUSION': 'vapor diffusion',
    'MICROBATCH': 'microbatch',
    'BATCH': 'batch',
    'DIALYSIS': 'dialysis',
    'FREE INTERFACE DIFFUSION': 'free interface diffusion',
    'LIPIDIC CUBIC PHASE': 'lipidic cubic phase',
    'LCP': 'lipidic cubic phase',
}


# ---------------------------------------------------------------------------
# Pre-processing
# ---------------------------------------------------------------------------

# Regex for screen-name prefix: "Index A5:", "Berkeley H5:", "Grid Salt Screen, C12:"
_SCREEN_PREFIX_RE = re.compile(
    r'^(?:[A-Za-z0-9 ,\-]+?)\s+[A-H]\d+\s*:\s*', re.IGNORECASE
)
# Also: "Morpheous D1:" style
_SCREEN_PREFIX2_RE = re.compile(
    r'^(?:[A-Za-z]+\s+[A-H]\d+)\s*:\s*', re.IGNORECASE
)

# Trailing tracking metadata: ". BrmaA.01375...", ". TrcrB...", ". plate ...", ". Puck:"
_TRACKING_RE = re.compile(
    r'\.\s+(?:[A-Z][a-z]{3}[A-Z]\.\d{5}|[Pp]late\s|[Pp]uck\s*:|[Pp]uck:)',
    re.IGNORECASE
)


def preprocess(raw: str) -> str:
    """Clean pdbx_details text before extraction."""
    if not raw:
        return ""

    s = raw.strip()

    # Strip screen-name prefix
    m = _SCREEN_PREFIX_RE.match(s)
    if m:
        s = s[m.end():]
    else:
        m2 = _SCREEN_PREFIX2_RE.match(s)
        if m2:
            s = s[m2.end():]

    # Strip trailing tracking metadata
    m = _TRACKING_RE.search(s)
    if m:
        s = s[:m.start()]

    # Uppercase for uniform matching
    s = s.upper()

    # Normalize concentration notation: "0.2M " → "0.2 M "
    s = re.sub(r'(\d+\.?\d*)(M)\b(?!M)', r'\1 M ', s)
    s = re.sub(r'(\d+\.?\d*)(MM)\b', r'\1 MM ', s)

    # Normalize percent: "%(w/v)", "% (w/v)", "%(v/v)" → "%"
    s = re.sub(r'%\s*\([WwVv]/[WwVv]\)', '%', s)
    s = re.sub(r'%\s+[WwVv]/[WwVv]', '%', s)

    # Normalize ranges: "15-20%" → "15%"  (take lower bound)
    s = re.sub(r'(\d+\.?\d*)-\d+\.?\d*\s*%', r'\1%', s)

    # Remove comma in PEG molecular weights: "PEG 4,000" → "PEG 4000"
    s = re.sub(r'(PEG\s*(?:MME\s*)?)(\d+),(\d+)', r'\1\2\3', s)

    # Collapse multiple spaces
    s = re.sub(r'\s+', ' ', s).strip()

    return s


# ---------------------------------------------------------------------------
# Concentration extraction
# ---------------------------------------------------------------------------

_MOLAR_RE = re.compile(r'(\d+\.?\d*)\s*M\b(?!\s*M)')
_MILLIMOLAR_RE = re.compile(r'(\d+\.?\d*)\s*MM\b')
_MICROMOLAR_RE = re.compile(r'(\d+\.?\d*)\s*(?:UM|MICROMOLAR)\b')
_PERCENT_RE = re.compile(r'(\d+\.?\d*)\s*%')


def extract_concentration(text_segment: str, last: bool = False) -> tuple[float | None, str | None]:
    """
    Extract a concentration + unit from a text segment.
    last=True returns the rightmost match (closest to a following chemical name).
    last=False returns the leftmost match.
    """
    all_matches: list[tuple[int, float, str]] = []
    for m in _MOLAR_RE.finditer(text_segment):
        all_matches.append((m.start(), float(m.group(1)), 'M'))
    for m in _MILLIMOLAR_RE.finditer(text_segment):
        all_matches.append((m.start(), float(m.group(1)), 'mM'))
    for m in _MICROMOLAR_RE.finditer(text_segment):
        all_matches.append((m.start(), float(m.group(1)), 'uM'))
    for m in _PERCENT_RE.finditer(text_segment):
        all_matches.append((m.start(), float(m.group(1)), '%'))
    if not all_matches:
        return None, None
    all_matches.sort(key=lambda x: x[0])
    _, val, unit = all_matches[-1] if last else all_matches[0]
    return val, unit


# ---------------------------------------------------------------------------
# PEG extraction
# ---------------------------------------------------------------------------

def extract_pegs(text: str) -> list[ChemComponent]:
    """Find all PEG components in text."""
    found = []
    for pattern, peg_type in PEG_PATTERNS:
        for m in pattern.finditer(text):
            mw = m.group(1).replace(',', '')
            canonical = f'{peg_type} {mw}'
            conc, unit = _concentration_in_segment(text, m.start(), m.end())
            found.append(ChemComponent(
                name=canonical,
                raw_name=m.group(0),
                concentration=conc,
                unit=unit,
                role='precipitant',
            ))
            # Mask matched region so inner patterns don't double-match
            text = text[:m.start()] + ' ' * (m.end() - m.start()) + text[m.end():]
    return found


# ---------------------------------------------------------------------------
# Salt/buffer extraction using token-based lookup
# ---------------------------------------------------------------------------

def _find_canonical(text: str, lookup: dict[str, str]) -> list[tuple[str, str, int, int]]:
    """
    Find all occurrences of lookup keys in text (longest match first).
    Returns list of (raw_match, canonical, start, end).
    """
    results = []
    keys_sorted = sorted(lookup.keys(), key=len, reverse=True)
    masked = text
    for key in keys_sorted:
        start = 0
        while True:
            idx = masked.find(key, start)
            if idx == -1:
                break
            # Confirm it's a word boundary (not mid-word)
            before = masked[idx-1] if idx > 0 else ' '
            after = masked[idx + len(key)] if idx + len(key) < len(masked) else ' '
            if not before.isalpha() and not after.isalpha():
                results.append((key, lookup[key], idx, idx + len(key)))
                # Mask to avoid sub-matches
                masked = masked[:idx] + ' ' * len(key) + masked[idx + len(key):]
            start = idx + 1
    return results


def _concentration_in_segment(text: str, match_start: int, match_end: int) -> tuple[float | None, str | None]:
    """
    Find concentration for a chemical match by scanning its comma-delimited segment.
    Prefers looking before the match within the same segment, then after.
    """
    # Find segment boundaries (comma or semicolon)
    seg_start = match_start
    for ch in reversed(text[:match_start]):
        if ch in ',;':
            break
        seg_start -= 1
    seg_start = max(0, seg_start)

    seg_end = match_end
    for i, ch in enumerate(text[match_end:]):
        if ch in ',;':
            break
        seg_end += 1
    seg_end = min(len(text), seg_end)

    segment = text[seg_start:seg_end]
    # Prefer concentration before the match within the segment
    before = text[seg_start:match_start]
    conc, unit = extract_concentration(before, last=True)  # closest = rightmost in prefix
    if conc is None:
        after = text[match_end:seg_end]
        conc, unit = extract_concentration(after, last=False)
    return conc, unit


# Regex: captures "X.X M/MM/% Chemical Name Words" within a single comma segment.
# Name capture stops at comma, semicolon, digit+unit (new component), or pH marker.
_CONC_CHEM_RE = re.compile(
    r'(\d+\.?\d*)\s*(M\b(?!M)|MM\b|%)\s+'
    r'((?:[A-Z][A-Z0-9\-/\(\)]*(?:\s+(?![\d,;]|PH\s))?){1,6})',
    re.IGNORECASE,
)

# Tokens that are NOT chemical names (stop words for name extraction)
_CHEM_STOP = {
    'PH', 'AT', 'AND', 'WITH', 'IN', 'FROM', 'TO', 'OR',
    'VAPOR', 'DIFFUSION', 'HANGING', 'SITTING', 'DROP', 'BATCH',
    'TEMPERATURE', 'SEEDED', 'DIRECT', 'CRYO',
}


def _clean_chem_name(raw_name: str) -> str:
    """Strip trailing stop words from extracted chemical name."""
    tokens = raw_name.strip().split()
    while tokens and tokens[-1].upper() in _CHEM_STOP:
        tokens.pop()
    return ' '.join(tokens).strip().rstrip(',;').strip()


def _classify_name(name_upper: str) -> tuple[str, str]:
    """
    Return (canonical_name, role) for a chemical name.
    Checks buffer table first, then salt table, then falls back to 'unknown'.
    """
    # Exact match first (longest-key first)
    for key in sorted(BUFFER_CANONICAL, key=len, reverse=True):
        if key in name_upper:
            return BUFFER_CANONICAL[key], 'buffer'
    for key in sorted(SALT_CANONICAL, key=len, reverse=True):
        if key in name_upper:
            return SALT_CANONICAL[key], 'precipitant'
    # No match — use the extracted name as-is (title-cased), role = precipitant
    return name_upper.title(), 'precipitant'


def extract_salts_and_buffers(text: str) -> tuple[list[ChemComponent], list[ChemComponent]]:
    """
    Pattern-first extraction: find all {concentration} {chemical_name} in text.
    Classifies using lookup tables; unknown chemicals default to precipitant.
    """
    salts: list[ChemComponent] = []
    buffers: list[ChemComponent] = []
    seen_names: set[str] = set()

    for m in _CONC_CHEM_RE.finditer(text):
        conc_str, unit_raw, chem_raw = m.group(1), m.group(2), m.group(3)
        conc = float(conc_str)
        unit = 'mM' if unit_raw.upper() == 'MM' else unit_raw.upper().rstrip('B')  # M, %

        chem_clean = _clean_chem_name(chem_raw)
        if not chem_clean or len(chem_clean) < 2:
            continue
        chem_upper = chem_clean.upper()

        # Skip duplicates
        if chem_upper in seen_names:
            continue
        seen_names.add(chem_upper)

        # Skip obvious noise tokens
        if chem_upper in _CHEM_STOP:
            continue

        canonical, role = _classify_name(chem_upper)

        comp = ChemComponent(
            name=canonical, raw_name=chem_clean,
            concentration=conc, unit=unit, role=role,
        )
        if role == 'buffer':
            buffers.append(comp)
        else:
            salts.append(comp)

    return salts, buffers


# ---------------------------------------------------------------------------
# Main normalizer
# ---------------------------------------------------------------------------

def normalize(pdbx_details: str | None) -> NormalizedCondition:
    """
    Parse a raw pdbx_details string into structured condition components.
    pH and temperature are NOT extracted here — use the structured API fields.
    """
    result = NormalizedCondition()
    if not pdbx_details:
        result.warnings.append("pdbx_details is null")
        return result

    cleaned = preprocess(pdbx_details)
    result.raw_cleaned = cleaned

    # PEG components first (before salt lookup can steal numbers)
    pegs = extract_pegs(cleaned)

    # Build a masked version for salt/buffer search (PEG regions masked)
    masked = cleaned
    for p in pegs:
        # Re-find in masked string (approximate — use name match)
        idx = masked.find(p.raw_name.upper())
        if idx != -1:
            masked = masked[:idx] + ' ' * len(p.raw_name) + masked[idx + len(p.raw_name):]

    salts, buffers = extract_salts_and_buffers(masked)

    # Tag role: if a buffer name also appears in salts (e.g. SODIUM CITRATE), keep as buffer
    buffer_names = {b.name for b in buffers}
    precipitants = [s for s in salts if s.name not in buffer_names] + pegs

    # Classify additives: known additive names or concentration < 5 mM
    final_precipitants = []
    final_additives = []
    for c in precipitants:
        upper_name = c.name.upper()
        if upper_name in ADDITIVE_NAMES or any(a in upper_name for a in ADDITIVE_NAMES):
            c.role = 'additive'
            final_additives.append(c)
        elif c.unit == 'mM' and c.concentration is not None and c.concentration < 5:
            c.role = 'additive'
            final_additives.append(c)
        else:
            final_precipitants.append(c)

    for b in buffers:
        if b.name.upper() in ADDITIVE_NAMES:
            b.role = 'additive'
            final_additives.append(b)

    result.precipitants = final_precipitants
    result.buffers = [b for b in buffers if b.name.upper() not in ADDITIVE_NAMES]
    result.additives = final_additives

    if not result.precipitants and not result.buffers:
        result.warnings.append("no components extracted")

    return result


# ---------------------------------------------------------------------------
# Human-readable summary
# ---------------------------------------------------------------------------

def summarize(cond: NormalizedCondition) -> str:
    lines = []
    if cond.precipitants:
        parts = []
        for p in cond.precipitants:
            if p.concentration is not None:
                parts.append(f"{p.concentration} {p.unit} {p.name}")
            else:
                parts.append(p.name)
        lines.append("Precipitants: " + " | ".join(parts))
    if cond.buffers:
        parts = []
        for b in cond.buffers:
            if b.concentration is not None:
                parts.append(f"{b.concentration} {b.unit} {b.name}")
            else:
                parts.append(b.name)
        lines.append("Buffers:      " + " | ".join(parts))
    if cond.additives:
        lines.append("Additives:    " + " | ".join(a.name for a in cond.additives))
    if cond.warnings:
        lines.append("Warnings:     " + "; ".join(cond.warnings))
    return "\n".join(lines) if lines else "(no components extracted)"


# ---------------------------------------------------------------------------
# Quick smoke test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    TESTS = [
        "3.0 M AMMONIUM SULFATE, 20 MM TRIS, 1MM EDTA, PH 9.0",
        "Index A5: 2.0 M Ammonium sulfate, 0.1 M HEPES pH 7.5. BrmaA.01375.b.B1.PS01744 at 31 mg/mL.",
        "Berkeley H5: 100 mM MgCl2, 100 mM Li2SO4, 25% PEG 400. BuxeA.00036.a.B2.PW39468 at 23.8 mg/mL.",
        "11% PEG8000, 15-20% ethylene glycol, 0.1 M MES, pH 6.5, 5% 6-aminohexanoic acid",
        "0.1 M BIS-TRIS at pH 6.5 and 28% w/v polyethylene glycol monomethyl ether 2,000.",
        "100 MM HEPES, pH 7.5, 30% (w/v) PEG 4000, 200 MM calcium chloride dihyrate",
        "100 mM MOPS, 1.25 M magnesium sulfate",
        "22% PEG 4,000, 200 mM Ammonium Sulfate, 100 mM Sodium Citrate pH 5.6",
        "Grid Salt Screen, C12: 3.40M Sodium Malonate, pH 5.0. TrcrB.01480.a.WW4.PS38791 at 16.9 mg/mL.",
        "pH 5.8",
        "0.1 M Bis-Tris propane pH 7.5, 0.2 M potassium thiocyanate, 12% (w/v) PEG3350",
        "Morpheous D1: 20%(v/v) PEG 500 MME, 10%(w/v) PEG 20000, 100 mM Imidazole/MES, pH 6.5. BaboA.18322.a.B2.PW39407 at 19.1 mg/mL.",
        "25 mM Citric acid, 75 mM Bis-tris propane, pH 7.4, 20 mM Cadmium chloride, 25% PEG400, seeded",
    ]

    for raw in TESTS:
        print(f"\nInput:  {raw!r}")
        cond = normalize(raw)
        print(f"Cleaned:{cond.raw_cleaned!r}")
        print(summarize(cond))
        print()
