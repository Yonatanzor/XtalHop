# Parser Spec: exptl_crystal_grow Field Audit

Generated from live RCSB GraphQL audit of 60 X-ray PDB entries with confirmed pH data.

---

## Field Coverage (60 records across 60 entries)

| Field | Populated | Rate | Notes |
|-------|-----------|------|-------|
| `pH` | 60/60 | 100% | Structured float — use directly, no parsing |
| `pdbx_pH_range` | 0/60 | 0% | Always null — ignore |
| `temp` | 48/60 | 80% | Structured float in **Kelvin** — convert: C = K − 273.15 |
| `temp_details` | 7/60 | 12% | Free text, sometimes has ramp rates e.g. `'338-293 at 0.4/hr'` |
| `method` | 49/60 | 82% | Semi-structured — see method variants below |
| `details` | 0/60 | 0% | Always null — ignore |
| `pdbx_details` | 60/60 | 100% | Free text — primary parsing target for precipitants/buffers |

**Key takeaway:** Only `pH`, `temp`, `method`, and `pdbx_details` are worth parsing.
`details` and `pdbx_pH_range` are dead fields — do not attempt to parse them.

---

## pH Field

- Type: `float` (e.g. `7.0`, `5.2`, `9.0`)
- Reliable — use directly for consensus pH calculations
- **Conflict risk**: 1 entry (108M) had `pH=7.0` but `pdbx_details` contained both `'PH 9.0'` and `'pH 7.0'`. Always prefer the structured `pH` field over text extraction.

---

## Temperature Field

- Type: `float`, always in **Kelvin**
- Must convert to Celsius: `temp_C = temp_K - 273.15`
- Common values: `277 K (4°C)`, `291 K (18°C)`, `293 K (20°C)`, `295 K (22°C)`, `298 K (25°C)`
- `temp_details` occasionally has ramp notation like `'338-293 at 0.4/hr'` — ignore for consensus

---

## Method Field Variants (semi-structured, 82% populated)

Observed values — normalize to canonical form:

| Raw value | Canonical |
|-----------|-----------|
| `'VAPOR DIFFUSION, HANGING DROP'` | hanging drop |
| `'VAPOR DIFFUSION, SITTING DROP'` | sitting drop |
| `'VAPOR DIFFUSION'` | vapor diffusion |
| `'MICROBATCH'` | microbatch |

---

## pdbx_details Field — Format Variations

This is the primary parsing target. **Wildly inconsistent** across depositors.
All patterns below were observed in the 60-entry audit.

### Pattern 1: Minimal (pH only)
```
'pH 5.8'
```
Extracts: pH from text only. No precipitant.

### Pattern 2: Method-encoded, no conditions
```
'pH 7.00, VAPOR DIFFUSION, HANGING DROP'
'pH 7.00, VAPOR DIFFUSION, SITTING DROP, temperature 277.00K'
```
Extracts: pH from text (redundant with field). Temperature occasionally in K.

### Pattern 3: All-caps, structured
```
'3.0 M AMMONIUM SULFATE, 20 MM TRIS, 1MM EDTA, PH 9.0'
'3.0 M AMMONIUM SULFATE, UNBUFFERED, pH 7.0'
'10% (W/V) PEG 8000 50 MM CALCIUM ACETATE 50 MM SODIUM CACODYLATE PH 6.5'
'100 MM AMMONIUM SULPHATE, 100 BIS-TRIS, PH 6.5, 30% PEG 3350, 300 UM SILVER NITRATE'
```
Notes:
- `MM` = mM (millimolar) in all-caps records — NOT megamolar
- `UM` = µM
- Missing unit possible: `'100 BIS-TRIS'` has no `MM` — context implies mM
- British spelling: `SULPHATE` = `SULFATE`
- `(W/V)` notation for percentages

### Pattern 4: Mixed-case, screen-name prefix
```
'Index A5: 2.0 M Ammonium sulfate, 0.1 M HEPES pH 7.5. BrmaA.01375...'
'Berkeley H5: 100 mM MgCl2, 100 mM Li2SO4, 25% PEG 400. BuxeA...'
'Grid Salt Screen, C12: 3.40M Sodium Malonate, pH 5.0. TrcrB...'
'Index F6 0.2M Ammonium sulfate 0.1M Bis-Tris pH 5.5 25% PEG 3350. TrbrA...'
'Morpheous D1: 20%(v/v) PEG 500 MME, 10%(w/v) PEG 20000, 100 mM Imidazole/MES...'
'PPX D8: 0.2M Magnesium Acetate, 0.1M MES pH 6.5, 15% PEG 6000. ...'
'Berkeley B3: 30%(v/v) PEG 3350, 100mM Bis-Tris pH 6.5...'
```
Notes:
- Screen name prefix before colon can be stripped: regex `^[A-Za-z ]+[A-F0-9]+:\s*`
- Trailing junk after period (`.`): plate numbers, puck IDs, cryo info — strip from `.` onward
  - Exception: `'pH 5.8'` has no trailing junk — only strip at `.` if followed by space + capital
- `%(v/v)` and `%(w/v)` are both % concentration variants
- `100mM` (no space) = `100 mM`
- `0.2M` (no space) = `0.2 M`

### Pattern 5: Clean mixed-case, no screen name
```
'0.1 M BIS-TRIS at pH 6.5 and 28% w/v polyethylene glycol monomethyl ether 2,000.'
'11% PEG8000, 15-20% ethylene glycol, 0.1 M MES, pH 6.5, 5% 6-aminohexanoic acid'
'100 MM HEPES, pH 7.5, 30% (w/v) PEG 4000, 200 MM calcium chloride dihyrate'
'0.1 M HEPES, pH 7.0, 10-12% PEG8000'
'0.1 M MES, pH 6.5, 24% PEG8000, 5% glycerol'
'25-29% PEG4000, 0.2 M CaCl2, 0.1 M Tris pH 8.5'
'22% PEG 4,000, 200 mM Ammonium Sulfate, 100 mM Sodium Citrate pH 5.6'
'100 mM sodium citrate pH 5.2, 300 mM Na/K tartrate, and 1.4 M ammonium sulfate'
'0.1 M HEPES, pH 7.5, 0.2 M NaOAc, 22% PEG 4000'
```
Notes:
- Long-form chemical name: `'polyethylene glycol monomethyl ether 2,000'` = PEG MME 2000
- Range percentages: `'15-20%'` and `'10-12%'` — take midpoint or lower bound for consensus
- `PEG4000` = `PEG 4000` = `PEG 4,000` (comma in number!)
- `NaOAc` = sodium acetate
- `dihyrate` = typo for dihydrate — don't trust spelling

### Pattern 6: Conflict / multiple conditions
```
'3.0 M AMMONIUM SULFATE, 20 MM TRIS, 1MM EDTA, PH 9.0, pH 7.0'
```
Entry 108M has TWO pH values in the text. Prefer `pH` structured field (7.0). Do not extract pH from pdbx_details text at all — only extract precipitants/buffers from text.

---

## Critical Parsing Rules

### Rule 1: Never parse pH from pdbx_details
Use the structured `pH` float field always. It is 100% populated and structured.
Text pH extraction will fail on cases like 108M (multiple values).

### Rule 2: Never parse temperature from pdbx_details
Use the structured `temp` float field. Convert K → C. Ignore `temp_details`.

### Rule 3: Primary parsing target = precipitant/buffer identification from pdbx_details

Extract these three things from `pdbx_details`:
1. **Precipitants** (salts, PEGs, organics with concentration)
2. **Buffers** (named buffer + concentration)
3. **Additives** (small molecules at low concentration)

### Rule 4: Pre-processing pipeline for pdbx_details
1. Uppercase the entire string
2. Strip screen-name prefix: remove `^[A-Z0-9 ]+[A-F0-9]+:\s*`
3. Strip trailing metadata: truncate at first `. ` followed by alphanumeric tracking ID (e.g. `. TrcrB`, `. BrmaA`, `. plate`)
4. Normalize concentration notation:
   - `(\d+\.?\d*)M\b` (no space) → `$1 M`
   - `(\d+\.?\d*)MM\b` → `$1 MM` → treat as mM
   - `%\(W/V\)` / `%\(V/V\)` / `% W/V` → `%`
   - `(\d+)-(\d+)%` (range) → take lower bound
5. Normalize PEG names → canonical chemotype:
   - `PEG\s?(\d+)` / `PEG\s?\d+,\d+` / `POLYETHYLENE GLYCOL.*?(\d+)` → `PEG {MW}`
   - `PEG MME`, `PEG MONOMETHYL ETHER` → `PEG MME {MW}`
6. Normalize salt/buffer names (see lookup table below)

### Rule 5: Case normalization
Always uppercase before matching. All chemical name lookups should use uppercase keys.

---

## Chemotype Lookup Table (from audit + domain knowledge)

### PEG variants → canonical
| Raw | Canonical |
|-----|-----------|
| PEG 3350, PEG3350 | PEG 3350 |
| PEG 4000, PEG4000, PEG 4,000 | PEG 4000 |
| PEG 8000, PEG8000 | PEG 8000 |
| PEG 400, PEG400 | PEG 400 |
| PEG 2000, PEG2000 | PEG 2000 |
| PEG 500 MME, PEG MME 500 | PEG MME 500 |
| PEG 20000, PEG20000 | PEG 20000 |
| POLYETHYLENE GLYCOL MONOMETHYL ETHER 2000 | PEG MME 2000 |
| PEG 200, PEG200 | PEG 200 (cryo agent — may appear as cryo, not precipitant) |

### Salt/precipitant variants → canonical
| Raw | Canonical |
|-----|-----------|
| AMMONIUM SULFATE, AMMONIUM SULPHATE, (NH4)2SO4 | AMMONIUM SULFATE |
| SODIUM MALONATE, NA MALONATE | SODIUM MALONATE |
| SODIUM FORMATE, NA FORMATE | SODIUM FORMATE |
| MAGNESIUM SULFATE, MGSO4 | MAGNESIUM SULFATE |
| LITHIUM SULFATE, LI2SO4, LISO4 | LITHIUM SULFATE |
| NA/K TARTRATE, SODIUM POTASSIUM TARTRATE | Na/K TARTRATE |
| CALCIUM CHLORIDE, CACL2 | CALCIUM CHLORIDE |
| MAGNESIUM ACETATE, MG(OAC)2 | MAGNESIUM ACETATE |
| MAGNESIUM CHLORIDE, MGCL2 | MAGNESIUM CHLORIDE |
| AMMONIUM ACETATE | AMMONIUM ACETATE |

### Buffer variants → canonical
| Raw | Canonical |
|-----|-----------|
| HEPES | HEPES |
| BIS-TRIS, BIS TRIS, BISTRIS | BIS-TRIS |
| MES | MES |
| MOPS | MOPS |
| TRIS, TRIS-HCL | TRIS |
| SODIUM CITRATE, NA CITRATE | SODIUM CITRATE |
| SODIUM CACODYLATE, CACODYLATE | SODIUM CACODYLATE |
| SODIUM ACETATE, NAOAC | SODIUM ACETATE |
| IMIDAZOLE | IMIDAZOLE |
| BIS-TRIS PROPANE, BTP | BIS-TRIS PROPANE |
| CITRIC ACID | CITRIC ACID |

### Additives (low concentration, do not include in precipitant consensus)
- EDTA, DTT, BME, GLYCEROL, ETHYLENE GLYCOL, MPD, PEG 200 (cryo)
- Any compound at < 5 mM concentration
- Anything after `CRYO:` in tracking-style entries

---

## Regex Patterns for Concentration Extraction

```python
# Molar concentrations (M)
MOLAR_RE = r'(\d+\.?\d*)\s*M\b(?!M)'  # matches "0.1 M", "2.0M" but not "MM"

# Millimolar concentrations (mM)  
MILLIMOLAR_RE = r'(\d+\.?\d*)\s*[Mm][Mm]\b'  # matches "100 MM", "100mM", "20 MM"

# Percent concentrations
PERCENT_RE = r'(\d+\.?\d*(?:-\d+\.?\d*)?)\s*%(?:\s*\([wWvV]/[wWvV]\))?'

# pH from text (use only for validation against structured field)
PH_RE = r'\bpH\s+(\d+\.?\d*)'

# Temperature from text in K
TEMP_K_RE = r'temperature\s+(\d+\.?\d*)\s*[Kk]\b'
```

---

## What the `details` field is (empty — for reference)

In older PDB entries (pre-2000), crystallization conditions were stored in REMARK 280 of mmCIF files. The `details` field was intended to capture this. In modern entries the depositor uses `pdbx_details`. In practice: **`details` is always null in this API** — all condition text is in `pdbx_details`.

---

## Benchmark Targets for Condition Normalizer

Based on this audit, a correctly implemented normalizer on `pdbx_details` should achieve:
- **Precipitant identification**: ≥ 85% precision on named salts/PEGs (pattern is clear)
- **Buffer identification**: ≥ 80% precision (name variants manageable with lookup table)
- **Concentration extraction**: ≥ 75% precision (unit normalization is the main challenge)
- **pH**: 100% (use structured field, not text)
- **Temperature**: 80% (structured field, just convert K→C)

Failure mode to test explicitly: entries with screen-name prefixes + trailing tracking metadata.
Use 108M (conflicting pH values) as a known adversarial test case.
