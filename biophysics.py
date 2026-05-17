"""
biophysics.py — Module 4
Local biophysical risk panel via Biopython ProteinAnalysis.
No external API calls — runs entirely on the input sequence.
"""

from dataclasses import dataclass, field
from Bio.SeqUtils.ProtParam import ProteinAnalysis


@dataclass
class BiophysicsResult:
    sequence_length: int
    pi: float                       # isoelectric point
    instability_index: float        # >40 = unstable by Guruprasad 1990
    gravy: float                    # hydropathicity; negative = hydrophilic
    molecular_weight_da: float
    aa_composition: dict[str, float]  # fraction 0.0–1.0 per residue
    # Risk flags
    risks: list[str] = field(default_factory=list)
    # Human-readable summary
    summary: str = ""


# Thresholds derived from crystallization literature
_INSTABILITY_WARN = 40.0     # Guruprasad cutoff
_GRAVY_HYDROPHOBIC = 0.0     # positive → likely membrane / aggregation-prone
_CYS_HIGH_PCT = 5.0          # > 5% Cys → disulfide complexity
_PRO_HIGH_PCT = 8.0          # > 8% Pro → rigid backbone, poor crystal contacts
_GLY_HIGH_PCT = 12.0         # > 12% Gly → disordered regions


def analyze(sequence: str) -> BiophysicsResult:
    """
    Run the biophysical risk panel on a protein sequence.

    Handles FASTA headers and strips whitespace. Replaces ambiguous
    residues (B→N, Z→Q, X→A, U→C) so ProteinAnalysis doesn't crash.
    """
    lines = sequence.strip().splitlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    # Sanitize ambiguous / non-standard residues
    seq = (seq.upper()
              .replace("B", "N")
              .replace("Z", "Q")
              .replace("X", "A")
              .replace("U", "C"))

    if not seq:
        raise ValueError("Empty sequence after parsing")

    pa = ProteinAnalysis(seq)

    pi = round(pa.isoelectric_point(), 2)
    instability = round(pa.instability_index(), 2)
    gravy = round(pa.gravy(), 3)
    mw = round(pa.molecular_weight(), 1)
    # amino_acids_percent returns percentages in 0–100 range
    comp = {aa: round(pct, 2) for aa, pct in pa.amino_acids_percent.items()}

    risks: list[str] = []

    if instability > _INSTABILITY_WARN:
        risks.append(
            f"Instability index {instability:.1f} > {_INSTABILITY_WARN} "
            f"— protein may be unstable in solution"
        )

    if gravy > _GRAVY_HYDROPHOBIC:
        risks.append(
            f"GRAVY {gravy:+.3f} (positive) — hydrophobic character; "
            f"aggregation or solubility issues possible"
        )

    cys_pct = comp.get("C", 0.0)
    if cys_pct > _CYS_HIGH_PCT:
        risks.append(
            f"High cysteine content ({cys_pct:.1f}%) — disulfide complexity; "
            f"consider reducing agent optimization"
        )

    pro_pct = comp.get("P", 0.0)
    if pro_pct > _PRO_HIGH_PCT:
        risks.append(
            f"High proline content ({pro_pct:.1f}%) — rigid backbone; "
            f"may reduce crystal contact surface"
        )

    gly_pct = comp.get("G", 0.0)
    if gly_pct > _GLY_HIGH_PCT:
        risks.append(
            f"High glycine content ({gly_pct:.1f}%) — flexible backbone; "
            f"disordered regions may impair crystallization"
        )

    pi_flag = ""
    if pi < 5.0:
        pi_flag = " (acidic — consider pH 4–6 buffers)"
    elif pi > 9.0:
        pi_flag = " (basic — consider pH 7–9 buffers)"
    else:
        pi_flag = " (near-neutral)"

    summary = (
        f"Length: {len(seq)} aa | "
        f"pI: {pi}{pi_flag} | "
        f"MW: {mw/1000:.1f} kDa | "
        f"Instability: {instability:.1f} | "
        f"GRAVY: {gravy:+.3f}"
    )
    if risks:
        summary += f" | Risks: {len(risks)}"

    return BiophysicsResult(
        sequence_length=len(seq),
        pi=pi,
        instability_index=instability,
        gravy=gravy,
        molecular_weight_da=mw,
        aa_composition=comp,
        risks=risks,
        summary=summary,
    )


# ---------------------------------------------------------------------------
# CLI test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    LYSOZYME = (
        "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRT"
        "PGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
    )

    print("=" * 60)
    print("biophysics.py — live test with hen lysozyme")
    print("=" * 60)

    result = analyze(LYSOZYME)

    print(f"\n{result.summary}")
    print(f"\nAmino acid composition (top 8 by frequency):")
    top = sorted(result.aa_composition.items(), key=lambda x: -x[1])[:8]
    for aa, pct in top:
        bar = "#" * int(pct * 2)   # scale 0–100% → 0–200 chars
        print(f"  {aa}: {pct:5.1f}%  {bar}")

    print(f"\nRisk flags ({len(result.risks)}):")
    if result.risks:
        for r in result.risks:
            print(f"  - {r}")
    else:
        print("  None — good crystallization candidate")
