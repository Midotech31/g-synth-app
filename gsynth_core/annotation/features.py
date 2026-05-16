"""Auto-annotation of common plasmid features.

A curated database of promoters, ORIs, resistance markers, affinity tags
and common cloning sequences. Features are detected by partial-match
(80% identity by default over the feature length) against the input.

For finer-grained annotation use pLannotate or SnapGene's database.
"""
from __future__ import annotations
from dataclasses import dataclass

from gsynth_core.io.records import Feature, FeatureLocation, SeqRecord
from gsynth_core.sequence.ops import reverse_complement


@dataclass(frozen=True, slots=True)
class FeatureHit:
    name: str
    type: str
    start: int
    end: int
    strand: int
    identity: float


# Sequences (5'→3' on the SENSE strand of the feature).
KNOWN_FEATURES: dict[str, tuple[str, str]] = {
    # promoters
    "T7_promoter":           ("promoter", "TAATACGACTCACTATAGGG"),
    "T3_promoter":           ("promoter", "ATTAACCCTCACTAAAGGGA"),
    "SP6_promoter":          ("promoter", "ATTTAGGTGACACTATAG"),
    "lac_promoter":          ("promoter", "TTTACACTTTATGCTTCCGGCTCG"),
    "CMV_promoter":          ("promoter", "CGCAAATGGGCGGTAGGCGTGTACGGTGGGAGGTCTATATAAGCAGAGCT"),
    "EF1a_promoter":         ("promoter", "GCTCCGGTGCCCGTCAGTGGGCAGAGCGCACATCGCCCACAGTCCCCGAG"),
    "SV40_promoter":         ("promoter", "GTGTGGAATGTGTGTCAGTTAGGGTGTGGAAAGTCCCCAGGCTCCCCAGCAGGCAGA"),
    "U6_promoter":           ("promoter", "GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCT"),

    # terminators
    "rrnB_T1":               ("terminator", "AAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTT"),
    "T7_terminator":         ("terminator", "CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG"),
    "bGH_polyA":             ("polyA_signal", "CTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCC"),
    "SV40_polyA":            ("polyA_signal", "AACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACAAATAAAGCATTTT"),

    # ribosome binding sites
    "RBS_strong":            ("RBS", "AAAGAGGAGAAATA"),
    "RBS_weak":              ("RBS", "GAAGGAGATATA"),
    "Kozak":                 ("misc_feature", "GCCGCCACCATGG"),

    # affinity tags (DNA)
    "6xHis_tag":             ("CDS", "CATCATCATCATCATCAT"),
    "FLAG_tag":              ("CDS", "GATTACAAGGATGACGATGACAAG"),
    "HA_tag":                ("CDS", "TATCCTTATGACGTTCCAGATTACGCT"),
    "Myc_tag":               ("CDS", "GAACAAAAACTCATCTCAGAAGAGGATCTG"),
    "Strep_II_tag":          ("CDS", "TGGAGCCACCCGCAGTTCGAAAAA"),
    "MBP_tag":               ("CDS", "ATGAAAATAAAAACAGGTGCACGCATCCTC"),
    "GST_tag":               ("CDS", "ATGTCCCCTATACTAGGTTATTGG"),

    # protease sites
    "TEV_site":              ("misc_feature", "GAAAACCTGTATTTTCAGGGC"),
    "PreScission_site":      ("misc_feature", "CTGGAAGTGCTGTTCCAGGGCCCA"),
    "Thrombin_site":         ("misc_feature", "CTGGTGCCGCGTGGTTCT"),
    "Factor_Xa_site":        ("misc_feature", "ATCGAAGGTCGT"),
    "Enterokinase_site":     ("misc_feature", "GACGACGACGACAAG"),

    # origins of replication
    "pUC_ori":               ("rep_origin", "TTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACC"),
    "pBR322_ori":            ("rep_origin", "AGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGAC"),
    "p15A_ori":              ("rep_origin", "GACGTCAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTT"),
    "f1_ori":                ("rep_origin", "ACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAG"),
    "SV40_ori":              ("rep_origin", "TGAGGCGGAAAGAACCAGCTGTGGAATGTGTGTCAGTTAGGGTGTGGAAAGTCCC"),
    "ColE1_ori":             ("rep_origin", "GACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCC"),

    # resistance markers (start of CDS, sufficient to identify)
    "AmpR_CDS":              ("CDS", "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATT"),
    "KanR_CDS":              ("CDS", "ATGATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGG"),
    "CmR_CDS":               ("CDS", "ATGGAGAAAAAAATCACTGGATATACCACCGTTGATATATCCC"),
    "TetR_CDS":              ("CDS", "ATGAAATCTAACAATGCGCTCATCGTCATCCTCGGCACCGTC"),

    # common reporters
    "GFP_CDS":               ("CDS", "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTC"),
    "RFP_CDS":               ("CDS", "ATGGCCTCCTCCGAGGACGTCATCAAGGAGTTCATGCGCTTCAAGGTG"),
    "luciferase_CDS":        ("CDS", "ATGGAAGACGCCAAAAACATAAAGAAAGGCCCGGCGCC"),
}


def annotate_features(record: SeqRecord, *,
                      min_identity: float = 0.85,
                      append_to_record: bool = True) -> list[FeatureHit]:
    """
    Scan a sequence for known plasmid features.

    Args:
        record: input :class:`SeqRecord`.
        min_identity: minimum fractional identity (default 0.85).
        append_to_record: if True, add detected features to `record.features`.

    Returns:
        List of :class:`FeatureHit` sorted by start position.
    """
    seq = record.sequence.upper()
    hits: list[FeatureHit] = []
    if not seq:
        return hits

    for name, (ftype, pattern) in KNOWN_FEATURES.items():
        pat = pattern.upper()
        for strand, target in ((1, pat), (-1, reverse_complement(pat))):
            L = len(target)
            if L == 0 or L > len(seq):
                continue
            best_pos = -1
            best_matches = -1
            # First try exact match
            i = seq.find(target)
            if i >= 0:
                best_pos = i
                best_matches = L
            else:
                # Fuzzy: scan and count matches
                threshold = int(L * min_identity)
                for i in range(0, len(seq) - L + 1):
                    matches = sum(1 for a, b in zip(seq[i:i+L], target) if a == b)
                    if matches >= threshold and matches > best_matches:
                        best_matches = matches; best_pos = i
            if best_pos >= 0:
                identity = best_matches / L
                hits.append(FeatureHit(
                    name=name, type=ftype,
                    start=best_pos, end=best_pos + L,
                    strand=strand, identity=identity,
                ))
                if append_to_record:
                    record.add_feature(Feature(
                        type=ftype,
                        location=FeatureLocation(start=best_pos, end=best_pos + L, strand=strand),
                        qualifiers={"label": name, "identity": f"{identity:.2%}"},
                    ))

    hits.sort(key=lambda h: (h.start, h.name))
    return hits
