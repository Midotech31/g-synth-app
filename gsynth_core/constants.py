"""Single source of truth for biological constants."""
from __future__ import annotations
from typing import Final, Mapping

DNA_BASES: Final[frozenset[str]] = frozenset("ACGT")

IUPAC_DNA: Final[Mapping[str, frozenset[str]]] = {
    "A": frozenset("A"), "C": frozenset("C"), "G": frozenset("G"), "T": frozenset("T"),
    "U": frozenset("T"),
    "R": frozenset("AG"), "Y": frozenset("CT"), "S": frozenset("GC"), "W": frozenset("AT"),
    "K": frozenset("GT"), "M": frozenset("AC"),
    "B": frozenset("CGT"), "D": frozenset("AGT"), "H": frozenset("ACT"), "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}
IUPAC_DNA_ALPHABET: Final[str] = "".join(IUPAC_DNA.keys())

IUPAC_COMPLEMENT: Final[Mapping[str, str]] = {
    "A": "T", "T": "A", "U": "A", "G": "C", "C": "G",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}
_C_UP = str.maketrans("".join(IUPAC_COMPLEMENT), "".join(IUPAC_COMPLEMENT.values()))
_C_LO = str.maketrans("".join(IUPAC_COMPLEMENT).lower(), "".join(IUPAC_COMPLEMENT.values()).lower())
COMPLEMENT_TABLE: Final[dict[int, int]] = {**_C_UP, **_C_LO}

GENETIC_CODE: Final[Mapping[str, str]] = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}
START_CODONS_STANDARD: Final[frozenset[str]] = frozenset({"ATG"})
START_CODONS_PROKARYOTIC: Final[frozenset[str]] = frozenset({"ATG", "GTG", "TTG"})
STOP_CODONS: Final[frozenset[str]] = frozenset({"TAA", "TAG", "TGA"})

AA_1_TO_3: Final[Mapping[str, str]] = {
    "A":"Ala","R":"Arg","N":"Asn","D":"Asp","C":"Cys","Q":"Gln","E":"Glu","G":"Gly",
    "H":"His","I":"Ile","L":"Leu","K":"Lys","M":"Met","F":"Phe","P":"Pro","S":"Ser",
    "T":"Thr","W":"Trp","Y":"Tyr","V":"Val","*":"Stop","X":"Xaa",
}
AA_3_TO_1: Final[Mapping[str, str]] = {v.upper(): k for k, v in AA_1_TO_3.items()}
STANDARD_AA: Final[str] = "ACDEFGHIKLMNPQRSTVWY"

AA_MONOISOTOPIC_MASS: Final[Mapping[str, float]] = {
    "A":71.03711,"R":156.10111,"N":114.04293,"D":115.02694,"C":103.00919,
    "Q":128.05858,"E":129.04259,"G":57.02146,"H":137.05891,"I":113.08406,
    "L":113.08406,"K":128.09496,"M":131.04049,"F":147.06841,"P":97.05276,
    "S":87.03203,"T":101.04768,"W":186.07931,"Y":163.06333,"V":99.06841,
}
WATER_MASS: Final[float] = 18.01056

KYTE_DOOLITTLE: Final[Mapping[str, float]] = {
    "A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,
    "K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,
    "T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3,
}

PKA: Final[Mapping[str, float]] = {
    "N_term":8.6,"C_term":3.6,"K":10.8,"R":12.5,"H":6.5,"D":3.9,"E":4.1,"C":8.5,"Y":10.1,
}

CodonTable = Mapping[str, Mapping[str, float]]

CODON_TABLE_ECOLI: Final[CodonTable] = {
    "F":{"TTT":0.58,"TTC":0.42},"L":{"TTA":0.14,"TTG":0.13,"CTT":0.10,"CTC":0.10,"CTA":0.04,"CTG":0.49},
    "I":{"ATT":0.49,"ATC":0.39,"ATA":0.11},"M":{"ATG":1.00},
    "V":{"GTT":0.28,"GTC":0.20,"GTA":0.17,"GTG":0.35},
    "S":{"TCT":0.17,"TCC":0.15,"TCA":0.14,"TCG":0.14,"AGT":0.16,"AGC":0.25},
    "P":{"CCT":0.18,"CCC":0.13,"CCA":0.20,"CCG":0.49},
    "T":{"ACT":0.19,"ACC":0.40,"ACA":0.17,"ACG":0.25},
    "A":{"GCT":0.18,"GCC":0.26,"GCA":0.23,"GCG":0.33},
    "Y":{"TAT":0.59,"TAC":0.41},"*":{"TAA":0.61,"TAG":0.09,"TGA":0.30},
    "H":{"CAT":0.57,"CAC":0.43},"Q":{"CAA":0.34,"CAG":0.66},
    "N":{"AAT":0.49,"AAC":0.51},"K":{"AAA":0.74,"AAG":0.26},
    "D":{"GAT":0.63,"GAC":0.37},"E":{"GAA":0.68,"GAG":0.32},
    "C":{"TGT":0.46,"TGC":0.54},"W":{"TGG":1.00},
    "R":{"CGT":0.36,"CGC":0.36,"CGA":0.07,"CGG":0.11,"AGA":0.07,"AGG":0.04},
    "G":{"GGT":0.34,"GGC":0.37,"GGA":0.13,"GGG":0.16},
}

CODON_TABLE_HUMAN: Final[CodonTable] = {
    "F":{"TTT":0.46,"TTC":0.54},"L":{"TTA":0.08,"TTG":0.13,"CTT":0.13,"CTC":0.20,"CTA":0.07,"CTG":0.39},
    "I":{"ATT":0.36,"ATC":0.47,"ATA":0.17},"M":{"ATG":1.00},
    "V":{"GTT":0.18,"GTC":0.24,"GTA":0.12,"GTG":0.46},
    "S":{"TCT":0.19,"TCC":0.22,"TCA":0.15,"TCG":0.05,"AGT":0.15,"AGC":0.24},
    "P":{"CCT":0.29,"CCC":0.32,"CCA":0.28,"CCG":0.11},
    "T":{"ACT":0.25,"ACC":0.36,"ACA":0.28,"ACG":0.11},
    "A":{"GCT":0.27,"GCC":0.40,"GCA":0.23,"GCG":0.11},
    "Y":{"TAT":0.44,"TAC":0.56},"*":{"TAA":0.30,"TAG":0.24,"TGA":0.46},
    "H":{"CAT":0.42,"CAC":0.58},"Q":{"CAA":0.27,"CAG":0.73},
    "N":{"AAT":0.47,"AAC":0.53},"K":{"AAA":0.43,"AAG":0.57},
    "D":{"GAT":0.46,"GAC":0.54},"E":{"GAA":0.42,"GAG":0.58},
    "C":{"TGT":0.46,"TGC":0.54},"W":{"TGG":1.00},
    "R":{"CGT":0.08,"CGC":0.18,"CGA":0.11,"CGG":0.20,"AGA":0.21,"AGG":0.21},
    "G":{"GGT":0.16,"GGC":0.33,"GGA":0.26,"GGG":0.25},
}

CODON_TABLE_YEAST: Final[CodonTable] = {
    "F":{"TTT":0.59,"TTC":0.41},"L":{"TTA":0.28,"TTG":0.29,"CTT":0.13,"CTC":0.06,"CTA":0.14,"CTG":0.11},
    "I":{"ATT":0.46,"ATC":0.26,"ATA":0.27},"M":{"ATG":1.00},
    "V":{"GTT":0.39,"GTC":0.21,"GTA":0.21,"GTG":0.19},
    "S":{"TCT":0.26,"TCC":0.16,"TCA":0.21,"TCG":0.10,"AGT":0.16,"AGC":0.11},
    "P":{"CCT":0.31,"CCC":0.15,"CCA":0.42,"CCG":0.12},
    "T":{"ACT":0.35,"ACC":0.22,"ACA":0.30,"ACG":0.13},
    "A":{"GCT":0.38,"GCC":0.22,"GCA":0.29,"GCG":0.11},
    "Y":{"TAT":0.56,"TAC":0.44},"*":{"TAA":0.47,"TAG":0.23,"TGA":0.30},
    "H":{"CAT":0.64,"CAC":0.36},"Q":{"CAA":0.69,"CAG":0.31},
    "N":{"AAT":0.59,"AAC":0.41},"K":{"AAA":0.58,"AAG":0.42},
    "D":{"GAT":0.65,"GAC":0.35},"E":{"GAA":0.71,"GAG":0.29},
    "C":{"TGT":0.63,"TGC":0.37},"W":{"TGG":1.00},
    "R":{"CGT":0.14,"CGC":0.06,"CGA":0.07,"CGG":0.04,"AGA":0.48,"AGG":0.21},
    "G":{"GGT":0.47,"GGC":0.19,"GGA":0.22,"GGG":0.12},
}

CODON_TABLE_PICHIA: Final[CodonTable] = {
    "F":{"TTT":0.54,"TTC":0.46},"L":{"TTA":0.16,"TTG":0.33,"CTT":0.17,"CTC":0.08,"CTA":0.11,"CTG":0.15},
    "I":{"ATT":0.52,"ATC":0.28,"ATA":0.20},"M":{"ATG":1.00},
    "V":{"GTT":0.42,"GTC":0.23,"GTA":0.16,"GTG":0.19},
    "S":{"TCT":0.29,"TCC":0.20,"TCA":0.19,"TCG":0.09,"AGT":0.14,"AGC":0.09},
    "P":{"CCT":0.39,"CCC":0.15,"CCA":0.33,"CCG":0.13},
    "T":{"ACT":0.40,"ACC":0.24,"ACA":0.23,"ACG":0.13},
    "A":{"GCT":0.45,"GCC":0.24,"GCA":0.22,"GCG":0.09},
    "Y":{"TAT":0.55,"TAC":0.45},"*":{"TAA":0.52,"TAG":0.21,"TGA":0.27},
    "H":{"CAT":0.56,"CAC":0.44},"Q":{"CAA":0.62,"CAG":0.38},
    "N":{"AAT":0.50,"AAC":0.50},"K":{"AAA":0.53,"AAG":0.47},
    "D":{"GAT":0.64,"GAC":0.36},"E":{"GAA":0.59,"GAG":0.41},
    "C":{"TGT":0.64,"TGC":0.36},"W":{"TGG":1.00},
    "R":{"CGT":0.16,"CGC":0.05,"CGA":0.09,"CGG":0.04,"AGA":0.47,"AGG":0.19},
    "G":{"GGT":0.43,"GGC":0.16,"GGA":0.28,"GGG":0.13},
}

CODON_TABLE_CHO: Final[CodonTable] = {
    "F":{"TTT":0.43,"TTC":0.57},"L":{"TTA":0.06,"TTG":0.13,"CTT":0.12,"CTC":0.21,"CTA":0.07,"CTG":0.41},
    "I":{"ATT":0.32,"ATC":0.51,"ATA":0.17},"M":{"ATG":1.00},
    "V":{"GTT":0.16,"GTC":0.25,"GTA":0.11,"GTG":0.48},
    "S":{"TCT":0.18,"TCC":0.23,"TCA":0.14,"TCG":0.06,"AGT":0.14,"AGC":0.25},
    "P":{"CCT":0.29,"CCC":0.33,"CCA":0.27,"CCG":0.11},
    "T":{"ACT":0.23,"ACC":0.38,"ACA":0.26,"ACG":0.13},
    "A":{"GCT":0.26,"GCC":0.41,"GCA":0.22,"GCG":0.11},
    "Y":{"TAT":0.40,"TAC":0.60},"*":{"TAA":0.25,"TAG":0.22,"TGA":0.53},
    "H":{"CAT":0.39,"CAC":0.61},"Q":{"CAA":0.23,"CAG":0.77},
    "N":{"AAT":0.42,"AAC":0.58},"K":{"AAA":0.39,"AAG":0.61},
    "D":{"GAT":0.44,"GAC":0.56},"E":{"GAA":0.40,"GAG":0.60},
    "C":{"TGT":0.42,"TGC":0.58},"W":{"TGG":1.00},
    "R":{"CGT":0.09,"CGC":0.20,"CGA":0.11,"CGG":0.22,"AGA":0.19,"AGG":0.19},
    "G":{"GGT":0.16,"GGC":0.36,"GGA":0.25,"GGG":0.23},
}

CODON_TABLES: Final[Mapping[str, CodonTable]] = {
    "e_coli": CODON_TABLE_ECOLI, "human": CODON_TABLE_HUMAN,
    "yeast": CODON_TABLE_YEAST, "pichia": CODON_TABLE_PICHIA, "cho": CODON_TABLE_CHO,
}
ORGANISM_LABELS: Final[Mapping[str, str]] = {
    "e_coli": "Escherichia coli K-12", "human": "Homo sapiens",
    "yeast": "Saccharomyces cerevisiae",
    "pichia": "Pichia pastoris", "cho": "CHO cells",
}

AFFINITY_TAGS_DNA: Final[Mapping[str, str]] = {
    "6xHis":  "CATCATCATCATCATCAT",
    "8xHis":  "CATCATCATCATCATCATCATCAT",
    "FLAG":   "GATTACAAGGATGACGATGACAAG",
    "HA":     "TATCCTTATGACGTTCCAGATTACGCT",
    "Myc":    "GAACAAAAACTCATCTCAGAAGAGGATCTG",
    "Strep-II": "TGGAGCCACCCGCAGTTCGAAAAA",
}

PROTEASE_SITES_DNA: Final[Mapping[str, str]] = {
    "Thrombin": "CTGGTGCCGCGTGGTTCT",
    "TEV": "GAAAACCTGTATTTTCAGGGC",
    "PreScission": "CTGGAAGTGCTGTTCCAGGGCCCA",
    "Factor_Xa": "ATCGAAGGTCGT",
    "Enterokinase": "GACGACGACGACAAG",
    "HRV_3C": "CTGGAAGTTCTGTTCCAGGGGCCC",
}

LINKERS_DNA: Final[Mapping[str, str]] = {
    "GS": "GGTAGC", "GGGGS": "GGTGGCGGTGGTAGC",
    "2xGGGGS": "GGTGGCGGTGGTAGC" * 2, "3xGGGGS": "GGTGGCGGTGGTAGC" * 3,
    "EAAAK": "GAAGCGGCGGCGAAG",
}
