"""File I/O — FASTA, GenBank, SnapGene, NCBI, UniProt."""
from gsynth_core.io.records import Feature, FeatureLocation, SeqRecord
from gsynth_core.io.fasta import parse_fasta, write_fasta
from gsynth_core.io.genbank import parse_genbank, write_genbank
from gsynth_core.io.snapgene import parse_snapgene
from gsynth_core.io.ncbi import fetch_ncbi, fetch_uniprot
__all__ = [
    "Feature", "FeatureLocation", "SeqRecord",
    "parse_fasta", "write_fasta",
    "parse_genbank", "write_genbank",
    "parse_snapgene",
    "fetch_ncbi", "fetch_uniprot",
]
