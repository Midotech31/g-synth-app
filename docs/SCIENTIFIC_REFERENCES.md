# Scientific references used by G-Synth

The software records its exact enzyme-table checksum in provenance. These
references explain the biological assumptions; current lot-specific reaction
conditions must still be checked against the relevant supplier datasheet.

1. Roberts RJ, Vincze T, Posfai J, Macelis D. **REBASE: a database for DNA
   restriction and modification: enzymes, genes and genomes.** *Nucleic Acids
   Research.* 2023;51(D1):D629–D630.
   <https://doi.org/10.1093/nar/gkac975>
2. New England Biolabs. **Cloning With Restriction Enzymes.** Selection guidance
   covers compatible ends, internal recognition sites, and methylation
   sensitivity.
   <https://www.neb.com/en-us/tools-and-resources/video-library/cloning-with-restriction-enzymes>
3. New England Biolabs. **Cleavage Close to the End of DNA Fragments.** The
   general guidance for enzymes without specific data is six terminal base
   pairs, while noting enzyme-specific results and asymmetric cleavage.
   <https://www.neb.com/en-ca/tools-and-resources/usage-guidelines/cleavage-close-to-the-end-of-dna-fragments>
4. New England Biolabs. **Traditional Cloning Quick Guide.** Directional
   cloning, PCR-added sites, insert purification, ligation, and transformation.
   <https://www.neb.com/en/tools-and-resources/usage-guidelines/cloning-guide>
5. SantaLucia J Jr. **A unified view of polymer, dumbbell, and oligonucleotide
   DNA nearest-neighbor thermodynamics.** *PNAS.* 1998;95(4):1460–1465.
   <https://doi.org/10.1073/pnas.95.4.1460>
6. New England Biolabs. **Factor Xa Protease.** The preferred cleavage site is
   Ile-Glu/Asp-Gly-Arg↓, establishing that recognition is peptide-based rather
   than tied to one synonymous DNA sequence.
   <https://www.neb.com/en-us/products/p8010-factor-xa-protease>
7. Kazusa DNA Research Institute. **Codon usage table: Escherichia coli
   W3110.** The table records GAA, GGT and CGT as substantially more frequent
   than less-preferred synonymous choices, especially AGG.
   <https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?aa=1&species=316407&style=GCG>
8. Nakamura Y, Gojobori T, Ikemura T. **Codon usage tabulated from the
   international DNA sequence databases: status for the year 2000.** *Nucleic
   Acids Research.* 2000;28(1):292.
   <https://doi.org/10.1093/nar/28.1.292>
9. Kazusa DNA Research Institute. **Species-wide codon usage tables.** The
   bundled host profiles use the archived nuclear coding-sequence frequencies
   for *Saccharomyces cerevisiae* (taxon 4932), *Pichia pastoris* / *Komagataella
   phaffii* (4922), *Bacillus subtilis* (1423), *Homo sapiens* (9606),
   *Cricetulus griseus* (10029), *Spodoptera frugiperda* (7108) and *Nicotiana
   benthamiana* (4100). Frequencies are normalized within each synonymous-codon
   family before use.
   <https://www.kazusa.or.jp/codon/>
10. Alexaki A, Kames J, Holcomb DD, et al. **Codon and Codon-Pair Usage Tables
    (CoCoPUTs): facilitating genetic variation analyses and recombinant gene
    design.** *Journal of Molecular Biology.* 2019;431(13):2434–2441.
    <https://doi.org/10.1016/j.jmb.2019.04.021>

## Implementation note

The restriction definitions are a versioned application dataset derived from
REBASE geometry through Biopython plus curated overrides. A reference URL alone
is not a reproducibility record because external databases change; therefore
the full canonical table is hashed into every generated provenance manifest.

Codon usage is species-, strain- and expression-context dependent. The bundled
profiles are reproducible defaults, not predictions of expression yield. For a
quantitative CAI claim or a specific strain, tissue or cell line, G-Synth should
be given an explicitly documented reference-gene set from that context.
