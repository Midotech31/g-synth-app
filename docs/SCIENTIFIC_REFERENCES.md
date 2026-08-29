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

## Implementation note

The restriction definitions are a versioned application dataset derived from
REBASE geometry through Biopython plus curated overrides. A reference URL alone
is not a reproducibility record because external databases change; therefore
the full canonical table is hashed into every generated provenance manifest.
