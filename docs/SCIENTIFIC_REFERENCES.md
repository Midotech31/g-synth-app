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
7. Sharp PM, Li WH. **The codon adaptation index—a measure of directional
   synonymous codon usage bias, and its potential applications.** *Nucleic
   Acids Research.* 1987;15(3):1281–1295.
   <https://doi.org/10.1093/nar/15.3.1281>
8. Athey J, Alexaki A, Osipova E, et al. **A new and updated resource for codon
   usage tables.** *BMC Bioinformatics.* 2017;18:391. HIVE-CUTs separates
   RefSeq and GenBank data, supports versioned tables and documents why older,
   sparsely sampled codon tables can alter optimisation decisions.
   <https://doi.org/10.1186/s12859-017-1793-7>
9. Alexaki A, Kames J, Holcomb DD, et al. **Codon and Codon-Pair Usage Tables
   (CoCoPUTs): facilitating genetic variation analyses and recombinant gene
   design.** *Journal of Molecular Biology.* 2019;431(13):2434–2441.
   <https://doi.org/10.1016/j.jmb.2019.04.021>
10. Ranaghan MJ, Li JJ, Laprise DM, et al. **Assessing optimal: inequalities in
    codon optimization algorithms.** *BMC Biology.* 2021;19:36. The study
    documents strong differences among algorithms and cautions that CAI alone
    is not a reliable predictor of soluble recombinant-protein yield.
    <https://doi.org/10.1186/s12915-021-00968-8>
11. Subramanian K, Payne B, Feyertag F, Alvarez-Ponce D. **The Codon Statistics
    Database: a database of codon usage bias.** *Molecular Biology and
    Evolution.* 2022;39(8):msac157. The database uses reference or
    representative RefSeq genomes and separately reports all nuclear genes and
    a high-expression proxy based on ribosomal-protein genes.
    <https://doi.org/10.1093/molbev/msac157>
12. NCBI. **Reference Sequence (RefSeq) Database.** RefSeq is a curated,
    non-redundant sequence collection with versioned releases.
    <https://www.ncbi.nlm.nih.gov/refseq/>
13. NCBI. **The Genetic Codes.** The reference tables distinguish translation
    assignments from initiator-codon interpretation; an initiator is translated
    as methionine even when an alternative start codon is used.
    <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
14. Ben-Bassat A, Bauer K, Chang SY, et al. **Processing of the initiation
    methionine from proteins: properties of the Escherichia coli methionine
    aminopeptidase and its gene structure.** *Journal of Bacteriology.*
    1987;169(2):751–757.
    <https://doi.org/10.1128/jb.169.2.751-757.1987>
15. Xiao Q, Zhang F, Nacev BA, Liu JO, Pei D. **Protein N-terminal processing:
    substrate specificity of Escherichia coli and human methionine
    aminopeptidases.** *Biochemistry.* 2010;49(26):5588–5599.
    <https://doi.org/10.1021/bi1005464>

## Implementation note

The restriction definitions are a versioned application dataset derived from
REBASE geometry through Biopython plus curated overrides. A reference URL alone
is not a reproducibility record because external databases change; therefore
the full canonical table is hashed into every generated provenance manifest.
The selectable table contains 109 non-redundant in-site cut geometries derived
from 289 commercially available enzyme names. Isoschizomers are retained as
aliases rather than duplicated as map sites. Enzymes with ambiguous recognition
sequences, out-of-site cleavage or two cuts are excluded because those mechanisms
require a different molecular model.

Codon usage is species-, strain- and expression-context dependent. G-Synth
commits all 64 raw codon counts for each of 15 hosts in
`gsynth_engine/data/codon_usage_hive_2021.json`. They were retrieved on
2026-09-02 from FDA HIVE service object 537, labelled September 2021, with
genomic RefSeq species aggregates preferred. *K. phaffii* and *N. benthamiana*
use the disclosed GenBank fallback because that snapshot returned no RefSeq
species record. The update script validates the taxon, 64-codon completeness
and total count before writing the deterministic snapshot; the API exposes its
SHA-256 hash.

Counts are divided by the most frequent synonymous codon for each amino acid to
obtain relative-adaptiveness weights. A score calculated from these species-wide
profiles is labelled **profile-relative CAI**: it is a reproducible design
metric, not a strict CAI against highly expressed genes and not a prediction of
expression yield. For a quantitative CAI claim or a particular strain, tissue
or cell line, G-Synth accepts an explicitly documented highly expressed
reference-gene set and identifies the resulting score as custom-context CAI.

Peptide back-translation separates the translated ORF from the mature product.
In automatic mode, an N-terminal M is treated as an existing initiator; a peptide
without M is preserved as a mature peptide whose translation initiation is
supplied by the subsequent design cassette. Users can override either inference.
If a non-M peptide is explicitly declared a complete ORF, exactly one initiator
Met/ATG is added. G-Synth does not predict post-translational methionine excision,
which depends on the expression system and N-terminal sequence.
