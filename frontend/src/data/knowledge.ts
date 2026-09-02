export type KnowledgeCategory = "Design" | "PCR" | "Cloning" | "Validation" | "Expression";

export type KnowledgeTopic = {
  id: string;
  title: string;
  category: KnowledgeCategory;
  summary: string;
  essentials: string[];
  caution: string;
  route: string;
  routeLabel: string;
  keywords: string[];
  references: { label: string; url: string }[];
};

const REBASE = {
  label: "Roberts et al., REBASE, Nucleic Acids Research (2023)",
  url: "https://doi.org/10.1093/nar/gkac975",
};
const NEB_CLONING = {
  label: "NEB, Traditional Cloning Quick Guide",
  url: "https://www.neb.com/en/tools-and-resources/usage-guidelines/cloning-guide",
};
const NEB_ENDS = {
  label: "NEB, Cleavage Close to the End of DNA Fragments",
  url: "https://www.neb.com/en-us/tools-and-resources/usage-guidelines/cleavage-close-to-the-end-of-dna-fragments",
};
const SANTALUCIA = {
  label: "SantaLucia, nearest-neighbour DNA thermodynamics, PNAS (1998)",
  url: "https://doi.org/10.1073/pnas.95.4.1460",
};

export const KNOWLEDGE_TOPICS: KnowledgeTopic[] = [
  {
    id: "restriction-cloning",
    title: "Restriction-enzyme cloning",
    category: "Cloning",
    summary: "Choose cutters by their actual sites, end compatibility and position in the vector and insert.",
    essentials: [
      "For directional cloning, two different compatible ends usually force one insert orientation and reduce vector self-ligation.",
      "Count every recognition site in the circular vector and the insert. A supposed single cutter that cuts twice changes the expected backbone and digest products.",
      "A site can cross the arbitrary origin of a circular sequence; it must still be counted and displayed.",
      "Recognition sequence alone is not the whole protocol: methylation sensitivity, buffer, temperature, heat inactivation and star activity remain enzyme- and supplier-specific.",
    ],
    caution: "G-Synth predicts sequence geometry. Confirm current reaction conditions and methylation constraints in the supplier datasheet for the exact enzyme lot and substrate.",
    route: "/clone", routeLabel: "Open Clone",
    keywords: ["restriction", "enzyme", "hindiii", "hindi", "ndei", "xhoi", "sticky", "overhang", "ligation", "directional", "digest"],
    references: [REBASE, NEB_CLONING],
  },
  {
    id: "tailed-primers",
    title: "Cloning primers and 5′ tails",
    category: "PCR",
    summary: "A cloning primer contains an annealing region plus a non-hybridising 5′ addition.",
    essentials: [
      "In cycle 1, only the primer's 3′ template-complementary region anneals. The 5′ clamp and restriction site remain single-stranded.",
      "Polymerase extends from the correctly paired 3′ end and copies the full primer into the product; later cycles can therefore contain the added tail on both strands.",
      "Use the annealing region—not the full tailed oligo—to estimate the initial annealing temperature.",
      "Extra 5′ bases upstream of a restriction site can improve end cleavage, but the required number and sequence are enzyme-specific.",
    ],
    caution: "A predicted Tm is a starting point. Polymerase, salt, Mg²⁺, additives, template complexity and primer concentration can shift the working annealing temperature.",
    route: "/pcr", routeLabel: "Open PCR",
    keywords: ["primer", "tail", "clamp", "anneal", "hybrid", "hybridization", "hybridisation", "tm", "cycle", "pcr"],
    references: [SANTALUCIA, NEB_ENDS],
  },
  {
    id: "pcr-controls",
    title: "PCR design and controls",
    category: "PCR",
    summary: "Primer specificity, matched annealing behaviour and controls matter more than a single ideal number.",
    essentials: [
      "The two 3′ ends define the amplified interval. Check that each annealing sequence occurs at the intended site and orientation.",
      "Avoid strong 3′ complementarity within or between primers because extension from a primer dimer can compete with the intended product.",
      "Include a no-template control; for diagnostic assays, include a known positive control and interpret it separately from the experimental sample.",
      "An expected amplicon size is a prediction. Identity still requires sequence- or locus-specific confirmation.",
    ],
    caution: "G-Synth does not model every secondary structure, reagent formulation or off-target in a complete genome. Review the real template context before ordering.",
    route: "/pcr", routeLabel: "Design PCR",
    keywords: ["pcr", "primer", "amplicon", "ntc", "control", "dimer", "specificity", "temperature", "annealing"],
    references: [SANTALUCIA],
  },
  {
    id: "diagnostic-gels",
    title: "Diagnostic digests and agarose gels",
    category: "Validation",
    summary: "Digest fragments must follow from actual cut coordinates and add up to the molecule length.",
    essentials: [
      "For a complete digest of a circular plasmid, n distinct cut positions produce n fragments; their sizes sum to the plasmid length.",
      "For a linear DNA molecule, n cuts produce n + 1 fragments.",
      "Choose a ladder that brackets the expected bands. Closely sized fragments may co-migrate and a very small fragment may be faint or leave the gel.",
      "Uncut plasmids can show several conformations and should not be interpreted as though migration depended only on base-pair length.",
    ],
    caution: "The G-Synth gel is an in-silico migration sketch, not an experimental image. Band intensity, partial digestion, conformation and detection limits require bench data.",
    route: "/clone", routeLabel: "Simulate a diagnostic digest",
    keywords: ["gel", "ladder", "marker", "digest", "band", "agarose", "fragment", "migration", "size"],
    references: [NEB_CLONING],
  },
  {
    id: "codon-optimisation",
    title: "Codon optimisation without changing protein",
    category: "Design",
    summary: "Synonymous codons can be changed for a host while the translated amino-acid sequence remains fixed.",
    essentials: [
      "Back-translate each amino acid with synonymous codons from the selected host table and verify the translated product afterwards.",
      "Bundled species profiles are reproducible defaults. A strain-, tissue- or cell-line-specific claim should use a documented reference-gene set from that expression context.",
      "Optimisation is multi-constraint: codon usage, GC content, repeats, homopolymers, unwanted restriction sites and difficult motifs can conflict.",
      "Preserving protein sequence is a hard gate. A favourable codon score is not allowed to hide an amino-acid change.",
      "Codon usage is only one determinant of expression; mRNA structure, transcription, toxicity, folding and culture conditions also matter.",
    ],
    caution: "A computationally optimised coding sequence does not guarantee soluble or functional protein expression.",
    route: "/optimise", routeLabel: "Open Optimise",
    keywords: ["codon", "cai", "optimise", "optimize", "expression", "host", "synonymous", "gc", "protein"],
    references: [{ label: "Kazusa codon-usage database", url: "https://www.kazusa.or.jp/codon/" }],
  },
  {
    id: "construct-annotation",
    title: "Construct annotation and ORF integrity",
    category: "Design",
    summary: "Annotations connect named biological features to exact coordinates, strand and sequence context.",
    essentials: [
      "Record promoter, operator, RBS, CDS, tags, cleavage sites, terminator, origins and resistance markers as separate coordinate-level features.",
      "For a CDS, verify the intended start, reading frame, stop policy and translated sequence across cloning junctions.",
      "A feature may cross the origin of a circular record; this is one biological span represented across the coordinate boundary.",
      "Exact motif detection is a proposal for review. A sequence match alone does not prove promoter activity, protein binding or cleavage in the biological context.",
    ],
    caution: "Imported annotations can be incomplete or wrong. Review feature names, coordinates, strand and translation before using them as experimental evidence.",
    route: "/projects", routeLabel: "Open Projects",
    keywords: ["annotation", "feature", "orf", "cds", "promoter", "rbs", "terminator", "tag", "origin", "insert"],
    references: [],
  },
  {
    id: "sanger-validation",
    title: "Post-sequencing validation",
    category: "Validation",
    summary: "Align quality-bearing reads to the designed reference and distinguish coverage from identity.",
    essentials: [
      "Inspect the ABIF or SCF chromatogram as well as base calls. Low-quality ends and overlapping peaks can create unsupported apparent variants.",
      "Place forward and reverse reads in their correct orientation, trim poor-quality sequence and review every mismatch, insertion and deletion against trace evidence.",
      "Report reference coverage separately from identity within covered bases. An uncovered interval is unknown, not confirmed.",
      "Confirm both cloning junctions, the complete insert and any sequence-critical regulatory or fusion region with adequate bidirectional evidence where possible.",
    ],
    caution: "G-Synth supports trace review and alignment; final acceptance still depends on the laboratory's validated sequencing and quality criteria.",
    route: "/verify", routeLabel: "Open Check",
    keywords: ["sanger", "sequencing", "ab1", "abif", "scf", "chromatogram", "trace", "alignment", "coverage", "mismatch", "variant"],
    references: [],
  },
  {
    id: "expression-cassette",
    title: "Expression-cassette logic",
    category: "Expression",
    summary: "Promoter, operator, RBS, start, ORF, tags, stop and terminator must form one coherent architecture.",
    essentials: [
      "Check feature order and strand first; a correctly spelled motif in the wrong orientation or position does not create a coherent cassette.",
      "A vector-encoded C-terminal tag is translated only when the insert remains in frame through the junction and does not stop beforehand.",
      "Protease sites are peptide motifs: synonymous DNA can encode the same cleavage sequence.",
      "Confirm that cloning has not duplicated a start codon, shifted the frame or inserted unintended residues at either junction.",
    ],
    caution: "Sequence coherence supports an expression hypothesis; it does not establish transcription level, protein folding, processing or biological activity.",
    route: "/projects", routeLabel: "Inspect an annotated construct",
    keywords: ["expression", "cassette", "pet21", "t7", "lac", "rbs", "his", "tag", "thrombin", "stop", "frame"],
    references: [],
  },
  {
    id: "oligo-assembly",
    title: "Oligo assembly and synthesis planning",
    category: "Design",
    summary: "A long design becomes buildable only when overlaps, strand pairing and assembly order remain explicit.",
    essentials: [
      "Every adjacent fragment needs the intended complementary overlap or compatible sticky end; verify both strands rather than inferring one from the other.",
      "Oligo count, maximum length, overlap length and terminal overhangs must be reported with the final sequence so the orderable material is auditable.",
      "Assembly plans should preserve the exact designed sequence after overlap collapse and flag any ambiguity before ordering.",
      "Synthetic oligos can contain synthesis errors; clone screening and sequence confirmation remain necessary.",
    ],
    caution: "In-silico assembly proves sequence consistency under the stated model, not physical assembly efficiency or clone correctness.",
    route: "/design", routeLabel: "Open Design",
    keywords: ["oligo", "assembly", "overlap", "ssd", "esd", "synthesis", "sticky", "ligation", "fragment"],
    references: [],
  },
  {
    id: "evidence-boundaries",
    title: "What in-silico validation can and cannot prove",
    category: "Validation",
    summary: "A sequence-derived pass is evidence about the design, not a substitute for a physical measurement.",
    essentials: [
      "Sequence checks can establish predicted length, coordinates, translation, restriction geometry, primer placement and expected digest fragments.",
      "Simulations cannot establish DNA yield, ligation efficiency, transformation success, expression, solubility, purity, potency or biological activity.",
      "Record inputs, software version, parameters and checksums so another scientist can reproduce the computational result.",
      "Keep predicted, observed and externally confirmed results visibly distinct in figures, tables and manuscripts.",
    ],
    caution: "Use experimental controls, traceable raw data and predefined acceptance criteria for every claim about a physical construct or product.",
    route: "/help", routeLabel: "Read G-Synth Help",
    keywords: ["evidence", "validation", "simulation", "prediction", "experimental", "proof", "reproducibility", "limitation"],
    references: [],
  },
];

export const KNOWLEDGE_CATEGORIES: ("All" | KnowledgeCategory)[] = [
  "All", "Design", "PCR", "Cloning", "Validation", "Expression",
];

const SEARCH_STOP_WORDS = new Set([
  "a", "an", "and", "are", "can", "do", "does", "for", "how", "in", "is",
  "it", "not", "of", "on", "or", "the", "to", "what", "when", "why", "with",
]);

export function searchKnowledge(query: string, category: "All" | KnowledgeCategory = "All") {
  const terms = query.toLowerCase()
    .replace(/[’′'"?!.,:;()\[\]{}]/g, " ")
    .trim()
    .split(/\s+/)
    .filter((term) => term.length > 1 && !SEARCH_STOP_WORDS.has(term));
  const ranked = KNOWLEDGE_TOPICS
    .filter((topic) => category === "All" || topic.category === category)
    .map((topic) => {
      const title = topic.title.toLowerCase();
      const haystack = [topic.title, topic.summary, ...topic.essentials, topic.caution, ...topic.keywords]
        .join(" ").toLowerCase();
      const score = terms.reduce((total, term) => total
        + (title.includes(term) ? 4 : 0)
        + (topic.keywords.some((keyword) => keyword.includes(term) || term.includes(keyword)) ? 3 : 0)
        + (haystack.includes(term) ? 1 : 0), 0);
      return { topic, score };
    })
    .sort((a, b) => b.score - a.score || a.topic.title.localeCompare(b.topic.title));
  if (!terms.length) return ranked.map(({ topic }) => topic);
  const strongest = ranked[0]?.score ?? 0;
  const threshold = Math.max(2, Math.ceil(strongest * 0.35));
  return ranked.filter(({ score }) => score >= threshold).map(({ topic }) => topic);
}
