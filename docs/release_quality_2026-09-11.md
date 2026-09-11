# Release quality review — 11 September 2026

Baseline: `e1a95762991bb141879036944e81b73b3efcb202` (PR #33).

## Target and acceptance criteria

The requested 9/10 is a quality target, not a scientific certification. A release
must pass the existing engine coverage gate (95%), API, interface, build,
dependency and browser workflow checks. Scores must not be substituted for
those results. Independent user testing and experimental validation are distinct
from software tests and are not claimed here.

## Confirmed problems and corrections

| Problem | Correction | Evidence |
|---|---|---|
| Sequence PATCH could leave the original output hash and design evidence attached to different DNA | Saved project sequences are immutable; changed sequences must be saved as a new project | API regression reproduces the old 200 response and now requires 400; stored sequence and provenance remain unchanged |
| General project writes bypassed annotation coordinate validation and accepted non-object JSON | Validate annotation arrays and their coordinates on general project writes as well as annotation edits | Invalid coordinates, arrays and strings rejected |
| Two tabs could overwrite each other's annotations | Viewer sends the loaded revision; API uses an atomic compare-and-swap and returns 409 for stale revisions | Two writes using one revision preserve the first edit |
| Late saves or scans could replace another project's visible state | Route lifecycle and request-version guards; editor and scan state reset on navigation | Deferred-response frontend tests |
| New SD context text exceeded the 300-character save limit | Evidence notes accept 4,000 characters through design, saving and GenBank import; lists of associated CDS names are bounded | Real detector output saves intact |
| GenBank CDS import ignored codon_start, and export ignored actual translation bounds | Import honors offsets on either strand; export uses the actual coding interval and codon_start=1 | Twelve independent Biopython extraction cases across both strands, all three offsets and circular-origin crossings; four engine export cases |
| Imported coding bounds were not transported when a vector was reversed or rotated | Mirror and move translation bounds with features; retain reverse codon phase on truncation; preserve unstranded direction | Involution and backbone sequence-extraction checks |
| Discontinuous GenBank features were flattened into misleading continuous spans | Reject unsupported discontinuous locations with a specific explanation | Split-feature import regression |
| Feature review preselected every detected candidate on a manual scan | Require an explicit choice; keep candidate status after review | Frontend selection regression |
| Annotation editing lacked an evidence field and save errors appeared behind the dialog | Add a source/evidence note and show save errors in the editor | Browser screenshot and editor workflow |
| Saved-viewer feature lists were hard to search and wrapped positions exceeded plasmid length | Search by feature name/type; show physical coordinates and “across origin” | UI build and existing wrapped-coordinate tests |

## Scientific basis

- [INSDC Feature Table Definition, version 11.4](https://www.insdc.org/submitting-standards/feature-table/): feature locations, complement/join and coding-frame qualifiers.
- [INSDC inference vocabulary](https://www.insdc.org/submitting-standards/inference-qualifiers/): computational identification is non-experimental evidence.
- The SD references and heuristic limits remain documented in
  [the sequence explorer review](sequence_explorer_review.md).

The independent checks use Biopython to extract the annotated DNA from exported
records and compare it with the intended strand-oriented interval. These are
software interoperability checks, not measured expression or wet-lab outcomes.
No reference vector bases were changed or newly copied from a supplier.

## Boundaries and remaining acceptance work

- Direct authenticated production UX review was blocked by the sign-in screen
  in the cloud browser; no account access was weakened. CI runs the isolated
  full browser workflow and captures screenshots for this change.
- This is a targeted release review, not a claim that every line of the
  application has been independently audited.
- Discontinuous features require a future multipart-coordinate model. They
  must not be flattened merely to make an import succeed.
- GenBank is the evidence-preserving export reviewed here. SBOL metadata parity,
  alternative genetic-code display and complex feature topology need a separate
  interoperability review before universal format support is claimed.
- Existing project records are not silently rewritten. Reimporting the original
  file is necessary to recover frame metadata that older imports discarded.
- Projects created by earlier software may lack a provenance hash; this change
  does not invent historical evidence for them.
- A defensible overall 9/10 still requires independent domain review and observed
  task completion by new users. No such human or experimental results are fabricated.

Final CI and deployment identifiers are linked in the pull request that delivers
this document; they identify the exact tested release.
