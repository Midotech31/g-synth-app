# Design QA — Design workspace redesign

## Comparison setup

- Source image: `frontend/design-reference/design-1.jpg` (1487 × 1058 px)
- Implementation capture: `frontend/design-qa/implementation-final.jpg` (1425 × 1013 px)
- Side-by-side evidence: `frontend/design-qa/comparison-final.jpg`
- Intended desktop viewport: 1440 × 1024 px
- Tested implementation state: signed-in researcher, default construct generated, verification passed

## Visual and interaction review

| Priority | Area | Finding | Resolution |
| --- | --- | --- | --- |
| P2 | Mobile header | The first implementation let the title shrink beside the action buttons, wrapping the title into a narrow column. | Fixed the header to reserve a full row for the title and progress steps; actions now form a separate equal-width row. |
| P2 | Mobile navigation | The horizontal workspace navigation displayed a browser scrollbar. | Kept the navigation swipeable while hiding the decorative scrollbar. |

No P0 or P1 findings remain. Desktop hierarchy, spacing, navy rail, teal active state, three-step header, verified banner, metric strip, construct map, and hybridisation panel match the selected reference closely. The implementation uses the product's existing icon component and real biological outputs rather than decorative mock data.

## Functional evidence

- Construct generation returned a verified 125 bp construct with 2 oligos.
- Export menu opened and exposed CSV, FASTA, protocol, and GenBank actions.
- Copy action changed to the `Copied` confirmation state.
- The verified result remained visible after navigating Design → Optimise → Design.
- Desktop and phone layouts had no horizontal page overflow.
- Browser console contained no application errors during the tested workflow.

## Validation history

1. Compared the selected reference and first desktop implementation side by side.
2. Exercised the real local Django API and generated a verified construct.
3. Tested navigation persistence and toolbar actions.
4. Found and corrected the mobile header and navigation issues.
5. Re-captured the desktop verified state and repeated the side-by-side review.

## Final Result

passed
