# Accessibility and mobile-layout audit — 2026-08-27

## Audit scope

The signed-in entry, Design, Clone, and Check workflows were exercised in the
real React application against a migrated local API. Desktop was checked at
1440 × 900 and mobile at 390 × 844. The review covered responsive reflow,
page-level overflow, visible labels and landmarks, live status text, text
contrast, control target size, reduced-motion behavior, and the ability to
distinguish empty, ready, review, and blocked states.

Target: WCAG 2.2 AA for the application-owned interface. This is an engineering
audit, not a conformance certification.

## Flow results

1. **Sign in — healthy.** Email and password have programmatic labels and
   appropriate autocomplete purposes; failed authentication is an alert.
2. **Dashboard — healthy.** A main landmark, named workspace navigation, clear
   heading order, and a first-focus skip link are present.
3. **Design input and preflight — healthy after fixes.** The form reflows to one
   column with no page-level horizontal overflow. Preflight uses text, icons,
   and stable codes rather than color alone.
4. **Clone setup — healthy after fixes.** The bundled-vector choice and insert
   controls remain readable at 390 px. Long navigation and scientific sequence
   material use bounded, intentional horizontal scrolling rather than widening
   the page.
5. **Sequencing check — healthy after fixes.** Empty state, workflow tabs, and
   disabled release action remain legible and correctly named at mobile width.

## Issues found and resolved

- Secondary text used `#78889b` on white, a measured contrast of 3.62:1.
  `--muted` and `--ink-faint` now use `#607185`, measured at 5.00:1.
- Repeated buttons, navigation links, selectors, and mode switches measured
  28–42 px high. Application-owned interactive targets now provide at least
  44 px in their clickable dimension; checkbox labels provide the 44 px target
  around the 20 px control.
- `--ink-faint` and `--paper-deep` were referenced but undefined. Both tokens
  are now explicit, preventing browser-dependent inherited colors/backgrounds.
- Reduced-motion mode previously kept the spinner moving, only more slowly.
  It now renders a static busy indicator while the adjacent live text continues
  to announce progress.

Post-fix checks on Design and Check reported no page-level horizontal overflow,
no sub-44 px application-owned button/link/text-input targets, and no sampled
small-text contrast failure.

## Confirmed strengths

- Skip link and focus styles are present.
- `main`, complementary navigation, and named navigation landmarks are present.
- Form controls have visible/programmatic labels.
- Errors use alerts; loading and successful design messages use polite live
  status regions.
- Toggle groups expose pressed state and accessible group names.
- Mobile sequence/table overflow is contained locally.
- Status meaning is not conveyed by color alone.

## Evidence limits

This pass did not use NVDA, JAWS, VoiceOver, TalkBack, Windows High Contrast,
or 200–400% browser zoom on physical devices. Automated DOM and contrast checks
cannot establish full WCAG conformance. Those assistive-technology/device checks
remain a release acceptance activity and must be recorded with browser, OS,
screen reader, version, and observed outcome.
