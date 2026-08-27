# Bench-scientist usability study

**Status:** protocol and capture sheets ready; independent participant sessions
have not yet been run. Do not report SUS, completion, or error-rate results until
signed participant records exist.

This protocol is the remaining human-validation gate. Automated tests can
prove consistency; they cannot prove that a scientist interprets a result
correctly under time pressure.

## Study question

Can a bench scientist who did not build G-Synth complete the workflow
Insert → Strategy → Preflight → Order → Clone → Verify without assistance,
and correctly distinguish a ready design, a review item, and a blocked design?

## Participants

Recruit 5–8 participants spanning one novice, several routine molecular-biology
users, and at least one senior reviewer. Do not use project contributors for
the primary result.

## Tasks

1. Design and export an NdeI/XhoI coding insert in Guided mode.
2. Switch to Expert mode and explain which choices changed.
3. Diagnose a construct with an internal restriction site.
4. Clone into a bundled vector and identify the expected backbone/insert bands.
5. Download the bench worksheet and locate its provenance checksum.
6. Classify five sequencing outcomes: fully verified, differences detected,
   partial match, reads unplaced, and not checked.
7. Recover an unfinished draft after refreshing the browser.

## Measures and pass criteria

- ≥90% task completion without facilitator intervention.
- 100% correct interpretation of blocked versus ready.
- 100% refusal to call partial coverage “verified.”
- Median System Usability Scale (SUS) ≥80.
- No critical error in enzyme choice, primer ordering, exported molecule, or
  sequencing verdict.
- Record time on task, wrong turns, help-page visits, and verbatim uncertainty.

## Session script

Ask participants to think aloud. Do not teach the interface during a task.
After each task ask: “What would you do next at the bench, and what on this
screen supports that decision?” Capture the screen and audio only with written
consent; use synthetic sequences and no personal or patient data.

## Issue triage

Classify findings as critical (could create/release the wrong molecule), major
(blocks task completion), moderate (recoverable confusion), or minor. A
critical finding blocks release until fixed and retested with a new participant.

## Per-participant capture sheet

Use an anonymous ID; keep consent records separately.

| Field | Value |
| --- | --- |
| Participant ID | |
| Experience band | novice / routine / senior reviewer |
| Independent of project | yes / no |
| Browser, OS, device | |
| Assistive technology (if any) | |
| Consent recorded | yes / no |

For every task record completion (independent / prompted / failed), elapsed
time, wrong turns, help visits, interpretation, next bench action, and a short
verbatim uncertainty quote. Any wrong ready/review/block decision is critical.

## SUS questionnaire and scoring

After all tasks, collect the standard ten SUS responses on a 1–5 scale. For odd
items subtract 1 from the response; for even items subtract the response from
5. Sum the ten contributions and multiply by 2.5. Report the median across
participants, the range, and the number of complete questionnaires. Do not
impute missing responses.

## Study result record

| Metric | Release threshold | Observed |
| --- | ---: | ---: |
| Independent task completion | ≥90% | pending |
| Correct ready vs blocked interpretation | 100% | pending |
| Partial coverage rejected as verified | 100% | pending |
| Median SUS | ≥80 | pending |
| Critical scientific errors | 0 | pending |

The study is complete only after 5–8 independent bench scientists have signed
consent, all raw task sheets are retained, issues are triaged, critical fixes
are retested with a new participant, and a reviewer signs the result record.
