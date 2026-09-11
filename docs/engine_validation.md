# Engine validation

Run `pytest gsynth_engine/tests -q --cov=gsynth_engine --cov-config=pyproject.toml`.

CI requires 100% line coverage across every production engine module. Test files are excluded from the coverage denominator. No production module or executable line is excluded. Coverage records are published as CI artifacts.

The suite includes independent Biopython alignment and GenBank checks, generated sequence round trips, malformed trace files, circular feature transport, partial consensus coverage, and explicit failure paths. Line coverage does not measure every possible input or establish experimental biological activity.

SBOL interchange preserves annotation evidence, candidate status, strand, coding bounds, and circular locations. Unsupported discontinuous locations and unsupported GenBank translation codes or exceptions are rejected explicitly.

References: [Biopython alignment](https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html), [SBOL specification](https://sbolstandard.org/datamodel-specification/), [INSDC feature table](https://www.insdc.org/submitting-standards/feature-table/).
