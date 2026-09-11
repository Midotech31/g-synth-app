import pytest
from django.urls import reverse

ENTA_FIXTURE = (
    "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCAACGGCGGC"
    "TTTGTGAACCAGCATCTGTGCGGCAGCCATCTGGTGGAAGCGCTGTACCTGGTGTGCGGCGAACGCGGC"
    "TTTTTTTACACCCCGAAAACCCGCCGCGAAGCGGAAGATCTGCAGGTGGGCCAGGTGGAACTGTAA"
)


def _circular_read(sequence: str, start: int, end: int) -> str:

    body = "".join(sequence[index % len(sequence)] for index in range(start, end))
    return "A" * 30 + body + "T" * 30


@pytest.mark.django_db
def test_entA_fixture_runs_from_clone_design_to_full_sequence_verification(auth_client):
    cloned = auth_client.post(
        reverse("design-clone"),
        {
            "sequence": ENTA_FIXTURE,
            "vector_key": "pET-21a",
            "left_enzyme": "NdeI",
            "right_enzyme": "XhoI",
            "name": "pGS-EntA-release-fixture",
            "save_as_project": True,
        },
    )
    assert cloned.status_code == 200, cloned.data
    construct = cloned.data
    assert construct["is_clonable"]
    assert construct["preflight"]["can_export"]
    assert construct["project_id"]
    assert construct["provenance"]["workflow"] == "cloning"

    primers = auth_client.post(
        reverse("design-primers"),
        {
            "template": construct["plasmid"],
            "target_start": construct["insert_start"],
            "target_end": construct["insert_end"],
            "name": "EntA-release-fixture",
        },
    )
    assert primers.status_code == 200, primers.data
    assert primers.data["covers_target"]
    assert {primer["direction"] for primer in primers.data["primers"]} == {1, -1}

    verified = auth_client.post(
        reverse("design-verify"),
        {
            "design": construct["plasmid"],
            "reads": {
                "T7-F": _circular_read(
                    construct["plasmid"],
                    construct["insert_start"] - 200,
                    construct["insert_end"] + 200,
                ),
            },
            "region_start": construct["insert_start"],
            "region_end": construct["insert_end"],
        },
        format="json",
    )
    assert verified.status_code == 200, verified.data
    assert verified.data["verification_state"] == "fully_verified"
    assert verified.data["coverage"] == 100.0
    assert verified.data["differences"] == []
    assert verified.data["preflight"]["can_export"]
    assert verified.data["provenance"]["workflow"] == "sequence_verification"
