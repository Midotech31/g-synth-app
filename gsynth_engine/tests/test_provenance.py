from gsynth_engine.provenance import build_provenance


def test_provenance_is_reproducible_and_does_not_duplicate_sequences():
    parameters = {
        "sequence": "ATGAAATAA",
        "left_enzyme": "NdeI",
        "right_enzyme": "XhoI",
    }
    first = build_provenance("ssd", parameters=parameters, output_sequence="TATGAAATAA")
    second = build_provenance("ssd", parameters=parameters, output_sequence="TATGAAATAA")
    assert first == second
    assert first["parameters"]["sequence"]["length"] == 9
    assert "ATGAAATAA" not in str(first)
    assert len(first["output_sha256"]) == 64
    assert len(first["enzyme_table"]["sha256"]) == 64


def test_provenance_changes_when_a_design_choice_changes():
    base = build_provenance(
        "ssd", parameters={"left_enzyme": "NdeI"}, output_sequence="ATG"
    )
    changed = build_provenance(
        "ssd", parameters={"left_enzyme": "NcoI"}, output_sequence="ATG"
    )
    assert base["parameters_sha256"] != changed["parameters_sha256"]
