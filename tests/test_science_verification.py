from pathlib import Path

import pytest

from analyses.reviewer_checks import science_verification as verification


def test_hash_mismatch_fails_closed(tmp_path):
    source = tmp_path / "input.tsv"
    source.write_text("frozen\n")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        verification._verify_hash(source, "0" * 64, "synthetic", {})


def test_missing_external_inventory_fails_clearly(tmp_path):
    with pytest.raises(FileNotFoundError, match="required science input is missing"):
        verification._load_inputs(tmp_path, Path(__file__).resolve().parents[1])


def test_verifier_refuses_existing_output_directory(tmp_path):
    out = tmp_path / "existing"
    out.mkdir()
    with pytest.raises(FileExistsError, match="choose a new --output-root"):
        verification.run_verification(level="2a", data_root=tmp_path,
            public_root=tmp_path, workers=1, output_root=out)


def test_reproduce_cli_routes_explicit_science_arguments(monkeypatch, capsys, tmp_path):
    import reproduce
    received = {}

    def fake_run(**kwargs):
        received.update(kwargs)
        return {"status": "PASS", "requested_level": kwargs["level"]}

    monkeypatch.setattr(verification, "run_verification", fake_run)
    monkeypatch.setattr("sys.argv", ["reproduce.py", "verify-science",
        "--level", "2b", "--workers", "4", "--batch-size", "3",
        "--data-root", str(tmp_path / "private"),
        "--output-root", str(tmp_path / "report")])
    assert reproduce.main() == 0
    assert received["level"] == "2b"
    assert received["workers"] == 4
    assert received["batch_size"] == 3
    assert received["data_root"] == tmp_path / "private"
    assert "Scientific verification: PASS" in capsys.readouterr().out


def test_prediction_parity_ignores_only_expected_structural_ties():
    import numpy as np

    assert verification._prediction_max_difference(
        "correct", np.array([1.0, 7.0, 2.0]),
        np.array([1.0 + 1e-13, np.nan, 2.0]),
        np.array([True, False, True]), "pair") == pytest.approx(1e-13)


def test_prediction_parity_rejects_nan_mask_mismatch():
    import numpy as np

    with pytest.raises(ValueError, match="finite-population mismatch"):
        verification._prediction_max_difference(
            "correct", np.array([1.0, np.nan]), np.array([1.0, 2.0]),
            np.array([True, False]), "pair")
