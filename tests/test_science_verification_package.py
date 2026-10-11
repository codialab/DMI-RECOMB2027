import hashlib
from pathlib import Path

import pytest

from analyses.reviewer_checks.science_verification import _package_file


def test_package_file_checks_sha256_and_size(tmp_path: Path):
    payload = tmp_path / "payload.bin"
    payload.write_bytes(b"frozen numerical input")
    manifest = {"files": {"payload.bin": {
        "sha256": hashlib.sha256(payload.read_bytes()).hexdigest(),
        "bytes": payload.stat().st_size,
    }}}
    hashes = {}

    assert _package_file(tmp_path, manifest, "payload.bin", hashes) == payload
    assert hashes["payload.bin"] == manifest["files"]["payload.bin"]["sha256"]

    manifest["files"]["payload.bin"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        _package_file(tmp_path, manifest, "payload.bin", {})


def test_package_file_rejects_paths_outside_bundle(tmp_path: Path):
    manifest = {"files": {"../outside.bin": {"sha256": "0" * 64, "bytes": 1}}}
    with pytest.raises(ValueError, match="escapes data root"):
        _package_file(tmp_path, manifest, "../outside.bin", {})
