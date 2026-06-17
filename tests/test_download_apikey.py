"""Tests for lsdviztools3.download.apikey — secure key resolution and masking."""

import pytest

from lsdviztools3.download import apikey


def test_env_var_takes_precedence(tmp_path, monkeypatch):
    key_file = tmp_path / "key.txt"
    key_file.write_text("FILEKEY\n")
    monkeypatch.setenv(apikey.ENV_OT_API_KEY, "ENVKEY")
    assert apikey.resolve_ot_api_key(str(key_file)) == "ENVKEY"


def test_reads_first_valid_line_of_key_file(tmp_path, monkeypatch):
    monkeypatch.delenv(apikey.ENV_OT_API_KEY, raising=False)
    key_file = tmp_path / "key.txt"
    key_file.write_text("# a comment\n\nMYKEY\notherstuff\n")
    assert apikey.resolve_ot_api_key(str(key_file)) == "MYKEY"


def test_no_key_raises(monkeypatch):
    monkeypatch.delenv(apikey.ENV_OT_API_KEY, raising=False)
    with pytest.raises(ValueError):
        apikey.resolve_ot_api_key()


def test_missing_key_file_raises(tmp_path, monkeypatch):
    monkeypatch.delenv(apikey.ENV_OT_API_KEY, raising=False)
    with pytest.raises(FileNotFoundError):
        apikey.resolve_ot_api_key(str(tmp_path / "nope.txt"))


def test_mask_key_in_url():
    url = "https://example.org/api?API_Key=SECRET&x=1"
    masked = apikey.mask_key_in_url(url, "SECRET")
    assert "SECRET" not in masked
    assert "***" in masked
