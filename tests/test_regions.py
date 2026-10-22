"""Tests for reading boxes pasted from Google Maps."""

from __future__ import annotations

import pytest

from lsdviztools3.regions import Box, box_from_corners, parse_lat_lon


def test_parse_lat_lon_reads_google_maps_text():
    assert parse_lat_lon("55.9476273103734, -3.158220317286066") == (55.9476273103734, -3.158220317286066)


def test_parse_lat_lon_ignores_brackets_and_spaces():
    assert parse_lat_lon("  (-33.8568, 151.2153) ") == (-33.8568, 151.2153)


@pytest.mark.parametrize("text", ["55.94", "", "abc", "1, 2, 3", "55°56'51\"N 3°09'29\"W"])
def test_parse_lat_lon_rejects_text_without_two_numbers(text):
    with pytest.raises(ValueError, match="Could not find"):
        parse_lat_lon(text)


def test_parse_lat_lon_rejects_out_of_range():
    with pytest.raises(ValueError, match="latitude"):
        parse_lat_lon("155.9, -3.15")
    with pytest.raises(ValueError, match="longitude"):
        parse_lat_lon("55.9, 190.0")


def test_box_from_corners_any_opposite_pair_gives_same_box(capsys):
    a = box_from_corners("55.90, -3.25", "55.99, -3.10")
    b = box_from_corners("55.99, -3.25", "55.90, -3.10", verbose=False)
    assert a == b == Box(west=-3.25, east=-3.10, south=55.90, north=55.99)
    west, east, south, north = a
    assert (west, east, south, north) == (-3.25, -3.10, 55.90, 55.99)
    assert a.bbox == (-3.25, 55.90, -3.10, 55.99)
    printed = capsys.readouterr().out
    assert "west -3.25000" in printed and "km east-west" in printed


def test_box_from_corners_rejects_zero_area():
    with pytest.raises(ValueError, match="different latitudes"):
        box_from_corners("55.9, -3.2", "55.9, -3.1")
