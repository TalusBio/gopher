"""Test that the annotations functions are working correctly"""
from gopher import annotations
import pytest
import re
import pandas as pd


def test_different_species():
    """Check that the annotations can load different species"""
    load = annotations.load_annotations("mgi", "cc")
    assert type(load) is pd.DataFrame


def test_error():
    """Check that annotations loaded function catches errors and returns correct error report"""
    annotations.load_annotations(species="human", aspect="bp")
    annotations.load_annotations(species="yeast", aspect="all")
    with pytest.raises(
        ValueError,
        match=re.escape(
            r"Expected apsect (b) to be one of 'cc', 'mf', 'bp', or 'all'."
        ),
    ):
        annotations.load_annotations(species="human", aspect="b")
