"""Test that the annotations functions are working correctly."""

import re

import pandas as pd
import pytest

from gopher import annotations


def test_different_species():
    """Check that the annotations can load different species."""
    load, _ = annotations.load_annotations("mgi", "cc")
    assert isinstance(load, pd.DataFrame)


def test_error():
    """Check that loading annotations catches errors and reports correctly."""
    annotations.load_annotations(species="human", aspect="bp")
    annotations.load_annotations(species="yeast", aspect="all")
    with pytest.raises(
        ValueError,
        match=re.escape(
            r"Expected apsect (b) to be one of 'cc', 'mf', 'bp', or 'all'."
        ),
    ):
        annotations.load_annotations(species="human", aspect="b")


def test_generate_annotations():
    """Test generate_annotations returns the expected dataframe."""
    prot = ["P10809", "P35527", "Q9UMS4", "P52907", "Q9NV31"]
    aspect = "mf"
    go_name = "temporary"
    data = {"uniprot_accession": prot, "aspect": aspect, "go_name": go_name}
    data_manual = pd.DataFrame.from_dict(data)
    result = annotations.generate_annotations(prot, aspect, go_name)
    assert len(result) == len(prot)
    assert data_manual.equals(result.drop(columns=["go_id"]))
    full, _ = annotations.load_annotations(species="human")
    final = pd.concat([result, full])
    assert len(final) == len(full) + len(result)
