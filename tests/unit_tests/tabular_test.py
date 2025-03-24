import pandas as pd
from io import StringIO
from pandas.testing import assert_frame_equal

from gopher.parsers.tabular import read_diann
import os

def test_read_diann_removes_metadata_and_sets_index():
    # Simulated DIANN output
    mock_data = StringIO(
        """Protein.Group\tProtein.Ids\tProtein.Names\tGenes\tFirst.Protein.Description\tIntensity.Sample1\tIntensity.Sample2
PG1\tP12345;P67890\tProtein A\tGENE1\tDescription A\t1000\t2000
PG2\tP23456\tProtein B\tGENE2\tDescription B\t1500\t2500
"""
    )

    # Expected DataFrame
    expected = pd.DataFrame({
        "Intensity.Sample1": [1000, 1500],
        "Intensity.Sample2": [2000, 2500],
    }, index=["P12345", "P23456"])
    expected.index.name = "Protein"

    result = read_diann(mock_data)

    assert_frame_equal(result, expected)