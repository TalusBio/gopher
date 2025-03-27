import pandas as pd
import pytest
from io import StringIO
from pandas.testing import assert_frame_equal

from gopher.parsers.tabular import read_diann


def test_read_diann_removes_metadata_and_sets_index():
    # Simulated DIANN output
    mock_data = StringIO(
        """Protein.Group\tProtein.Ids\tProtein.Names\tGenes\tFirst.Protein.Description\tIntensity.Sample1\tIntensity.Sample2
PG1\tP12345;P67890\tProtein A\tGENE1\tDescription A\t1000\t2000
PG2\tP23456\tProtein B\tGENE2\tDescription B\t1500\t2500
"""
    )

    # Expected DataFrame
    expected = pd.DataFrame(
        {
            # The real diann data has float values in the intensities.
            "Intensity.Sample1": [1000.0, 1500.0],
            "Intensity.Sample2": [2000.0, 2500.0],
        },
        index=["P12345", "P23456"],
    )
    expected.index.name = "Protein"

    result = read_diann(mock_data)
    assert_frame_equal(result, expected)


def test_read_diann_faile_with_gg():
    # Simulated DIANN output
    mock_data = StringIO(
        """Genes\tFirst.Protein.Description\tIntensity.Sample1\tIntensity.Sample2
GENE1\tDescription A\t1000\t2000
GENE2\tDescription B\t1500\t2500
"""
    )

    with pytest.raises(ValueError) as e:
        result = read_diann(mock_data)

    assert "Expected columns" in str(e.value.args[0])
