import pytest
from pathlib import Path
import pandas as pd
from gopher import normalize
import numpy as np


@pytest.fixture
def real_data(tmp_path):
    """Test using small files."""
    fasta_df = Path("../data/small-yeast.fasta")
    quant = pd.read_csv("../data/yeast_small.csv")
    quant = quant.set_index("Protein")

    return quant, fasta_df


def test_normalize(real_data):
    """Check that the normalization returns a dataframe."""
    quant, fasta = real_data
    df = normalize.normalize_values(quant, fasta)
    assert isinstance(df, pd.DataFrame)


def test_normalize_one_prot(real_data):
    """Check that the normalization returns the correct result."""
    quant, fasta = real_data
    # Create a subset of data that is only 1 protein
    single_prot_quant = quant.iloc[:1, :]
    # Get the fasta
    fasta_df = normalize.read_fasta(fasta)
    # Perform the manual calculation
    vals = single_prot_quant.values
    mass = fasta_df[fasta_df["Protein"] == single_prot_quant.index[0]][
        "Mass"
    ].iloc[0]
    manual_result = vals / vals / mass
    # Get the calculation from the function and compare the results
    result = normalize.normalize_values(single_prot_quant, fasta).values
    np.testing.assert_array_equal(result, manual_result)
