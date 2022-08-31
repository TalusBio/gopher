import pytest
from pathlib import Path
import pandas as pd
from gopher import normalize


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
    df = normalize.normalize(quant, fasta)
    assert isinstance(df, pd.DataFrame)
