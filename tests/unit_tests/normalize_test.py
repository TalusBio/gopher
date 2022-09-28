import pandas as pd
from gopher import normalize


def test_normalize(real_data):
    """Check that the normalization returns a dataframe."""
    quant, fasta = real_data
    df = normalize.normalize(quant, fasta)
    assert isinstance(df, pd.DataFrame)
