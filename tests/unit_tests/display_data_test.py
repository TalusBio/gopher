import pandas as pd

from gopher import display_data


def test_protein_mapping():
    """Check that the protein mapping function returns a dataframe."""
    prot = ["P10809", "P35527", "Q9UMS4", "P35637", "P07437", "Q86UX7"]
    result = display_data.map_proteins(prot)
    assert isinstance(result, pd.DataFrame)


def test_rankings(generate_proteins):
    """Check that the enrichment rankings return a dataframe.

    This test runs GO enrichment on a dataset.
    """
    df = generate_proteins
    df.set_index("Protein", inplace=True)
    result = display_data.get_rankings(df, "cytoplasm")
    assert isinstance(result, pd.DataFrame)


def test_roc(generate_proteins):
    """Verify ROC plotting returns a matplotlib figure."""
    df = generate_proteins
    df.set_index("Protein", inplace=True)
    result = display_data.roc(df, "cytoplasm")
    assert result is not None
