"""Parse EncyclopeDIA results"""
import pandas as pd


def read_encyclopedia(proteins_txt):
    """Read results from EncyclopeDIA.

    Parameters
    ----------
    proteins_txt : str
        The EncyclopeDIA protein output.

    Returns
    -------
    pandas.DataFrame
        The EncyclopeDIA results in a format for gopher.
    """
    proteins = pd.read_table(proteins_txt)
    accessions = proteins["Protein"].str.extract(f"\|(.+?)\|", expand=False)

    proteins = proteins.set_index(accessions)
    return proteins.drop(
        columns=["Protein", "NumPeptides", "PeptideSequences"]
    )
