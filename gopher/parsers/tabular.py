"""Parse tabular result files from common tools"""
import pandas as pd


def read_encyclopedia(proteins_txt: str) -> pd.DataFrame:
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
    accessions = proteins["Protein"].str.extract(r"\|(.+?)\|", expand=False)

    proteins = proteins.set_index(accessions)
    return proteins.drop(
        columns=["Protein", "NumPeptides", "PeptideSequences"]
    )


def read_metamorpheus(proteins_txt: str) -> pd.DataFrame:
    """Read results from Metamorpheus.

    Parameters
    ----------
    proteins_txt : str
        The Metamorpheus protein output file.

    Returns
    -------
    pandas.DataFrame
    """
    proteins = pd.read_table(proteins_txt, low_memory=False)
    accessions = proteins["Protein"].str.extract(
        r"^(.*?)(\||$)", expand=False
    )[0]
    accessions.name = "Protein"
    int_cols = [c for c in proteins.columns if c.startswith("Intensity")]
    proteins = proteins.set_index(accessions)
    proteins = (
        proteins.loc[
            proteins["Protein Decoy/Contaminant/Target"] == "T", int_cols
        ]
        .apply(pd.to_numeric)
        .fillna(0)
    )
    return proteins
