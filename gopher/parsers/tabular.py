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
        The Metamorpheus results in a format for gopher.
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

def read_diann(proteins_tsv: str) -> pd.DataFrame:
    """
    Reads a DIANN-generated TSV file containing protein information, processes
    it, and returns a cleaned Pandas DataFrame with relevant data.

    The function:
    - Extracts the first protein accession from the "Protein.Ids" column to use
        as the DataFrame index.
    - Renames the index axis to "Protein".
    - Drops unnecessary metadata columns.

    Args:
        proteins_tsv (str): Path to the DIANN-generated TSV file.

    Returns:
        pd.DataFrame: A DataFrame with the processed protein data, indexed by
            the first protein accession.
            The returned DataFrame excludes the following columns:
            ["Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
            "First.Protein.Description"].
    """
    proteins = pd.read_table(proteins_tsv)
    accessions = proteins["Protein.Ids"].str.split(";").str[0]

    proteins = proteins.set_index(accessions)
    proteins = proteins.rename_axis("Protein", axis="index")
    return proteins.drop(
        columns=[
            "Protein.Group",
            "Protein.Ids",
            "Protein.Names",
            "Genes",
            "First.Protein.Description",
        ]
    )