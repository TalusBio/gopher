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
    # If the proteins.txt includes multiple proteins in the protein column
    # encyclopedia formats it like this: sp|Q5T1J5|CHCH9_HUMAN;sp|Q9Y6H1|CHCH2_HUMAN
    # Therefore we need to explode the column after splitting it by ';'.
    proteins = proteins.assign(
        **{
            "Protein": proteins["Protein"].apply(
                lambda row: row.split(";") if isinstance(row, str) else row
            )
        }
    ).explode(column="Protein", ignore_index=True)
    accessions = proteins["Protein"].str.extract(f"\|(.+?)\|", expand=False)

    proteins = proteins.set_index(accessions)
    return proteins.drop(
        columns=["Protein", "NumPeptides", "PeptideSequences"], axis=1
    )
