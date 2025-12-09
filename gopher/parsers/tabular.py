"""Parse tabular result files from common tools"""

import os
import io
import pandas as pd
import numpy as np
from cloudpathlib import AnyPath


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


def _read_colnames(file: os.PathLike | io.TextIOBase) -> list[str]:
    with open(AnyPath(file)) as f:
        firstcol = f.readline()

    return firstcol.strip().split("\t")


def read_diann(proteins_tsv: os.PathLike) -> pd.DataFrame:
    """
    Reads a DIANN-generated TSV file (pg_matrix) containing protein information.

    Also processes it, and returns a cleaned Pandas DataFrame with relevant data.

    The function:
    - Extracts the first protein accession from the "Protein.Ids" column to use
        as the DataFrame index.
    - Renames the index axis to "Protein".
    - Drops unnecessary metadata columns.

    Args:
        proteins_tsv (str): Path to the DIANN-generated TSV file.
            Expected columns:
                'Protein.Group',
                'Protein.Ids',
                'Protein.Names',
                'Genes',
                'First.Protein.Description',
                <several Intensity columns>


    Returns:
        pd.DataFrame: A DataFrame with the processed protein data, indexed by
            the first protein accession.
            The returned DataFrame has the "Protein.Ids" column as the
            index and all columns are the MSR columns.
    """

    columns = _read_colnames(proteins_tsv)

    expect = [
        "Protein.Group",
        "Protein.Ids",
        "Protein.Names",
        "Genes",
        "First.Protein.Description",
    ]

    if not all(c in columns for c in expect):
        msg = f"Expected columns {expect}, got {columns}, make sure you are"
        msg += " using the 'diann_report.pg_matrix.tsv' output."
        raise ValueError(msg)

    schema: dict[str, type] = {k: float for k in columns if k not in expect}
    schema["Protein.Ids"] = str

    proteins = pd.read_table(
        AnyPath(proteins_tsv), dtype=schema, usecols=list(schema)
    )
    proteins["Protein.Ids"] = proteins["Protein.Ids"].str.split(";").str[0]

    proteins = proteins.set_index("Protein.Ids", drop=True)
    proteins = proteins.rename_axis("Protein", axis="index")

    return proteins
