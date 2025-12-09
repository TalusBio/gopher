"""Parse tabular result files from common tools."""

from io import StringIO
from pathlib import Path

import pandas as pd
from cloudpathlib import CloudPath


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


def _load_table(path: str) -> pd.DataFrame:
    """Load a DIANN-style TSV from local or cloud storage."""
    path_str = str(path)
    cloud_prefixes = ("s3://", "gs://", "az://", "http://", "https://")
    if path_str.startswith(cloud_prefixes):
        cp = CloudPath(path_str)
        text = cp.read_text()
        return pd.read_csv(StringIO(text), sep="\t")

    return pd.read_csv(Path(path_str), sep="\t")


def read_diann(path: str) -> pd.DataFrame:
    """Read DIA-NN protein report.

    Parameters
    ----------
    path : str
        Local path or S3 URL to a DIA-NN `*.pg_mat.tsv` report.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by Protein accession with only intensity columns.

    Raises
    ------
    ValueError
        If required intensity columns are not present.

    """
    proteins = _load_table(path)

    intensity_cols = [
        c for c in proteins.columns if c.startswith("Intensity.")
    ]
    if not intensity_cols or "Protein.Ids" not in proteins.columns:
        raise ValueError(
            "Expected columns beginning with 'Intensity.' and 'Protein.Ids' "
            "in DIA-NN report."
        )

    accessions = proteins["Protein.Ids"].str.split(";").str[0]
    accessions.name = "Protein"
    proteins = proteins.set_index(accessions)
    proteins = proteins[intensity_cols].apply(pd.to_numeric).astype(float)
    return proteins
