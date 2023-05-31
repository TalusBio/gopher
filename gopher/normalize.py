import pandas as pd
from Bio import SeqIO, SeqUtils
from loguru import logger
from pathlib import Path


def normalize_values(proteins: pd.DataFrame, fasta: Path):
    """Normalize intensity values.

    Normalize using the proteomic ruler approach outlined by Wiśniewski et al.
    (doi: https://doi.org/10.1074/mcp.M113.037309)

    Parameters
    ----------
    proteins : pandas.DataFrame
        A dataframe where the indices are UniProt accessions and each column is
        an experiment to test. The values in this dataframe raw protein abundance
    fasta : Path
        Use the FASTA file to generate molecular weights for normalization

    Returns
    -------
    pandas.DataFrame
        The normalized intensities for every protein in each sample.
    """
    fasta_df = read_fasta(fasta)
    fasta_df = fasta_df.set_index("Protein")
    proteins = proteins.apply(col_norm, axis=0)
    df = proteins.join(fasta_df)
    df = df.drop(columns=["Sequence"])
    df = df.apply(mass_norm, axis=1)
    df = df.drop(columns=["Mass"])

    return df


def read_fasta(fasta: Path):
    """Read a fasta file into a dataframe with relevant columns

    Parameters
    ----------
    fasta : Path
        Use the FASTA file to read and generate molecular weights

    Returns
    -------
    pandas.DataFrame
        The fasta dataframe with relevant columns.
    """
    fasta_df = []
    for entry in SeqIO.parse(open(fasta), "fasta"):
        name = entry.id.split("|")[1]
        try:
            mass = SeqUtils.molecular_weight(entry.seq, seq_type="protein")
        except:
            logger.warning("Ambiguous peptide in {}".format(name))
        temp = pd.DataFrame(
            {"Protein": name, "Sequence": str(entry.seq), "Mass": mass},
            index=[0],
        )
        fasta_df.append(temp)

    fasta_df = pd.concat(fasta_df)
    return fasta_df


def col_norm(col):
    """Calculate the ratio of protein/total protein"""
    return col / col.sum()


def mass_norm(row):
    """Do mass part of normalization"""
    mass = row.loc["Mass"]
    return row / mass
