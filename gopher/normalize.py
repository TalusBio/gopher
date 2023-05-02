import pandas as pd
from Bio import SeqIO, SeqUtils


def normalize_values(proteins, fasta):
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

    proteins = proteins.apply(col_norm, axis=0)
    df = pd.merge(fasta_df, proteins, on="Protein")
    df = df.set_index("Protein")
    df = df.drop(columns=["Sequence"])
    df = df.apply(mass_norm, axis=1)
    df = df.drop(columns=["Mass"])

    return df


def read_fasta(fasta):
    fasta_df = []
    for entry in SeqIO.parse(open(fasta), "fasta"):
        name = entry.id.split("|")[1]
        mass = SeqUtils.molecular_weight(entry.seq, seq_type="protein")
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
