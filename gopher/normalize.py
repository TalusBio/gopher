import pandas as pd
from Bio import SeqIO


def normalize(proteins, fasta):
    """Normalize intensity values.
    Parameters
    ----------
    proteins : pandas.DataFrame
        A dataframe where the indices are UniProt accessions and each column is
        an experiment to test. The values in this dataframe should be some
        measure of protein abundance: these could be the raw measurement if
        originating from a single sample or a fold-change/p-value if looking at
        the difference between two conditions.
    fasta : fasta file
        Use the FASTA file to generate molecular weights for normalization
    Returns
    -------
    pandas.DataFrame
        The normalized intensities for every protein in each sample.
    """
    fasta_df = pd.DataFrame(columns=["Protein", "Sequence", "Mass"])
    fasta_sequences = SeqIO.parse(open(fasta), "fasta")
    for f in fasta_sequences:
        name = f.id.split("|")[1]
        mass = calculate_mass(f.seq)
        temp = {"Protein": name, "Sequence": str(f.seq), "Mass": mass}
        temp = pd.DataFrame(temp, index=[0])
        fasta_df = pd.concat([fasta_df, pd.DataFrame(temp)], ignore_index=True)

    proteins = proteins.apply(col_norm, axis=0)
    df = pd.merge(fasta_df, proteins, on="Protein")
    df = df.set_index("Protein")
    df = df.drop(columns=["Sequence"])
    df = df.apply(mass_norm, axis=1)
    df = df.drop(columns=["Mass"])

    return df


def col_norm(col):
    """Calculate the ratio of protein/total protein"""
    return col / col.sum()


def mass_norm(row):
    """Do mass part of normalization"""
    mass = row.loc["Mass"]
    return row / mass


def calculate_mass(seq):
    """Calculate the weights of each protein by its sequence"""
    weights = {
        "A": 89.10,
        "R": 174.20,
        "N": 132.12,
        "D": 133.11,
        "C": 121.16,
        "E": 147.13,
        "Q": 146.15,
        "G": 75.07,
        "H": 155.16,
        "I": 131.18,
        "L": 131.18,
        "K": 146.19,
        "M": 149.21,
        "F": 165.19,
        "P": 115.13,
        "S": 105.09,
        "T": 119.12,
        "W": 204.23,
        "Y": 181.19,
        "V": 117.15,
    }

    mw = 0
    for aa in seq:
        mw += weights[aa]

    return mw
