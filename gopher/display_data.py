from .annotations import load_annotations
import pandas as pd
import logging
from .stats import rankdata
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns


def map_proteins(
    protein_list,
    aspect="all",
    species="human",
    release="current",
    fetch=False,
):
    """Map the proteins to the GO terms

    Parameters
    ----------
    protein_list : list
        A list of UniProt accessions.
    aspect : str, {"cc", "mf", "bp", "all"}, optional
        The Gene Ontology aspect to use. Use "cc" for "Cellular Compartment",
        "mf" for "Molecular Function", "bp" for "Biological Process", or "all"
        for all three.
    species : str, {"human", "yeast", ...}, optional.
        The species for which to retrieve GO annotations. If not "human" or
        "yeast", see
        [here](http://current.geneontology.org/products/pages/downloads.html).
    release : str, optional
        The Gene Ontology release version. Using "current" will look up the
        most current version.
    fetch : bool, optional
        Download the GO annotations even if they have been downloaded before?

    Returns
    -------
    pandas.DataFrame
        Dataframe with protein accessions and GO terms
    """
    # Put the protein list into a dataframe so it can be used to get the annotations
    proteins = pd.DataFrame(index=protein_list)
    # Get the annotations
    annot = get_annotations(proteins, aspect, species, release, fetch)
    # Return relevant columns
    return annot[["uniprot_accession", "go_id", "go_name"]]


def get_rankings(
    proteins,
    go_term,
    aspect="all",
    species="human",
    release="current",
    fetch=False,
):
    """Rank the proteins and show whether proteins are in a specified term

    Parameters
    ----------
    proteins : pd.DataFrame
        Dataframe of protein quant data
    go_term : str
        String of specified GO term name
    aspect : str, {"cc", "mf", "bp", "all"}, optional
        The Gene Ontology aspect to use. Use "cc" for "Cellular Compartment",
        "mf" for "Molecular Function", "bp" for "Biological Process", or "all"
        for all three.
    species : str, {"human", "yeast", ...}, optional.
        The species for which to retrieve GO annotations. If not "human" or
        "yeast", see
        [here](http://current.geneontology.org/products/pages/downloads.html).
    release : str, optional
        The Gene Ontology release version. Using "current" will look up the
        most current version.
    fetch : bool, optional
        Download the GO annotations even if they have been downloaded before?

    Returns
    -------
    pandas.DataFrame
        Dataframe with protein rankings and whether or not the protein is in the specified term
    """
    # Get the annotations
    annot = get_annotations(proteins, aspect, species, release, fetch)
    # Rank the data and format it as a dataframe
    ranked = rankdata(proteins.to_numpy())
    ranked = pd.DataFrame(ranked, columns=proteins.columns)
    ranked = ranked.set_index(proteins.index)
    # Add column for if proteins are in given term and return data
    ranked = in_term(ranked, go_term, annot).drop_duplicates()
    return ranked


def get_annotations(
    proteins,
    aspect="all",
    species="human",
    release="current",
    fetch=False,
    go_subset=None,
):
    """Gets the annotations for proteins in a dataset

    Parameters
    ----------
    proteins : pd.DataFrame
        Dataframe of proteins and quantifications
    aspect : str, {"cc", "mf", "bp", "all"}, optional
        The Gene Ontology aspect to use. Use "cc" for "Cellular Compartment",
        "mf" for "Molecular Function", "bp" for "Biological Process", or "all"
        for all three.
    species : str, {"human", "yeast", ...}, optional.
        The species for which to retrieve GO annotations. If not "human" or
        "yeast", see
        [here](http://current.geneontology.org/products/pages/downloads.html).
    release : str, optional
        The Gene Ontology release version. Using "current" will look up the
        most current version.
    fetch : bool, optional
        Download the GO annotations even if they have been downloaded before?
    go_subset: list of str, optional
        The go terms of interest. Should consists of the go term names such
        as 'nucleus' or 'cytoplasm'.
    Returns
    -------
    pandas.DataFrame
        Dataframe with protein annotations
    """
    # Load the annotation file
    annot = load_annotations(
        species=species,
        aspect=aspect,
        release=release,
        fetch=fetch,
    )

    # If a subset of terms was given, filter for the subset
    if go_subset:
        in_names = annot["go_name"].isin(go_subset)
        in_ids = annot["go_id"].isin(go_subset)
        annot = annot.loc[in_names | in_ids, :]

    # Filter for proteins seen in input dataframe
    accessions = pd.DataFrame(
        list(proteins.index),
        columns=["uniprot_accession"],
    )
    annot = accessions.merge(annot, how="inner")
    return annot


def in_term(proteins, go_term, annot):
    """See if proteins are associated with a specific term

    Parameters
    ----------
    proteins : pd.DataFrame
        Dataframe of proteins and quantifications
    go_term : str
        String of specified GO term name
    annot : pd.DataFrame
        Annotation file for the dataset

    Returns
    -------
    pandas.DataFrame
        Dataframe with protein quant and if protein is in the given term
    """
    # Get columns of terms
    term = annot[annot["go_name"] == go_term]
    proteins["in_term"] = proteins.index.isin(
        term["uniprot_accession"].unique()
    )
    return proteins


def roc(
    proteins,
    go_term,
    aspect="all",
    species="human",
    release="current",
    fetch=False,
):
    """Plot the ROC curve for a go term in each sample

    Parameters
    ----------
    proteins : pd.DataFrame
        Dataframe of proteins and quantifications
    go_term : str
        String of specified GO term name
    aspect : str, {"cc", "mf", "bp", "all"}, optional
        The Gene Ontology aspect to use. Use "cc" for "Cellular Compartment",
        "mf" for "Molecular Function", "bp" for "Biological Process", or "all"
        for all three.
    species : str, {"human", "yeast", ...}, optional.
        The species for which to retrieve GO annotations. If not "human" or
        "yeast", see
        [here](http://current.geneontology.org/products/pages/downloads.html).
    release : str, optional
        The Gene Ontology release version. Using "current" will look up the
        most current version.
    fetch : bool, optional
        Download the GO annotations even if they have been downloaded before?
    Returns
    -------
    plot
        Plot of ROC curve for a GO term
    """
    # Get a list of the samples
    samples = proteins.columns
    # Get annotations
    annot = get_annotations(proteins, aspect, species, release, fetch)
    # Rank the data based on the go term
    proteins = in_term(proteins, go_term, annot)
    # Get the number of positive and negative cases
    n_pos = sum(proteins["in_term"])
    n_neg = len(proteins) - n_pos

    # Set up plot
    fig, axs = plt.subplots(1, len(samples), figsize=(13, 4.5))
    i = 0

    # Graph the ROC curve for each sample
    for sample in samples:

        # Sort values
        sorted = proteins.sort_values(sample, ascending=False)

        # Calculate TPR and FPR
        sorted["tpr"] = sorted["in_term"].cumsum() / sorted["in_term"].sum()
        sorted["fpr"] = (~sorted["in_term"]).cumsum() / (
            ~sorted["in_term"]
        ).sum()

        # Need to add a point at 0, 0:
        tpr = np.append(0, sorted["tpr"])
        fpr = np.append(0, sorted["fpr"])

        # Graph ROC curve
        ax = axs[i]
        p = sns.lineplot(x=fpr, y=tpr, drawstyle="steps-pre", ci=None, ax=ax)
        sns.lineplot(
            x=(0, 1), y=(0, 1), color="black", linestyle="dashed", ax=ax
        )
        ax.title.set_text(sample)
        ax.set_xlabel("False Positive Rate (FPR)")
        ax.set_ylabel("True Positive Rate (TPR)")

        # Calculate AUC and put on graph in lower right corner
        U, _ = stats.mannwhitneyu(
            proteins.loc[proteins["in_term"], sample],
            proteins.loc[~proteins["in_term"], sample],
        )
        auc = U / (n_pos * n_neg)
        auc = "AUC = " + str(round(auc, 3))
        p.annotate(auc, xy=(0.75, 0))

        i += 1

    # Format and return plot
    for ax in axs:
        ax.set(adjustable="box", aspect="equal")

    plt.tight_layout()

    return plt
