"""Calculate the enrichments for a collection of experiments."""
import logging

import pandas as pd
from statsmodels.stats import multitest
from tqdm.auto import tqdm
from .stats import mannwhitneyu
from .tree_search import tree_search

from .annotations import load_annotations

LOGGER = logging.getLogger(__name__)


def test_enrichment(
    proteins,
    desc=True,
    aspect="all",
    species="human",
    release="current",
    go_subset=None,
    contaminants_filter=None,
    fetch=False,
    progress=False,
    annotations=None,
    mapping=None,
    aggregate_terms=True,
):
    """Test for the enrichment of Gene Ontology terms from protein abundance.

    The Mann-Whitney U Test is applied to each column of proteins dataframe and
    for each Gene Ontology (GO) term. The p-values are then corrected for
    multiple hypothesis testing across all of the columns using the
    Benjamini-Hochberg procedure.

    Parameters
    ----------
    proteins : pandas.DataFrame
        A dataframe where the indices are UniProt accessions and each column is
        an experiment to test. The values in this dataframe should be some
        measure of protein abundance: these could be the raw measurement if
        originating from a single sample or a fold-change/p-value if looking at
        the difference between two conditions.
    desc : bool, optional
        Rank proteins in descending order?
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
    go_subset: list of str, optional
        The go terms of interest. Should consists of the go term names such
        as 'nucleus' or 'cytoplasm'.
    contaminants_filter: List[str], optional
        A list of uniprot accessions for common contaminants such as
        Keratin to filter out.
    fetch : bool, optional
        Download the GO annotations even if they have been downloaded before?
    progress : bool, optional
        Show a progress bar during enrichment tests?
    annotations: pandas.DataFrame, optional
        A custom annotations dataframe.
    mapping: defaultdict, optional
        A custom mapping of the GO term relationships.
    aggregate_terms : bool, optional
        Aggregate the terms and do the tree search.
    Returns
    -------
    pandas.DataFrame
        The adjusted p-value for each tested GO term in each sample.
    """
    LOGGER.info("Retrieving GO annotations...")

    if annotations is not None:
        annot = annotations
    else:
        annot, map = load_annotations(
            species=species,
            aspect=aspect,
            release=release,
            fetch=fetch,
        )
        if not mapping:
            mapping = map

    if go_subset:
        if aggregate_terms and mapping:
            annot = tree_search(mapping, go_subset, annot)

        in_names = annot["go_name"].isin(go_subset)
        in_ids = annot["go_id"].isin(go_subset)
        annot = annot.loc[in_names | in_ids, :]

    accessions = pd.DataFrame(
        list(proteins.index),
        columns=["uniprot_accession"],
    )
    if contaminants_filter:
        accessions = accessions[
            ~accessions["uniprot_accession"].isin(contaminants_filter)
        ]

    # Get the GO terms and proteins
    annot = accessions.merge(annot, how="inner")
    n_prot = proteins.shape[1]
    proteins = pd.DataFrame(proteins).loc[annot["uniprot_accession"], :]
    lost = n_prot - proteins.shape[1]
    if lost:
        LOGGER.warning("%i proteins not found in GO annotations.", lost)

    if not desc:
        proteins = -proteins

    results = []
    grp_cols = ["go_id", "go_name", "aspect"]

    LOGGER.info("Testing enrichment...")
    for term, accessions in tqdm(
        annot.groupby(grp_cols), disable=not progress
    ):

        in_term = proteins.index.isin(accessions["uniprot_accession"].unique())
        in_vals = proteins[in_term].to_numpy()
        out_vals = proteins[~in_term].to_numpy()
        res = mannwhitneyu(in_vals, out_vals, alternative="greater")
        if res != None:
            results.append(list(term) + list(res[1]))

    cols = ["GO ID", "GO Name", "GO Aspect"] + list(proteins.columns)
    results = pd.DataFrame(results, columns=cols)
    results.loc[:, proteins.columns] = results.loc[:, proteins.columns].apply(
        adjust_pvals, raw=True
    )

    return results


def adjust_pvals(pvals):
    """Compute BH adjusted p-values.

    Paramerters
    -----------
    pvals : numpy.ndarray
        A 1D numpy array of p-values.

    Returns
    -------
    numpy.ndarray
        The FDR adjusted p-values.
    """
    return multitest.fdrcorrection(pvals)[1]
