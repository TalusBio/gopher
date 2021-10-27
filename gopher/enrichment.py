"""Calculate the enrichments for a collection of experiments."""
import logging

import pandas as pd
from scipy import stats
from statsmodels.stats import multitest
from tqdm.auto import tqdm

from .annotations import load_annotations

LOGGER = logging.getLogger(__name__)


def test_enrichment(
    proteins,
    desc=True,
    aspect=None,
    species="human",
    release="current",
    fetch=False,
    go_filters=None,
    progress=True,
):
    """Test for the enrichment of Gene Ontology terms from protein abundance.

    The Mann-Whitney U Test is applied to each column of proteins dataframe and
    for each Gene Ontology (GO) term. The p-values are then corrected for
    multiple hypothesis testing accross all of the columns using the
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
    aspect : str, {"c", "f", "p"} or None, optional
        The Gene Ontology aspect to use. Use "c" for "Cellular Compartment",
        "f" for "Molecular Function", or "p" for "Biological Process". ``None``
        uses all of the them.
    species : str, {"human", "yeast", ...}, optional.
        The species for which to retrieve GO annotations. If not "humnan" or
        "yeast", see
        [here](http://current.geneontology.org/products/pages/downloads.html).
    release : str, optional
        The Gene Ontology release version. Using "current" will look up the
        most current version.
    go_filters: List[str], optional
        The go terms of interest. Should consists of the go term names such
        as 'nucleus' or 'cytoplasm'.
    fetch : bool, optional
        Download the annotations even if the file already exists?

    Returns
    -------
    pandas.DataFrame
        The adjusted p-value for each tested GO term in each sample.
    """
    annot = load_annotations(
        species=species,
        aspect=aspect,
        release=release,
        fetch=fetch,
    )

    if go_filters:
        annot = annot[annot["go_name"].isin(go_filters)]

    accessions = pd.DataFrame(
        list(proteins.index),
        columns=["uniprot_accession"],
    )
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
    for term, accessions in tqdm(
        annot.groupby(grp_cols), disable=not progress
    ):
        in_term = proteins.index.isin(accessions["uniprot_accession"].unique())
        in_vals = proteins[in_term].values
        out_vals = proteins[~in_term].values
        res = stats.mannwhitneyu(in_vals, out_vals, alternative="greater")
        results.append(list(term) + list(res[1]))

    cols = ["GO Accession", "GO Name", "GO Aspect"] + list(proteins.columns)
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
