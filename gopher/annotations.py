"""Get GO annotations."""
import requests
import pandas as pd

from . import utils, config, ontologies

SPECIES = {
    "yeast": "sgd",
    "saccharomyces cerevisiae": "sgd",
    "human": "goa_human",
    "homo sapiens": "goa_human",
}


def download_annotations(stem, release="current", fetch=False):
    """Download the annotation file.

    See http://current.geneontology.org/annotations/index.html for details.

    Parameters
    ----------
    stem : str
        The stem of the annotation file name to retrieve.
    release: str
        The Gene Ontology release version. Using "current" will look up the
        most current version.
    fetch : bool
        Download the file even if it already exists?

    Returns
    -------
    Path
        The path to downloaded file.
    """
    fname = stem.split(".")[0] + ".gaf.gz"
    if release == "current":
        base = "http://current.geneontology.org/"
        meta_url = base + "metadata/release-date.json"
        release = requests.get(meta_url).json()["date"]
        url = base + "annotations/"
    else:
        url = f"http://release.geneontology.org/{release}/annotations/"

    out_file = config.get_data_dir() / "annotations" / release / fname
    if out_file.exists() and not fetch:
        return out_file

    out_file.parent.mkdir(exist_ok=True, parents=True)
    utils.http_download(url + fname, out_file)
    return out_file


def load_annotations(species, aspect=None, release="current", fetch=False):
    """Load the Gene Ontology (GO) annotations for a species.

    Parameters
    ----------
    species : str, {"human", "yeast", ...}
        The species for which to retrieve GO annotations. If not "humnan" or
        "yeast", see
        [here](http://current.geneontology.org/products/pages/downloads.html).
    aspect : str, {"c", "f", "p"} or None
        The Gene Ontology aspect to use. Use "c" for "Cellular Compartment",
        "f" for "Molecular Function", or "p" for "Biological Process". ``None``
        uses all of the them.
    release: str
        The Gene Ontology release version. Using "current" will look up the
        most current version.
    fetch : bool
        Download the file even if it already exists?

    Returns
    -------
    dict of str: list of str
        A mapping of GO terms (keys) to Uniprot accessions with that
        annotation.
    """
    aspect = None if aspect is None else aspect.upper()
    if aspect not in {"C", "F", "P", None}:
        raise ValueError(
            f"Expected apsect ({aspect}) to be one of 'c', 'f', 'p', or None."
        )

    cols = [
        "db",
        "uniprot_accession",
        "db_object_symbol",
        "qualifier",
        "go_id",
        "db_reference",
        "evidence_code",
        "with_or_from",
        "aspect",
        "db_object_name",
        "db_object_synonym",
        "db_object_type",
        "taxon",
        "date",
        "assigned_by",
        "annotation_extension",
        "gene_product_form_id",
    ]

    terms = ontologies.load_ontology()
    species = SPECIES.get(species.lower(), species.lower())
    annot_file = download_annotations(species, release=release, fetch=fetch)
    annot = pd.read_table(
        annot_file,
        comment="!",
        header=None,
        names=cols,
        low_memory=False,
    )

    if aspect is not None:
        annot = annot.loc[annot["aspect"] == aspect, :]

    uniprot = annot["db"] == "UniProtKB"
    if not uniprot.all():
        annot.loc[~uniprot, "uniprot_accession"] = annot.loc[
            ~uniprot, "gene_product_form_id"
        ].str.extract(r"UniProtKB:(.+)", expand=False)

    keep = ["uniprot_accession", "go_id", "aspect"]
    annot = annot.loc[:, keep].drop_duplicates()
    annot["go_name"] = annot["go_id"].map(terms)
    return annot
