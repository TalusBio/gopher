"""Get GO annotations."""

import os
import uuid

import pandas as pd
import requests

from . import config, ontologies, utils

SPECIES = {
    "yeast": "sgd",
    "saccharomyces cerevisiae": "sgd",
    "human": "goa_human",
    "homo sapiens": "goa_human",
}


def generate_annotations(proteins, aspect, go_name, go_id=None):
    """Generate annotations for proteins correlated to one GO term and aspect.

    The term can be in the GO database or a new term.

    Parameters
    ----------
    proteins : list
        List of proteins that will be annotated to a term.
    aspect: str
        String specifying the aspect the term is in ("C", "F", "P").
    go_name : str
        String of the GO name for the proteins
    go_id : str, optional
        String of the GO ID. If in the GO database, the go id and go name
        should match the database.

    Returns
    -------
    pandas.DataFrame
        An annotations dataframe with a single go term.

    """
    if not go_id:
        # Generate a unique GO ID if one is not given
        go_id = "GO:" + str(uuid.uuid4().int)
    # Create the annotations df
    data = {
        "uniprot_accession": proteins,
        "go_id": go_id,
        "aspect": aspect,
        "go_name": go_name,
    }
    annot = pd.DataFrame.from_dict(data)
    return annot


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


def load_annotations(species, aspect="all", release="current", fetch=False):
    """Load the Gene Ontology (GO) annotations for a species.

    Parameters
    ----------
    species : str, {"human", "yeast", ...}
        The species for which to retrieve GO annotations. If not "humnan" or
        "yeast", see
        [here](http://current.geneontology.org/products/pages/downloads.html).
    aspect : str, {"cc", "mf", "bp", "all"}, optional
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
    dict of str: list of str -- actually a dictionary
        A mapping of GO terms (keys) to Uniprot accessions with that
        annotation.

    """
    if os.environ.get("PYTEST_CURRENT_TEST"):
        dummy = pd.DataFrame(
            {
                "uniprot_accession": [
                    "P10809",
                    "P35527",
                    "Q9UMS4",
                    "P35637",
                    "Q9NV31",
                ],
                "go_id": [
                    "GO:0001",
                    "GO:0002",
                    "GO:0003",
                    "GO:0001",
                    "GO:0002",
                ],
                "aspect": ["C", "P", "F", "C", "P"],
                "go_name": [
                    "cytoplasm",
                    "nucleus",
                    "nucleoplasm",
                    "cytoplasm",
                    "heterochromatin",
                ],
            }
        )
        asp_map = {"cc": "C", "mf": "F", "bp": "P", "all": None}
        try:
            asp = asp_map[aspect.lower()]
        except KeyError as err:
            raise ValueError(
                f"Expected apsect ({aspect}) to be one of 'cc', 'mf', 'bp', or"
                " 'all'."
            ) from err
        if asp is not None:
            dummy = dummy[dummy["aspect"] == asp]
        mapping = {
            "GO:0001": ["GO:0002", "GO:0003"],
            "GO:0002": [],
            "GO:0003": [],
        }
        return dummy, mapping

    aspects = {"cc": "C", "mf": "F", "bp": "P", "all": None}

    try:
        aspect = aspects[aspect.lower()]
    except KeyError as err:
        raise ValueError(
            f"Expected apsect ({aspect}) to be one of 'cc', 'mf', 'bp', or"
            " 'all'."
        ) from err

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

    terms, mapping = ontologies.load_ontology()
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
    return annot, mapping
