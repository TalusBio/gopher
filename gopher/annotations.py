"""Get GO annotations."""
from pathlib import Path
import requests
import pandas as pd
import polars as pl

from . import utils, config, ontologies
import uuid

SPECIES = {
    "yeast": "sgd",
    "saccharomyces cerevisiae": "sgd",
    "human": "goa_human",
    "homo sapiens": "goa_human",
}


def generate_annotations(proteins, aspect, go_name, go_id=None):
    """Generate an annotation file for a list of proteins that are correlated to a single term and aspect.

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
        String of the GO ID. If in the GO database, the go id and go name should match the database.

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
    if (not out_file.exists()) or fetch:
        out_file.parent.mkdir(exist_ok=True, parents=True)
        utils.http_download(url + fname, out_file)

    # Decompress file if it ends with gz
    target = Path(str(out_file).replace(".gz", ""))
    if not target.exists():
        out_file = utils.decompress(out_file, target)
    else:
        out_file = target

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

    keep = [
        "uniprot_accession",
        "go_id",
        "aspect",
        "db",
        "gene_product_form_id",
    ]
    # I decompressed the file manually just to see how much faster or efficient it would be
    polars_annot = pl.scan_csv(
        annot_file,
        separator="\t",
        comment_char="!",
        has_header=False,
        new_columns=cols,
    )
    polars_annot = polars_annot.select(keep)
    # polars_annot = pl.read_csv( annot_file, separator="\t", comment_char="!", has_header=False, new_columns=new_cols, columns=cols_keep)

    if aspect is not None:
        polars_annot = polars_annot.filter(pl.col("aspect") == aspect)

    polars_annot = polars_annot.with_columns(
        [
            pl.when(pl.col("db").eq("UniProtKB"))
            .then(pl.col("uniprot_accession"))
            .otherwise(
                pl.col("gene_product_form_id")
                .fill_null("")
                .str.extract(r"UniProtKB:(.+)", 1)
            )
            .alias("uniprot_accession")
        ],
    )

    polars_annot = (
        polars_annot.select(["uniprot_accession", "go_id", "aspect"])
        .unique()
        .with_columns(pl.col("go_id").map_dict(terms).alias("go_name"))
    )

    annot = polars_annot.collect().to_pandas()

    return annot, mapping


if __name__ == "__main__":
    load_annotations("human", aspect="all", release="current", fetch=False)
