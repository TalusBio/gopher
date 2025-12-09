"""The command line entry point for gopher-enrich."""

import logging
from argparse import ArgumentParser

from .enrichment import test_enrichment
from .parsers import read_encyclopedia

LOGGER = logging.getLogger(__name__)


def parse_args():
    """Get the command line arguments.

    Returns
    -------
    Namespace
        A namespace populated with the parsed arguments.

    """
    desc = """
    gopher: Gene ontology enrichment analysis using protein expression. For
     more details see TalusBio.github.io/gopher
    """
    parser = ArgumentParser(description=desc)

    parser.add_argument(
        "proteins",
        type=str,
        help="""
        The quantified proteins in each sample. Currently, only results from
         EncyclopeDIA are supported.
        """,
    )

    parser.add_argument(
        "-o",
        "--output",
        help="The name of the tab-delimited output file.",
    )

    parser.add_argument(
        "-a",
        "--aspect",
        choices=["cc", "mf", "bp", "all"],
        default="all",
        help="""
        The Gene Ontology aspect to use. Use "cc" for "Cellular Compartment",
         "mf" for "Molecular Function", and "bp" for "Biological Process", or
         "all" for all three.
        """,
    )

    parser.add_argument(
        "-s",
        "--species",
        type=str,
        default="human",
        help="""
        The species for which to retrieve GO annotations. If not "human" or
         "yeast", see http://current.geneontology.org/products/pages/
        downloads.html
        """,
    )

    parser.add_argument(
        "-r",
        "--release",
        type=str,
        default="current",
        help="""
        The Gene Ontology release version. Using "current" will look up the
         most current version.
        """,
    )

    parser.add_argument(
        "-g",
        "--go_filters",
        type=str,
        help="""
        The GO terms of interest, separated by commas (","). May consist of
         GO term names (ex: "nucleus,cytosol") or GO term accessions (ex:
         "GO:0005634,GO:0005829").
        """,
    )

    parser.add_argument(
        "-f",
        "--fetch",
        action="store_true",
        help="""
        Download the GO annotations even if they have been downloaded before.
        """,
    )

    parser.add_argument(
        "-p",
        "--progress",
        action="store_true",
        help="""
        Show a progress bar during enrichment tests.
        """,
    )

    return parser.parse_args()


def main():
    """The main command line function."""
    logging.basicConfig(
        level=logging.INFO, format="[%(levelname)s] %s(message)s"
    )

    args = parse_args()
    proteins = read_encyclopedia(args.proteins)
    if args.go_filters is not None:
        args.go_filters = args.go_filters.split(",")

    results = test_enrichment(
        proteins=proteins,
        desc=True,
        aspect=args.aspect,
        species=args.species,
        release=args.release,
        go_filters=args.go_filters,
        fetch=args.fetch,
        progress=args.progress,
    )

    results.to_csv(args.output, index=False, sep="\t")


if __name__ == "__main__":
    main()
