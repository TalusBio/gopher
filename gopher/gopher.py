"""The command line entry point for gopher-enrich"""
import logging
from argparse import ArgumentParser


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
        The quanitfied proteins in each sample. Currently, only results from
        EncyclopeDIA are supported.
        """,
    )

    parser.add_argument(
        "-a",
        "--aspect",
        choices=["cc", "mf", "bp", "all"],
        default="all",
        help="""
        The Gene Ontology aspect to use. Use "cc" for  "Cellular Compartment",
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
        The species for which to retrieve GO annotations
        """,
    )


def main():
    """The command line function"""


if __name__ == "__main__":
    main()
