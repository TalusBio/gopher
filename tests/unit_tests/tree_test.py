import pandas as pd

from gopher import tree_search
from gopher import enrichment


def test_tree_mapping(generate_mapping):
    """Test the tree mapping functionality maps the terms and returns the new mapping in a dictionary."""
    mapping = generate_mapping
    go_subset = ["a", "y", "z"]
    mapped = tree_search.new_tree_map(mapping, go_subset)
    mapped_true = {
        "a": ["b", "c", "d", "e", "f", "g", "h"],
        "y": ["z", "l", "m", "n"],
        "z": ["l", "m", "n"],
    }
    assert mapped == mapped_true


def test_tree_search(generate_annotations, generate_mapping):
    """Test the tree search algorithm gives the information that would be expected."""
    annot = generate_annotations
    mapping = generate_mapping
    go_subset = ["y", "z"]
    new = [
        [25, "y", "c", "Y"],
        [11, "y", "c", "Y"],
        [12, "y", "c", "Y"],
        [13, "y", "c", "Y"],
        [11, "z", "c", "Z"],
        [12, "z", "c", "Z"],
        [13, "z", "c", "Z"],
    ]
    additional_df = pd.DataFrame(
        new, columns=["uniprot_accession", "go_id", "aspect", "go_name"]
    )
    manual_result = pd.concat([annot, additional_df], ignore_index=True)
    result = tree_search.tree_search(mapping, go_subset, annot)
    assert result.equals(manual_result)
    assert len(result) > len(annot)


def test_enrichment_tree_search(generate_proteins):
    """Check that the test enrichment function works with the incorporated tree algorithm
    and gives you different results from the enrichment without the tree algorithm.
    """
    df = generate_proteins
    df.set_index("Protein", inplace=True)
    terms = [
        "nucleus",
        "nuclear chromosome",
        "nucleoplasm",
        "euchromatin",
        "heterochromatin",
        "cytoplasm",
    ]
    result_orig = enrichment.test_enrichment(df, go_subset=terms)
    result_tree_search = enrichment.test_enrichment(
        df, go_subset=terms, aggregate_terms=True
    )
    assert result_orig.equals(result_tree_search)
