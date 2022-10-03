from gopher import tree_search
from gopher import enrichment


def test_tree_mapping(generate_mapping):
    """Test the tree mapping functionality maps the terms and returns the new mapping in a dictionary."""
    mapping = generate_mapping
    go_subset = ["a", "y", "z"]
    mapped = tree_search.new_tree_map(mapping, go_subset)
    assert isinstance(mapped, dict)


def test_tree_search(generate_annotations, generate_mapping):
    """Test the tree search functionality adds the information to the file before returning it."""
    annot = generate_annotations
    mapping = generate_mapping
    go_subset = ["a", "y", "z"]
    result = tree_search.tree_search(mapping, go_subset, annot)
    assert len(result) > len(annot)


def test_enrichment_tree_search(generate_proteins):
    """Check that the test enrichment function works with the incorporated tree algorithm
    and gives you different results from the enrichment without the tree algorithm."""
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
        df, go_subset=terms, tree_alg=True
    )
    assert not result_orig.equals(result_tree_search)
