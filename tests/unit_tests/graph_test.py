import pandas as pd

from gopher import enrichment, graph_search


def test_graph_mapping(generate_mapping):
    """Ensure graph mapping expands terms and returns the updated mapping."""
    mapping = generate_mapping
    go_subset = ["a", "y", "z"]
    mapped = graph_search.new_map(mapping, go_subset)
    mapped_true = {
        "a": ["b", "c", "d", "e", "f", "g", "h"],
        "y": ["z", "l", "m", "n"],
        "z": ["l", "m", "n"],
    }
    assert mapped == mapped_true


def test_graph_search(generate_annotations, generate_mapping):
    """Ensure graph search returns expected term expansions."""
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
    result = graph_search.graph_search(mapping, go_subset, annot)
    assert result.equals(manual_result)
    assert len(result) > len(annot)


def test_enrichment_graph_search(generate_proteins):
    """Check enrichment works with graph algorithm and changes the results."""
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
    result_orig = enrichment.test_enrichment(
        df, go_subset=terms, aggregate_terms=False
    )
    result_graph_search = enrichment.test_enrichment(
        df, go_subset=terms, aggregate_terms=True
    )
    assert not result_orig.equals(result_graph_search)
