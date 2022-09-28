from gopher import tree_search
from gopher import annotations


def test_tree_search():
    mapping = {
        "a": ["b", "c", "d"],
        "b": ["e"],
        "e": ["f"],
        "d": ["g", "h"],
        "x": ["y", "z"],
        "y": ["z"],
        "z": ["1", "2", "3"],
    }
    terms = ["a", "y"]
    mapped = tree_search.tree_map(mapping, terms)
    assert isinstance(mapped, dict)


def test_tree_mapping():
    annot, mapping = annotations.load_annotations(
        species="human",
        aspect="all",
        release="current",
        fetch=False,
    )
    go_subset = ["GO:0005737", "GO:0005654", "GO:0005634", "GO:0009986"]
    result = tree_search.tree_search(mapping, go_subset, annot)
    assert len(result) > len(annot)
