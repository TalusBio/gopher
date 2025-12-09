"""Wrapper helpers for GO graph traversal."""

from . import tree_search


def new_map(mapping, go_subset):
    """Return expanded mapping for the subset of GO terms."""
    return tree_search.new_tree_map(mapping, go_subset)


def graph_search(mapping, go_subset, annot):
    """Expand annotations using GO child relationships."""
    return tree_search.tree_search(mapping, go_subset, annot)
