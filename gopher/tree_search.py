import copy
import pandas as pd


def tree_search(mapping, go_subset, annot):
    """Incorporates the tree search to get all children from parent node.

    First, get the GO IDs for the terms of interest. Run the tree algorithm to
    get the updated tree mapping. Update the tree to reflect the new mapping.

    Parameters
    ----------
    mapping : dict
        A dictionary with the mapping of the terms of interest as keys and children as values.
    go_subset : list
        List of terms of interest as their GO IDs.
    annot : pd.DataFrame
        Dataframe of annotations of each protein and the term(s) it is associated with.
    Returns
    -------
    pandas.DataFrame
        The dataframe with the new annotations of the mapping incorporated.
    """
    # Get the go ids for the terms of interest
    in_names = annot["go_name"].isin(go_subset)
    in_ids = annot["go_id"].isin(go_subset)
    subset = annot.loc[in_names | in_ids, :]["go_id"].unique().tolist()
    # Get the new mapping for the subset of terms with all the children associated with each term
    mapped = tree_map(mapping, subset)
    # Incorporate the new mapping into the graph and return the updated graph
    new_annot = update_tree(mapped, annot)
    return new_annot


def tree_map(mapping, subset):
    """Tree search algorithm that gets all children nodes (or terms) from the specified parent.

    Parameters
    ----------
    mapping : dict
        A dictionary with the mapping of the terms of interest as keys and
        children terms as values.
    subset : list
        List of terms of interest as their GO IDs.

    Returns
    -------
    dictionary
        Dictionary with terms in the subset as the keys and all children of the term as values
    """
    subset_mapping = {}
    # For every term, get all the children of that node
    for item in subset:
        result = map(item, mapping)
        if result:
            subset_mapping[item] = result
        else:
            subset_mapping[item] = []
    # Return the new mapping
    return subset_mapping


def map(term, mapping):
    """Tree search recursive helper function that gets all children nodes (or terms) from the specified parent.

    Parameters
    ----------
    term: str
        Current term
    mapping : dict
        A dictionary with the mapping of the terms of interest as keys and
        children terms as values.

    Returns
    -------
    list
        list of all terms that relate to the term of interest
    """
    # Base case: if there are no children of the current term, return
    if term not in mapping.keys():
        return
    # Get all children of the current term
    children = mapping[term]
    result = copy.copy(children)
    # Iterate through each child, get their children, and add them to the list of children
    for child in children:
        res = map(child, mapping)
        if res:
            result += res
    # Return the list of child nodes
    return result


def update_tree(mapping, annot):
    """Tree search algorithm that gets all children nodes (or terms) from the specified parent.

    Parameters
    ----------
    mapping : dict
        Dictionary with terms in the subset as the keys and all children of the term as values.
    annot : pd.DataFrame
        Dataframe of annotations of each protein and the term(s) it is associated with.
    Returns
    -------
    pandas.DataFrame
        The dataframe with the new annotations of the mapping incorporated.
    """
    # For every term in the new mapping, get the list of children and the go name
    for key in mapping.keys():
        list = mapping[key]
        name = annot.loc[annot["go_id"] == key].iloc[0]["go_name"]
        # For every child, get the annotations for that child. Change the id and name to match the parent.
        for value in list:
            new = annot.loc[annot["go_id"] == value]
            new["go_id"] = key
            new["go_name"] = name
            annot = pd.concat([annot, new], ignore_index=True)
    # Drop duplicate rows and return the updated annotation dataframe
    annot = annot.drop_duplicates()
    return annot
