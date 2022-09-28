import copy
import pandas as pd


def tree_search(mapping, go_subset, annot):
    in_names = annot["go_name"].isin(go_subset)
    in_ids = annot["go_id"].isin(go_subset)
    subset = annot.loc[in_names | in_ids, :]["go_id"].unique().tolist()
    mapped = tree_map(mapping, subset)
    result = reflect_mapping(mapped, annot)
    return result


def tree_map(mapping, subset):
    subset_mapping = {}
    for item in subset:
        result = map(item, mapping)
        if result:
            subset_mapping[item] = result
        else:
            subset_mapping[item] = []
    return subset_mapping


def map(term, mapping):
    if term not in mapping.keys():
        return
    list = mapping[term]
    result = copy.copy(list)
    for item in list:
        x = map(item, mapping)
        if x:
            result += x
    return result


def reflect_mapping(mapping, annot):
    for key in mapping.keys():
        list = mapping[key]
        name = annot.loc[annot["go_id"] == key].iloc[0]["go_name"]
        for value in list:
            new = annot.loc[annot["go_id"] == value]
            new["go_id"] = key
            new["go_name"] = name
            annot = pd.concat([annot, new], ignore_index=True)
    annot = annot.drop_duplicates()
    return annot
