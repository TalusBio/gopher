"""Fixtures used for testing."""

import random
import string

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def generate_proteins():
    """Generate a random list of protein data."""
    prot = [
        "P10809",
        "P35527",
        "Q9UMS4",
        "P35637",
        "P07437",
        "Q86UX7",
        "Q8N766",
        "P52907",
        "Q9NV31",
        "P50570",
        "Q9BZJ0",
        "Q9GZM5",
        "P08758",
        "Q9BRU9",
        "Q8NBJ7",
        "Q15532",
        "Q14241",
        "Q14562",
        "Q9NZM5",
        "P52926",
    ]
    sample1 = []
    sample2 = []
    sample3 = []
    for _ in range(len(prot)):
        sample1.append(random.randint(100000, 1000000))
        sample2.append(random.randint(100000, 1000000))
        sample3.append(random.randint(100000, 1000000))
    data = {
        "Protein": prot,
        "Sample 1": sample1,
        "Sample 2": sample2,
        "Sample 3": sample3,
    }
    return pd.DataFrame(data)


@pytest.fixture
def generate_arrays():
    """Generate arrays with random numbers."""
    rng = np.random.default_rng(1)  # A random number generator
    list1 = rng.integers(1, 10, size=(20, 5)).tolist()
    list2 = rng.integers(5, 15, size=(20, 5)).tolist()
    return list1, list2


@pytest.fixture
def generate_mapping():
    """Create mock mapping."""
    mapping = {
        "a": ["b", "c", "d"],
        "b": ["e"],
        "d": ["g", "h"],
        "e": ["f"],
        "i": ["j", "k", "y"],
        "x": ["y", "z"],
        "y": ["z"],
        "z": ["l", "m", "n"],
    }
    return mapping


@pytest.fixture
def generate_annotations():
    """Create annotations for mock data."""
    annot = {
        "uniprot_accession": list(range(0, 26)),
        "go_id": list(string.ascii_lowercase),
        "aspect": ["c"] * 26,
        "go_name": list(string.ascii_uppercase),
    }
    return pd.DataFrame.from_dict(annot)


@pytest.fixture
def generate_fake_proteins():
    """Generate a random list of protein data."""
    prot = list(range(0, 26))
    sample1 = []
    sample2 = []
    sample3 = []
    for _ in range(len(prot)):
        sample1.append(random.randint(1, 10))
        sample2.append(random.randint(1, 10))
        sample3.append(random.randint(1, 10))
    data = {
        "Protein": prot,
        "Sample 1": sample1,
        "Sample 2": sample2,
        "Sample 3": sample3,
    }
    return pd.DataFrame(data)
