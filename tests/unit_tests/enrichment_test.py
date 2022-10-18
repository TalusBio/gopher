"""Test that the enrichment functions are working correctly"""
from gopher import enrichment
import gopher
import pandas as pd
import random
from scipy import stats
import numpy as np


def test_entire_enrichment_analysis(generate_proteins):
    """Check that the test enrichment function returns a dataframe.
    This test will run the full GO enrichment on a dataset."""
    df = generate_proteins
    df.set_index("Protein", inplace=True)
    result = enrichment.test_enrichment(df)
    assert isinstance(result, pd.DataFrame)


def test_subset_enrichment_analysis(generate_proteins):
    """Check that the test enrichment function returns a dataframe.
    This test will run the GO enrichment on a subset of terms."""
    df = generate_proteins
    df.set_index("Protein", inplace=True)
    terms = [
        "nucleus",
        "nucleoplasm",
        "euchromatin",
        "heterochromatin",
        "protein-DNA complex",
        "lysosome",
        "cytoplasm",
    ]
    result = enrichment.test_enrichment(df, go_subset=terms)
    terms_found = result["GO Name"].unique()
    assert all(t in terms for t in terms_found)


def test_mannwhitneyu_small(generate_arrays):
    """Test the Mann-Whitney U numba function returns the same values
    as the scipy function on a small dataset for a two-sided test"""
    list1, list2 = generate_arrays
    res_scipy = stats.mannwhitneyu(list1, list2)
    res_numba = gopher.stats.mannwhitneyu(list1, list2)
    sp = np.round(res_scipy[1], 5)
    num = np.round(res_numba[1], 5)
    np.testing.assert_allclose(sp, num)


def test_mannwhitneyu_greater(generate_arrays):
    """Test the Mann-Whitney U numba function returns the same values
    as the scipy function on a small dataset for a one-sided, greater test"""
    list1, list2 = generate_arrays
    res_scipy = stats.mannwhitneyu(list1, list2, alternative="greater")
    res_numba = gopher.stats.mannwhitneyu(list1, list2, alternative="greater")
    sp = np.round(res_scipy[1], 5)
    num = np.round(res_numba[1], 5)
    np.testing.assert_allclose(sp, num)


def test_mannwhitneyu_less(generate_arrays):
    """Test the Mann-Whitney U numba function returns the same values
    as the scipy function on a small dataset for a one-sided, lesser test"""
    list1, list2 = generate_arrays
    res_scipy = stats.mannwhitneyu(list1, list2, alternative="less")
    res_numba = gopher.stats.mannwhitneyu(list1, list2, alternative="less")
    sp = np.round(res_scipy[1], 5)
    num = np.round(res_numba[1], 5)
    np.testing.assert_allclose(sp, num)


def test_mannwhitneyu_big():
    """Test the Mann-Whitney U numba function returns the same values
    as the scipy function on a larger dataset for a two-sided test"""
    list1 = []
    list2 = []
    random.seed(10)
    for i in range(0, 10):
        temp1 = []
        temp2 = []
        for i in range(0, 100):
            n = random.randint(1, 100)
            temp1.append(n)
            n = random.randint(50, 150)
            temp2.append(n)
        list1.append(temp1)
        list2.append(temp2)
    res_scipy = stats.mannwhitneyu(list1, list2)
    res_numba = gopher.stats.mannwhitneyu(list1, list2)
    sp = np.round(res_scipy[1], 10)
    num = np.round(res_numba[1], 10)
    np.testing.assert_allclose(sp, num)
