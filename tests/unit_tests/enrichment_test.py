"""Test that the enrichment functions are working correctly"""
from gopher import enrichment, mannwhitneyu
import pandas as pd
import random
from scipy import stats
import numpy as np


def test_entire_enrichment_analysis():
    """Check that the test enrichment function returns a dataframe.
    This test will run the full GO enrichment on a dataset and takes
    roughly 3.5 minutes to run."""
    df = pd.read_csv("../data/df.csv")
    df.set_index("Protein", inplace=True)
    result = enrichment.test_enrichment(df)
    assert type(result) is pd.DataFrame


def test_mannwhitneyu_small():
    """Test the Mann-Whitney U numba function returns the same values
    as the scipy function on a small dataset for a two-sided test"""
    list1 = []
    list2 = []
    random.seed(1)
    for i in range(0, 20):
        temp1 = []
        temp2 = []
        for i in range(0, 5):
            n = random.randint(1, 10)
            temp1.append(n)
            n = random.randint(5, 15)
            temp2.append(n)
        list1.append(temp1)
        list2.append(temp2)
    res_scipy = stats.mannwhitneyu(list1, list2)
    res_numba = mannwhitneyu.numba_mannwhitneyu(list1, list2)
    sp = np.round(res_scipy[1], 5)
    num = np.round(res_numba[1], 5)
    assert np.array_equal(sp, num)


def test_mannwhitneyu_greater():
    """Test the Mann-Whitney U numba function returns the same values
    as the scipy function on a small dataset for a one-sided, greater test"""
    list1 = []
    list2 = []
    random.seed(1)
    for i in range(0, 20):
        temp1 = []
        temp2 = []
        for i in range(0, 5):
            n = random.randint(1, 10)
            temp1.append(n)
            n = random.randint(5, 15)
            temp2.append(n)
        list1.append(temp1)
        list2.append(temp2)
    res_scipy = stats.mannwhitneyu(list1, list2, alternative="greater")
    res_numba = mannwhitneyu.numba_mannwhitneyu(
        list1, list2, alternative="greater"
    )
    sp = np.round(res_scipy[1], 5)
    num = np.round(res_numba[1], 5)
    assert np.array_equal(sp, num)


def test_mannwhitneyu_less():
    """Test the Mann-Whitney U numba function returns the same values
    as the scipy function on a small dataset for a one-sided, lesser test"""
    list1 = []
    list2 = []
    random.seed(1)
    for i in range(0, 20):
        temp1 = []
        temp2 = []
        for i in range(0, 5):
            n = random.randint(1, 10)
            temp1.append(n)
            n = random.randint(5, 15)
            temp2.append(n)
        list1.append(temp1)
        list2.append(temp2)
    res_scipy = stats.mannwhitneyu(list1, list2, alternative="less")
    res_numba = mannwhitneyu.numba_mannwhitneyu(
        list1, list2, alternative="less"
    )
    sp = np.round(res_scipy[1], 5)
    num = np.round(res_numba[1], 5)
    assert np.array_equal(sp, num)


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
    res_numba = mannwhitneyu.numba_mannwhitneyu(list1, list2)
    sp = np.round(res_scipy[1], 10)
    num = np.round(res_numba[1], 10)
    assert np.array_equal(sp, num)
