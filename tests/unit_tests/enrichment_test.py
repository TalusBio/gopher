"""Test that the enrichment functions are working correctly."""

import random

import numpy as np
import pandas as pd
from scipy import stats

import gopher
from gopher import annotations, enrichment


def test_entire_enrichment_analysis(generate_proteins):
    """Check that the test enrichment function returns a dataframe.

    This test runs the full GO enrichment on a dataset.
    """
    df = generate_proteins
    df.set_index("Protein", inplace=True)
    result = enrichment.test_enrichment(df)
    assert isinstance(result, pd.DataFrame)


def test_subset_enrichment_analysis(generate_proteins):
    """Check that the test enrichment function returns a dataframe.

    This test runs the GO enrichment on a subset of terms.
    """
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


def test_custom_annotations(generate_proteins):
    """Test custom annotations change results and include new terms."""
    prot = generate_proteins
    custom = annotations.generate_annotations(
        proteins=prot["Protein"][: len(prot) // 2], aspect="C", go_name="test"
    )
    full, _ = annotations.load_annotations(species="human")
    annot = pd.concat([custom, full])
    prot.set_index("Protein", inplace=True)
    terms = [
        "nucleus",
        "nucleoplasm",
        "euchromatin",
        "heterochromatin",
        "protein-DNA complex",
        "lysosome",
        "cytoplasm",
        "test",
    ]
    custom_result = enrichment.test_enrichment(
        prot, annotations=annot, go_subset=terms, aggregate_terms=False
    )
    normal_result = enrichment.test_enrichment(
        prot, go_subset=terms, aggregate_terms=False
    )
    assert not custom_result.equals(normal_result)
    terms_found = custom_result["GO Name"].unique()
    assert all(t in terms for t in terms_found)
    assert "test" in terms_found


def test_custom_mapping_annotations(
    generate_fake_proteins, generate_annotations, generate_mapping
):
    """Test the custom mapping and custom annotations on fake data."""
    subset = ["b", "i", "z"]
    result = enrichment.test_enrichment(
        generate_fake_proteins,
        annotations=generate_annotations,
        mapping=generate_mapping,
        go_subset=subset,
    )
    assert result["GO ID"].values.tolist() == subset


def test_mannwhitneyu_small(generate_arrays):
    """Test Mann-Whitney U numba matches SciPy on a small two-sided test."""
    list1, list2 = generate_arrays
    res_scipy = stats.mannwhitneyu(list1, list2)
    res_numba = gopher.stats.mannwhitneyu(list1, list2)
    sp = np.round(res_scipy[1], 5)
    num = np.round(res_numba[1], 5)
    np.testing.assert_allclose(sp, num)


def test_mannwhitneyu_greater(generate_arrays):
    """Test Mann-Whitney U numba matches SciPy on a greater one-sided test."""
    list1, list2 = generate_arrays
    res_scipy = stats.mannwhitneyu(list1, list2, alternative="greater")
    res_numba = gopher.stats.mannwhitneyu(list1, list2, alternative="greater")
    sp = np.round(res_scipy[1], 5)
    num = np.round(res_numba[1], 5)
    np.testing.assert_allclose(sp, num)


def test_mannwhitneyu_less(generate_arrays):
    """Test Mann-Whitney U numba matches SciPy on a lesser one-sided test."""
    list1, list2 = generate_arrays
    res_scipy = stats.mannwhitneyu(list1, list2, alternative="less")
    res_numba = gopher.stats.mannwhitneyu(list1, list2, alternative="less")
    sp = np.round(res_scipy[1], 5)
    num = np.round(res_numba[1], 5)
    np.testing.assert_allclose(sp, num)


def test_mannwhitneyu_big():
    """Test Mann-Whitney U numba matches SciPy on a larger two-sided test."""
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
