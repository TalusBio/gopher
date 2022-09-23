"""Numba Mann-Whitney U test"""
import numpy as np
from scipy import stats
import numba as nb


@nb.njit(parallel=True)
def tiecorrect(rankvals):
    """Parallelized version of scipy.stats.tiecorrect"""
    tc = np.ones(rankvals.shape[1], dtype=np.float64)
    for j in nb.prange(rankvals.shape[1]):
        arr = np.sort(np.ravel(rankvals[:, j]))
        is_not_tie = np.concatenate(
            (np.array([True]), arr[1:] != arr[:-1], np.array([True]))
        )
        idx = np.nonzero(is_not_tie)[0]
        cnt = np.diff(idx).astype(np.float64)

        size = np.float64(arr.size)
        if size >= 2:
            tc[j] = 1.0 - (cnt**3 - cnt).sum() / (size**3 - size)

    return tc


@nb.njit(parallel=True)
def rankdata(data):
    """Parallelized version of scipy.stats.rankdata"""
    ranked = np.empty(data.shape, dtype=np.float64)
    for j in nb.prange(data.shape[1]):
        arr = np.ravel(data[:, j])
        sorter = np.argsort(arr)

        arr = arr[sorter]
        obs = np.concatenate((np.array([True]), arr[1:] != arr[:-1]))

        dense = np.empty(obs.size, dtype=np.int64)
        dense[sorter] = obs.cumsum()

        # cumulative counts of each unique value
        count = np.concatenate((np.nonzero(obs)[0], np.array([len(obs)])))
        ranked[:, j] = 0.5 * (count[dense] + count[dense - 1] + 1)

    return ranked


def mannwhitneyu(x, y, alternative="two-sided", use_continuity=True):
    """Version of Mann-Whitney U-test that runs in parallel on 2d arrays.
    This is the asymptotic algo only."""
    x = np.asarray(x)
    y = np.asarray(y)
    assert x.shape[1] == y.shape[1]

    n1 = x.shape[0]
    n2 = y.shape[0]

    # rank the data for the test
    ranked = rankdata(np.concatenate((x, y)))
    rankx = ranked[0:n1, :]  # get the x-ranks
    r = np.sum(rankx, axis=0)
    u1 = n1 * n2 + (n1 * (n1 + 1)) / 2.0 - r  # calc U for x
    u2 = n1 * n2 - u1  # remainder is U for y

    # check for ties in the rankings
    T = tiecorrect(ranked)

    # if *everything* is identical we'll raise an error, not otherwise
    if np.all(T == 0):
        raise ValueError("All numbers are identical")

    # get mean and standard deviation
    sd = np.sqrt(T * n1 * n2 * (n1 + n2 + 1) / 12.0)
    meanrank = n1 * n2 / 2.0 + 0.5 * use_continuity

    # do a one-sided or two-sided test depending on user specification
    if alternative == "greater":
        U, f = u2, 1
    elif alternative == "less":
        U, f = u1, 1
    else:
        U, f = np.maximum(u1, u2), 2

    # calculate z and use it to find the p-values
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (U - meanrank) / sd

    p = stats.norm.sf(z)
    p *= f
    p = np.clip(p, 0, 1)

    auc = U / (n1 * n2)

    return U, p
