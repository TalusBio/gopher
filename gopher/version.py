"""Get the package version"""


def _get_version():
    """Retreive the package version.

    Using setuptools-scm dynamically sets the pacakge version depending on
    the latest git release and commits since that point. This function
    returns the current version number.

    Returns
    -------
    str or None
        The package version number. If not version is found, returns None.
    """
    try:
        # Fast, but only works in Python 3.8+
        from importlib.metadata import PackageNotFoundError, version

        try:
            return version("gopher-enrich")
        except PackageNotFoundError:
            return None

    except ImportError:
        # Slow, but works for all Python 3+
        from pkg_resources import DistributionNotFound, get_distribution

        try:
            return get_distribution("gopher-enrich").version
        except DistributionNotFound:
            return None
