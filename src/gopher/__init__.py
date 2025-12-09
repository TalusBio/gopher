"""See the README for detailed documentation and examples."""

try:
    from importlib.metadata import PackageNotFoundError, version

    try:
        __version__ = version(__name__)
    except PackageNotFoundError:
        __version__ = None

except ImportError:
    from pkg_resources import DistributionNotFound, get_distribution

    try:
        __version__ = get_distribution(__name__).version
    except DistributionNotFound:
        __version__ = None

from . import graph_search
from .config import get_data_dir, set_data_dir
from .enrichment import test_enrichment
from .parsers import read_diann, read_encyclopedia, read_metamorpheus
from .version import _get_version

# Fall back to version helper if metadata lookup failed
if __version__ is None:
    __version__ = _get_version()
