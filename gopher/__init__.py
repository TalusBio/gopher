"""See the README for detailed documentation and examples."""
try:
    from importlib.metadata import version, PackageNotFoundError

    try:
        __version__ = version(__name__)
    except PackageNotFoundError:
        pass

except ImportError:
    from pkg_resources import get_distribution, DistributionNotFound

    try:
        __version__ = get_distribution(__name__).version
    except DistributionNotFound:
        pass

from .config import get_data_dir, set_data_dir
from .parsers import read_encyclopedia, read_metamorpheus
from .enrichment import test_enrichment
from .display_data import map_proteins
