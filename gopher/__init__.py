"""See the README for detailed documentation and examples."""
try:
    from importlib.metadata import PackageNotFoundError, version

    try:
        __version__ = version(__name__)
    except PackageNotFoundError:
        pass

except ImportError:
    from pkg_resources import DistributionNotFound, get_distribution

    try:
        __version__ = get_distribution(__name__).version
    except DistributionNotFound:
        pass

from .annotations import generate_annotations, load_annotations
from .config import get_data_dir, set_data_dir
from .display_data import (
    get_annotations,
    get_rankings,
    in_term,
    map_proteins,
    roc,
)
from .enrichment import test_enrichment
from .normalize import normalize_values
from .parsers import read_encyclopedia, read_metamorpheus, read_diann
