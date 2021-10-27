"""Test that setuptools-scm is working correctly"""
import gopher

def test_version():
    """Check that the version is not None"""
    assert gopher.__version__ is not None
