"""Utility functions"""
import socket
import gzip
from pathlib import Path

import requests


def http_download(url, path):
    """Download a file using GET.

    Parameters
    ----------
    url : str
        The URL of the file to download.
    path : Path
        The downloaded file path.
    """
    path = Path(path)
    with requests.get(url, stream=True) as res:
        res.raise_for_status()
        try:
            with path.open("wb") as out_ref:
                for chunk in res.iter_content(chunk_size=8192):
                    out_ref.write(chunk)

        except (
            ConnectionRefusedError,
            ConnectionResetError,
            socket.timeout,
            socket.gaierror,
            socket.herror,
            EOFError,
            OSError,
        ) as err:
            path.unlink()
            raise err


def decompress(gz_file, target_file):
    """Decompress a gzipped file.

    Parameters
    ----------
    gz_file : Path
        The path to the gzipped file.
    target_file : Path
        The path to the decompressed file.
    """
    with gzip.open(gz_file, "rb") as gz_ref:
        with target_file.open("wb") as out_ref:
            out_ref.write(gz_ref.read())

    return target_file