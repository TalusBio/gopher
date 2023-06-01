"""Utility functions"""
import socket
from pathlib import Path

import requests


def http_download(url: str, path: Path):
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
