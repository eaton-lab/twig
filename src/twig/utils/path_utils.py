#!/usr/bin/env python

"""Simple utilities for file paths.

"""

from pathlib import Path



def expand_multiple_paths(paths: list[Path]) -> list[Path]:
    """Return a list of Paths expanded from >=1 paths using regex.
    
    Returns full paths with user expanded.
    """
    fpaths = []
    for path in paths:
        if path.is_file():
            fpaths.append(path.expanduser().absolute())
        else:
            for path in path.parent.glob(path.name):
                if path.exists():
                    fpaths.append(path.expanduser().absolute())
    return fpaths


if __name__ == "__main__":
    pass
