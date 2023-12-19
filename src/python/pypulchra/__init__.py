from __future__ import annotations

import argparse
import sys

from ._pypulchra import pulchra as pulchra
from ._pypulchra import run_from_pdb_str as run_from_pdb_str

try:
    from ._version import version as __version__
    from ._version import version_tuple as version_tuple
except ImportError:
    __version__ = "unknown version"
    version_tuple = (0, 0, "unknown version")


def parseargs(raw_args: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pulchra protein pdb/cif CA -> fullatom algorithm."
    )
    parser.add_argument(
        "-i", "--input-pdb", required=True, help="The pdb file you want as input"
    )
    parser.add_argument(
        "-o", "--output-pdb", required=True, help="The pdb file you want as output"
    )
    return parser.parse_args(raw_args)


def commandline() -> None:
    args = parseargs(sys.argv[1:])
    with open(args.input_pdb) as fh:
        pdb_text = fh.read()
    ret = run_from_pdb_str(pdb_text)
    if not ret.pdb_str:
        raise RuntimeError("Failure to generate any pdb text")
    with open(args.output_pdb, "w") as fh:
        fh.write(ret.pdb_str)
