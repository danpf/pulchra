import sys
from typing import List
import argparse
import glob
from pathlib import Path

print(__file__)
print(list(sorted(glob.glob(str(Path(__file__).parent / "*" )))))
print(list(sorted(glob.glob(str(Path(__file__).parent.parent / "*" )))))

from ._pypulchra import pulchra as pulchra, run_from_pdb_str as run_from_pdb_str


def parseargs(raw_args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pulchra protein pdb/cif CA -> fullatom algorithm.")
    parser.add_argument("-i", "--input-pdb", required=True, help="The pdb file you want as input")
    parser.add_argument("-o", "--output-pdb", required=True, help="The pdb file you want as output")
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
