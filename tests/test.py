
from pathlib import Path

from pypulchra import run_from_pdb_str


def main():
    thisdir = Path(__file__).parent.resolve()
    with open(Path(thisdir, "2ily_CA.pdb")) as fh:
        input_pdb_text = fh.read()
    assert " N " not in input_pdb_text
    ret = run_from_pdb_str(input_pdb_text)
    assert " N " in ret.pdb_str


if __name__ == "__main__":
    main()
