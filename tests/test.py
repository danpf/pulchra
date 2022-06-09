#!/usr/bin/env python
print("start?")

from pathlib import Path
from subprocess import check_output

from pypulchra import run_from_pdb_str

print("post import?")

def main():
    print("main?")
    thisdir = Path(__file__).parent.resolve()
    pdb_fn = Path(thisdir) / "2ily_CA.pdb"
    print("s", pdb_fn)
    # with open(pdb_fn) as fh:
    #     input_pdb_text = fh.read()
    # assert " N " not in input_pdb_text
    # ret = run_from_pdb_str(input_pdb_text)
    # assert " N " in ret.pdb_str
    # print("halfway")
    #
    # out_pdb_fn = Path("__out.pdb")
    # if out_pdb_fn.is_file():
    #     check_output(["rm", str(out_pdb_fn)])
    # print("2halfway")
    # assert not out_pdb_fn.is_file()
    # check_output(["pypulchra", "-i", str(pdb_fn), "-o", str(out_pdb_fn)])
    # print("3halfway")
    # assert out_pdb_fn.is_file()
    # check_output(["rm", str(out_pdb_fn)])
    # print("4halfway")


print("about to main import?")
if __name__ == "__main__":
    main()
