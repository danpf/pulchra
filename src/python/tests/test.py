#!/usr/bin/env python

from pathlib import Path
from subprocess import check_output
import tempfile
import unittest

from pypulchra import run_from_pdb_str

PYPULCHRA_TEST_DIR = Path(__file__).parent

class TestPyPulchra(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_basic_functions(self):
        pdb_fn = PYPULCHRA_TEST_DIR / "2ily_CA.pdb"
        with open(pdb_fn, encoding="utf8") as fh:
            input_pdb_text = fh.read()
        self.assertTrue(" N " not in input_pdb_text)
        ret = run_from_pdb_str(input_pdb_text)
        self.assertTrue(" N " in ret.pdb_str)

        out_pdb_fn = Path(self.temp_dir.name) / "pypulchra_output.pdb"
        check_output(["pypulchra", "-i", str(pdb_fn), "-o", str(out_pdb_fn)])
        self.assertTrue(out_pdb_fn.is_file())
