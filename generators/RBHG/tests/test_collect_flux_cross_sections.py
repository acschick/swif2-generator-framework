import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "collect_flux_cross_sections.py"
SPEC = importlib.util.spec_from_file_location("collect_flux_cross_sections", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


SAMPLE = """
 Command-line args: seed= 1349198110 job= 0 nevents= 50000
 Runtime energies: E_lo= 11.338766000000000 E_hi= 11.338865999999999 E_coherent= 11.338816000000000
 Configured x range: 0.10000000149011612 to 0.89999997615814209
 *** CROSS_MAX FOUND: 628143371.56418991 at theta1= 1.5 theta2= 1.53
 test events 100000 Egamma 11.338816 total cross section nb 27930.337590069026 percent change 100.0
 test events 1000000 Egamma 11.338816 total cross section nb 28135.390404389374 percent change 0.72880742500149620
 test events 10000000 Egamma 11.338816 total cross section nb 28773.290757766426 percent change 2.2169878264788716
 test events 100000000 Egamma 11.338816 total cross section nb 28357.483593286484 percent change -1.4663048754381811
 Starting event generation...
"""


class ParseFluxCrossSectionTest(unittest.TestCase):
    def test_parses_final_integrated_cross_section(self):
        with tempfile.TemporaryDirectory() as directory:
            path = (
                Path(directory)
                / "CPP_FluxStudy3"
                / "generation"
                / "FLUX_TAGH_FFN_2205_135DEG_DBLRAD_CPP_ee_Heitler_v113"
                / "stdout_workflow_TAGH_001.out"
            )
            path.parent.mkdir(parents=True)
            path.write_text(SAMPLE)
            row = MODULE.parse_stdout(path)

        self.assertEqual(row["study"], "CPP_FluxStudy3")
        self.assertEqual(row["run_period"], "2205")
        self.assertEqual(row["polarization"], "135DEG")
        self.assertEqual(row["detector"], "TAGH")
        self.assertEqual(row["counter"], 1)
        self.assertEqual(row["job_index"], 0)
        self.assertEqual(row["requested_events"], 50000)
        self.assertAlmostEqual(row["E_coherent_GeV"], 11.338816)
        self.assertAlmostEqual(row["cross_section_nb"], 28357.483593286484)
        self.assertEqual(row["integration_test_events"], 100_000_000)
        self.assertAlmostEqual(row["final_percent_change"], -1.4663048754381811)
        self.assertEqual(row["integration_estimate_count"], 4)
        self.assertEqual(row["status"], "cross_section_complete")

    def test_marks_partial_integration(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "stdout_workflow_TAGM_007.out"
            path.write_text("\n".join(SAMPLE.splitlines()[:-2]))
            row = MODULE.parse_stdout(path)

        self.assertEqual(row["detector"], "TAGM")
        self.assertEqual(row["counter"], 7)
        self.assertEqual(row["status"], "partial_integration")


if __name__ == "__main__":
    unittest.main()
