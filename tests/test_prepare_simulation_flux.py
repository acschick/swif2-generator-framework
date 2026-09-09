import importlib.util
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "hddm_scripts"))
sys.modules.setdefault("hddm_s", types.SimpleNamespace())

SPEC = importlib.util.spec_from_file_location(
    "prepareSimulation", ROOT / "prepareSimulation.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class FluxSimulationPreparationTest(unittest.TestCase):
    def make_config(self, vectors_dir):
        return {
            "submission_status": {
                "workflow_name": "gen_CPP_FluxStudy_FLUX_FFN_2205_135DEG_TAGH_x"
            },
            "rbhg_config": {
                "generation_info": {"study_name": "CPP_FluxStudy", "nametag": "FLUX_TAGH"},
                "physics_settings": {
                    "run_period": "2205",
                    "polarization_deg": "135",
                },
                "directory_paths": {
                    "rbhg_generation": {
                        "vectors_hddm_directory": str(vectors_dir),
                    }
                },
                "mcwrapper_settings": {
                    "characteristic_run": 101582,
                    "batch_system": "swif2",
                },
                "flux_study": {
                    "enabled": True,
                    "detector": "TAGH",
                    "counter_count": 2,
                    "events_per_counter": 50000,
                },
                "job_parameters": [
                    {"job_i": 0, "counter": 1, "detector": "TAGH", "job_label": "TAGH_001", "nevents": "50000"},
                    {"job_i": 1, "counter": 2, "detector": "TAGH", "job_label": "TAGH_002", "nevents": "50000"},
                ],
            },
        }

    def test_prepares_one_canonical_run_request_per_counter(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            vectors = root / "vectors"
            vectors.mkdir()
            (vectors / "vectors_0.hddm").touch()
            (vectors / "vectors_1.hddm").touch()
            config = self.make_config(vectors)
            config_path = root / "rbhg_config.json"
            config_path.write_text(json.dumps(config))

            with mock.patch("hddm_event_utils.count_hddm_events", return_value=50000), \
                 mock.patch.object(MODULE, "create_mcwrapper_config", return_value=True) as create, \
                 mock.patch.object(MODULE, "update_completion_status"):
                result = MODULE.process_flux_individual_directory(
                    config, str(config_path), str(vectors)
                )

            self.assertTrue(result["success"])
            self.assertEqual(result["event_count"], 100000)
            self.assertEqual(create.call_count, 2)
            for call in create.call_args_list:
                self.assertEqual(call.kwargs["run_range_override"], 101582)
                self.assertEqual(call.kwargs["events_per_file_override"], 50000)
                self.assertEqual(call.kwargs["batch_mode_override"], 1)

            master = Path(result["submission_script"]).read_text()
            self.assertIn("TAGH_001/gluex_MC_sub.csh", master)
            self.assertIn("TAGH_002/gluex_MC_sub.csh", master)
            self.assertIn(
                "swif2 run SIM_CPP_FluxStudy_FLUX_FFN_2205_135DEG_TAGH_x",
                master,
            )

    def test_counter_config_uses_one_file_one_run_and_one_tag(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            hddm_file = root / "vectors_0.hddm"
            hddm_file.write_bytes(b"test")
            config = self.make_config(root)
            rbhg = config["rbhg_config"]
            rbhg["physics_settings"].update({
                "experiment": "CPP",
                "lepton_type": "ee",
                "target": "pb208",
                "form_factor": "FFN",
            })
            rbhg["directory_paths"].update({
                "base_paths": {
                    "framework_home": str(ROOT),
                    "farm_out_base": str(root / "farm_out"),
                },
                "templates": {
                    "mcwrapper_template_directory": str(ROOT / "template_MCWrapper")
                },
                "rbhg_generation": {
                    "vectors_hddm_directory": str(root),
                    "output_directory": str(root / "generation"),
                },
            })
            rbhg["mcwrapper_settings"].update({
                "recon_env": "/tmp/recon.xml",
                "analysis_env": "/tmp/analysis.xml",
                "sim_env": "/tmp/sim.xml",
                "background_version": "NONE",
                "jana_config": str(
                    ROOT / "template_MCWrapper/jana_configs/jana_CPPlaunch1_epem.config"
                ),
                "variation": "mc",
            })
            output = root / "MCWrapper" / "TAGH_001"
            output.mkdir(parents=True)

            success = MODULE.create_mcwrapper_config(
                config,
                str(output),
                str(hddm_file),
                50000,
                counter_label="TAGH_001",
                batch_mode_override=1,
                events_per_file_override=50000,
                run_range_override=101582,
            )

            self.assertTrue(success)
            mc_config = (output / "MC.config").read_text()
            submit = (output / "gluex_MC_sub.csh").read_text()
            self.assertIn(f"GENERATOR=file:{hddm_file}", mc_config)
            self.assertIn("TAG=TAGH_001", mc_config)
            self.assertIn("/Simulation/TAGH_001", mc_config)
            self.assertNotIn("RBHGSPATHTOHDDMFILE", mc_config)
            self.assertIn("101582 50000 per_file=50000 batch=1", submit)

    def test_event_count_mismatch_blocks_master_submission(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            vectors = root / "vectors"
            vectors.mkdir()
            (vectors / "vectors_0.hddm").touch()
            (vectors / "vectors_1.hddm").touch()
            config = self.make_config(vectors)
            config_path = root / "rbhg_config.json"
            config_path.write_text(json.dumps(config))

            with mock.patch("hddm_event_utils.count_hddm_events", side_effect=[50000, 49999]), \
                 mock.patch.object(MODULE, "create_mcwrapper_config", return_value=True), \
                 mock.patch.object(MODULE, "update_completion_status"):
                with self.assertRaisesRegex(RuntimeError, "TAGH_002"):
                    MODULE.process_flux_individual_directory(
                        config, str(config_path), str(vectors)
                    )

            blocked_script = root / "MCWrapper" / "submit_all_flux_counters.sh"
            self.assertTrue(blocked_script.exists())
            self.assertIn("preparation is incomplete", blocked_script.read_text())

    def test_existing_master_detector_pattern_finds_actual_flux_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            child = root / "2205_135DEG_DBLRAD_CPP_ee_Heitler_v113"
            child.mkdir()
            individual = child / "rbhg_config.json"
            individual.write_text("{}")
            master_path = root / "run_period_master_config.json"
            master = {
                "run_period_master_config": {
                    "study_info": {
                        "dataset_name": "FLUX_TAGH",
                        "preset": "2205",
                    },
                    "batch_processing": {
                        "individual_configs": [
                            "2205_135DEG_TAGH_*/rbhg_config.json"
                        ]
                    },
                }
            }

            with mock.patch.object(MODULE, "update_completion_status"), \
                 mock.patch.object(MODULE, "load_and_validate_config", return_value=({}, "individual", 1)), \
                 mock.patch.object(MODULE, "process_individual_directory", return_value={"success": True, "event_count": 50000, "merged_hddm_file": None}):
                result = MODULE.process_run_period_master(
                    master, str(master_path)
                )

            self.assertTrue(result["success"])
            self.assertEqual(len(result["processed_configs"]), 1)
            self.assertEqual(result["processed_configs"][0]["config_path"], str(individual))


if __name__ == "__main__":
    unittest.main()
