import pytest
import logging
from pathlib import Path
import os
import pandas as pd
from b3p.cli.app_state import AppState
from b3p.cli.build_app import BuildApp
from b3p.cli.ccblade_app import CCBladeApp, has_ccblade

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)


@pytest.fixture(scope="session")
def run_ccblade(temp_example_dir):
    """Fixture to build the blade and run the CCBlade process."""
    if not has_ccblade:
        pytest.skip("CCBlade functionality is not available")

    original_dir = os.getcwd()
    os.chdir(temp_example_dir)
    try:
        state = AppState()
        yml_path = Path("blade_test.yml")

        logger.info("Building blade with YAML: %s", yml_path)
        build_app = BuildApp(state, yml_path)
        build_app.build()

        logger.info("Running CCBlade analysis")
        ccblade_app = CCBladeApp(state, yml_path)
        ccblade_app.ccblade()

        workdir = temp_example_dir / "temp_blade"
        logger.info("CCBlade workdir: %s", workdir)
        yield {
            "workdir": workdir,
            "temp_dir": temp_example_dir.parent,
        }
    finally:
        os.chdir(original_dir)


def test_ccblade_output(run_ccblade):
    """Test if CCBlade generates the expected output table matching the reference."""
    workdir = run_ccblade["workdir"]
    temp_dir = run_ccblade["temp_dir"]

    generated_file = workdir / "ccblade_output.csv"
    logger.info("Checking for generated CCBlade output: %s", generated_file)
    assert (
        generated_file.exists()
    ), f"CCBlade should generate {generated_file} in the working directory"

    reference_file = temp_dir / "data" / "ccblade_output.csv"
    logger.info("Reference CCBlade output: %s", reference_file)
    assert (
        reference_file.exists()
    ), f"Reference file {reference_file} not found in temp data directory"

    generated_df = pd.read_csv(generated_file)
    reference_df = pd.read_csv(reference_file)

    assert generated_df.equals(
        reference_df
    ), "Generated CCBlade output does not match reference file"
