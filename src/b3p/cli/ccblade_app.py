import logging
from pathlib import Path
from b3p.bem.ccblade_run import ccblade_run  # Import the updated class

logger = logging.getLogger(__name__)

try:
    from b3p.models.config import BladeConfig  # Ensure models are used
    has_ccblade = True
except ImportError:
    logger.warning("Could not import ccblade_run. Functionality will be disabled.")
    has_ccblade = False


class CCBladeApp:
    def __init__(self, state, yml: Path):
        self.state = state  # AppState instance
        self.yml = yml

    def ccblade(self):
        if has_ccblade:
            # Use the loaded config and yml_dir from state
            ccblade = ccblade_run(self.state.config, self.state.yml_dir)
            ccblade.run()
        else:
            logger.error("ccblade_run is not available.")
