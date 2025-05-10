import logging
from pathlib import Path

logger = logging.getLogger(__name__)

try:
    from b3p.bem import ccblade_run

    has_ccblade = True
except ImportError:
    logger.warning("Could not import ccblade_run. Functionality will be disabled.")
    has_ccblade = False


class CCBladeApp:
    def __init__(self, state, yml: Path):
        self.state = state
        self.yml = yml

    def ccblade(self):
        if has_ccblade:
            ccblade = ccblade_run.ccblade_run(self.yml)
            ccblade.run()
        else:
            logger.error("ccblade_run is not available.")
