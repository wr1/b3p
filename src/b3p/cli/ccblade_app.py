import logging
from pathlib import Path

logger = logging.getLogger(__name__)

try:
    from b3p.models.config import BladeConfig  # Ensure models are used
    from b3p.bem.ccblade_run import ccblade_run  # Import the updated class

    has_ccblade = True

except ImportError:
    logger.warning("Could not import ccblade_run. Functionality will be disabled.")
    has_ccblade = False


from treeparse import cli, command, option


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


def ccblade_callback(yml: Path):
    state = AppState.get_instance()
    app = CCBladeApp(state, yml)
    app.ccblade()


ccblade_cli = cli(
    name="ccblade",
    help="Run CCBlade analysis",
    line_connect=True,
    show_types=True,
    show_defaults=True,
    options=[
        option(
            flags=["--yml", "-y"],
            arg_type=Path,
            required=True,
            help="Path to YAML config file",
        ),
    ],
)

ccblade_cli.commands.append(
    command(
        name="run",
        help="Run CCBlade",
        callback=ccblade_callback,
        arguments=[],
    )
)
