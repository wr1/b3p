from pathlib import Path

try:
    from b3p.bem import ccblade_run

    has_ccblade = True
except ImportError:
    print("** Could not import ccblade_run. Functionality will be disabled.")
    has_ccblade = False


class CCBladeApp:
    def __init__(self, state):
        self.state = state

    def ccblade(self, yml: Path):
        if has_ccblade:
            ccblade = ccblade_run.ccblade_run(yml)
            ccblade.run()
        else:
            print("** ccblade_run is not available.")
