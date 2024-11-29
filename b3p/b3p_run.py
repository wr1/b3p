
import json
import os 

class b3p_run:
    def __init__(self, config_file: str, state_file: str = "b3prun_state.json"):
        self.config_file = config_file
        self.state_file = state_file
        self.state = {
            "geometry": False,
            "mesh": False,
            "laminates": False,
            "output_files": {},
        }

    def save_state(self):
        with open(self.state_file, "w") as f:
            json.dump(self.state, f)
        print(f"State saved to {self.state_file}.")

    def load_state(self):
        if os.path.exists(self.state_file):
            with open(self.state_file, "r") as f:
                self.state = json.load(f)
            print(f"State loaded from {self.state_file}.")
        else:
            print(f"No state file found. Starting fresh.")

    def check_prerequisites(self, step: str):
        prerequisites = {
            "geometry": [],
            "mesh": ["geometry"],
            "laminates": ["mesh"],
        }
        missing_steps = [
            prereq for prereq in prerequisites.get(step, [])
            if not self.state.get(prereq, False)
        ]
        if missing_steps:
            raise RuntimeError(
                f"Missing prerequisite steps for '{step}': {', '.join(missing_steps)}"
            )

    def build_geometry(self):
        print("Building geometry...")
        # Simulate output
        self.state["geometry"] = True
        self.state["output_files"]["geometry"] = "geometry_output.json"
        self.save_state()

    def build_mesh(self):
        self.check_prerequisites("mesh")
        print("Building mesh...")
        # Simulate output
        self.state["mesh"] = True
        self.state["output_files"]["mesh"] = "mesh_output.json"
        self.save_state()

    def apply_laminates(self):
        self.check_prerequisites("laminates")
        print("Applying laminates...")
        # Simulate output
        self.state["laminates"] = True
        self.state["output_files"]["laminates"] = "laminates_output.json"
        self.save_state()

    def run_step(self, step: str):
        steps = {
            "geometry": self.build_geometry,
            "mesh": self.build_mesh,
            "laminates": self.apply_laminates,
        }
        if step not in steps:
            raise ValueError(f"Unknown step: {step}")
        steps[step]()
