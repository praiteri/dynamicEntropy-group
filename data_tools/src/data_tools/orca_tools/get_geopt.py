import logging

from ..orca_tools.extract_from_orca_output import extract_from_orca_output
from ..logger import formatting


class get_geopt(extract_from_orca_output):
    """
    Extract the geometry optimisation log
    """

    def read(self):
        idx = self.search_string("Geometry Optimization Run")
        # Not a geometry optimisation run
        if idx is None:
            self.result = None
            return

        energies_indices = self.search_string("FINAL SINGLE POINT ENERGY")
        energies = self.parse_single_line(energies_indices)

        self.read_blocks("|Geometry convergence|", [""], 0)
        convergence = self.result
        try:
            ncyles = len(convergence)
        except:
            ncyles = 0

        converged = self.search_string("THE OPTIMIZATION HAS CONVERGED") is not None

        self.result = {
            "cycles": ncyles,
            "converged": converged,
            "energies": energies,
            "convergence": convergence,
        }
        return

    def write(self, **kwargs):
        if self.result is None:
            return
        self.logger.info(formatting().dashes)
        if self.result["converged"]:
            self.logger.info(self.fmt.format("Geometry optimisation summary", "CONVERGED"))
        else:
            self.logger.info(self.fmt.format("Geometry optimisation summary", "NOT CONVERGED"))
        self.logger.info(self.fmt.format("Number of cycles", self.result["cycles"]))

        for i in range(len(self.result["energies"]) - 1):
            self.logger.verbose(
                self.fmt.format(f"Step {i} energy (Eh)", self.result["energies"][i].split()[-1])
            )
        # Final energy
        self.logger.info(
            self.fmt.format("Final energy (Eh)", self.result["energies"][-1].split()[-1])
        )

        # Convergence data for each step
        if self.logger.level <= logging.DEBUG:
            for i in range(len(self.result["convergence"])):
                for ll in self.result["convergence"][i]:
                    self.logger.debug(">>> " + ll)
        else:
            # Final convergence data
            n = len(self.result["convergence"][-1])
            for i in range(2, n):
                ll = self.result["convergence"][-1][i]
                self.logger.verbose(">>> " + ll)

        # self.logger.info(formatting().dashes)
        return

    def show_xyz(self):
        if self.result is None:
            return
        self.logger.info("Final coordinates")
        xyz = self.result["coordinates"][-1]
        self.logger.info(f"{len(xyz)}")
        self.logger.info("")
        for ll in xyz:
            self.logger.info(f"{ll}")
        return

    def write_xyz(self, filname="final_coordinates.xyz"):
        if self.result is None:
            return
        xyz = self.result["coordinates"][-1]
        with open(filname, "w") as f:
            f.write(f"{len(xyz)}\n")
            f.write("\n")
            for ll in xyz:
                f.write(f"{ll}\n")
        return
