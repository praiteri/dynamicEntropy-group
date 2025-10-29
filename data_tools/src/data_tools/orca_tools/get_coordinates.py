import numpy as np

from ..orca_tools.extract_from_orca_output import extract_from_orca_output
from ..orca_tools.hierarchical_clustering import hierarchical_clustering


class get_coordinates(extract_from_orca_output):
    """
    Extract the coordinates
    """

    def read(self):
        string = "CARTESIAN COORDINATES (ANGSTROEM)"
        escape_string = [""]
        self.read_blocks(string, escape_string, offset=2)

    def print_formula(self):
        """
        Print the chemical formula
        """
        last_frame = self.result[-1]
        labels = [x.split()[0] for x in last_frame]
        from collections import Counter

        counter = Counter(labels)
        formula = ""
        for item, count in counter.items():
            if count > 1:
                formula += f"{item}{count} "
            else:
                formula += f"{item} "
        self.logger.info(self.fmt.format("Chemical formula", formula.strip()))

    def check_cluster_integrity(self):
        """
        Check the cluster integrity
        """
        last_frame = self.result[-1]
        cutoff_distance = (
            4  # maximum distance between two points to be considered in the same cluster
        )
        method = "single"  # linkage method for hierarchical clustering
        coordinates = np.ndarray(shape=(len(last_frame), 3))

        for i in range(len(last_frame)):
            coordinates[i, :] = np.array(last_frame[i].split()[1:4])

        cluster_labels = hierarchical_clustering(coordinates, cutoff_distance, method)
        if 2 in (cluster_labels):
            self.logger.warning(self.fmt.format("Cluster integrity check", "FAILED"))
            self.logger.verbose("  Atom clustering")
            self.logger.verbose(cluster_labels)
        else:
            self.logger.info(self.fmt.format("Cluster integrity check", "PASSED"))
            (self.fmt.format("Cluster integrity check", "PASSED"))

    def write(self, **kwargs):
        if self.result is None:
            return

        if "frame" in kwargs:
            nframe = kwargs["frame"]
        else:
            nframe = "last"

        if nframe is None:
            raise Exception("No frame specified")
        elif nframe == "last" or nframe == -1:
            nframe = -1
            self.logger.info("Writing final coordinates (XYZ)")
        elif nframe == "first":
            nframe = 0
            self.logger.info("Writing initial coordinates")
        else:
            self.logger.info(f"Coordinates for frame {nframe}")

        xyz = self.result[nframe]
        self.logger.info(f"Number of atoms = {len(xyz)}")
        self.logger.info("")
        for ll in xyz:
            self.logger.info(f"{ll}")
