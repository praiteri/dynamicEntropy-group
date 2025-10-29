import copy
import numpy as np

from ..orca_tools.extract_from_orca_output import extract_from_orca_output
from ..logger import formatting


class get_vibrational_modes(extract_from_orca_output):
    """
    Extract the input commands
    """

    def read(self):
        string = "VIBRATIONAL FREQUENCIES"
        escape_string = ""
        self.read_blocks(string, escape_string, offset=5)

    def write(self, **kwargs):
        if self.result is None:
            return

        self.logger.info(formatting().dashes)
        modes = self.result[0]
        nmodes = len(modes)
        self.logger.info(self.fmt.format("Number of vibrational modes", nmodes))
        val = np.array([float(modes[i].split(":")[-1].split()[0]) for i in range(0, 6)])
        if not (val == 0).all():
            self.logger.error("The first 6 vibrational modes are not zero!")

        val = np.array([float(modes[i].split(":")[-1].split()[0]) for i in range(6, nmodes)])
        negative_freqs = np.count_nonzero(val < 0)

        if negative_freqs > 0:
            self.logger.warning(f"There are {negative_freqs} negative frequencies!")
            self.logger.warning("The imaginary frequencies are:")
            for v in modes[6 : 6 + negative_freqs]:
                self.logger.warning(v.split(":")[-1])
        else:
            self.logger.info("There are no negative frequencies!")
            self.logger.verbose("The first 5 non-zero frequencies are:")
            for v in modes[6:11]:
                self.logger.verbose(v.split(":")[-1])
            self.logger.debug("The remaining frequencies are:")
            for v in modes[11:]:
                self.logger.debug(v.split(":")[-1])

    def get_normal_modes(self):
        text_block = self.parse_blocks("NORMAL MODES", [""], 7)[0]

        number_of_modes = len(self.result[0])
        normal_modes = np.zeros(shape=(number_of_modes, number_of_modes))
        number_of_blocks = int(number_of_modes / 6) + min(1, number_of_modes % 6)

        iLine = 0
        for _ in range(number_of_blocks):
            mode_indices = [int(x) for x in text_block[iLine].split()]
            for kk in range(number_of_modes):
                iLine += 1
                ll = text_block[iLine].split()
                jdx = int(ll[0])

                for k, idx in enumerate(mode_indices):
                    normal_modes[idx, jdx] = float(ll[k + 1])
            iLine += 1

        return normal_modes

    def animate(self, **kwargs):
        if self.result is None:
            return
        self.logger.info(formatting().dashes)

        if "mode" in kwargs:
            if not isinstance(kwargs["mode"], list):
                kwargs["mode"] = [kwargs["mode"]]
            list_of_modes = kwargs["mode"]
        else:
            raise Exception("No vibrational mode specified; use -a mode=12")
        strings = [f"{x}" for x in list_of_modes]
        self.logger.info(self.fmt.format("Vibrational modes animated", f"{','.join(strings)}"))

        string = "CARTESIAN COORDINATES (ANGSTROEM)"
        escape_string = [""]
        last_frame = self.parse_blocks(string, escape_string, offset=2)[-1]
        number_of_atoms = len(last_frame)
        labels = [x.split()[0] for x in last_frame]
        coordinates = np.array(
            [last_frame[i].split()[1:4] for i in range(number_of_atoms)], dtype=float
        )

        normal_modes = self.get_normal_modes()

        nframes = 41
        if "nframes" in kwargs:
            nframes = int(kwargs["nframes"])
        if nframes % 2 == 0:
            nframes += 1
        self.logger.info(self.fmt.format("Number of frames (always an odd number)", nframes))

        max_disp = 0.1
        if "max_disp" in kwargs:
            max_disp = float(kwargs["max_disp"])
        self.logger.info(self.fmt.format("Maximum displacement (a.u.)", max_disp))

        for mode in list_of_modes:
            output_file = f"freq_{mode}.xyz"
            if "output" in kwargs:
                x = kwargs["output"].split(".")
                output_file = f"{x[0]}_{mode}.{x[1]}"

            displacement = normal_modes[int(mode), :].reshape(number_of_atoms, 3)

            output_xyz = open(output_file, "w")
            self.logger.info(self.fmt.format("Coordinates file", output_file))

            c = copy.deepcopy(coordinates)
            for mm in np.linspace(-max_disp, max_disp, nframes):
                output_xyz.write(f"{number_of_atoms}\n")
                output_xyz.write("\n")
                for i in range(number_of_atoms):
                    pos = c[i]
                    delta = displacement[i] * mm
                    string = "{:6s} {:10.7f} {:10.7f} {:10.7f}\n".format(labels[i], *(pos + delta))
                    output_xyz.write(string)

    def follow_negative_modes(self, **kwargs):
        if self.result is None:
            return

        frequencies = self.result[0]
        number_of_frequencies = len(frequencies)

        normal_modes = self.get_normal_modes()

        self.logger.info(formatting().dashes)
        if "mode" in kwargs:
            if not isinstance(kwargs["mode"], list):
                kwargs["mode"] = [kwargs["mode"]]
            list_of_modes = kwargs["mode"]
        else:
            self.logger.info("Following all negative vibrational modes")
            list_of_modes = [
                x for x in range(6, number_of_frequencies) if float(frequencies[x].split()[1]) < 0
            ]

        string = ",".join([str(x) for x in list_of_modes])
        self.logger.info(self.fmt.format("Vibrational modes followed", string))

        string = "CARTESIAN COORDINATES (ANGSTROEM)"
        escape_string = [""]
        last_frame = self.parse_blocks(string, escape_string, offset=2)[-1]
        number_of_atoms = len(last_frame)
        labels = [x.split()[0] for x in last_frame]
        coordinates = np.array(
            [last_frame[i].split()[1:4] for i in range(number_of_atoms)], dtype=float
        )

        factor = 1.0
        if "scale" in kwargs:
            factor = float(kwargs["scale"])
        self.logger.info(self.fmt.format("Scaling factor", factor))

        for mode in list_of_modes:
            displacement = normal_modes[int(mode), :].reshape(number_of_atoms, 3) * factor
            for i in range(number_of_atoms):
                coordinates[i, :] += displacement[i]

        output_file = "new_coordinates.xyz"
        if "output" in kwargs:
            output_file = kwargs["output"]
        self.logger.info(self.fmt.format("New coordinates file", output_file))

        output_xyz = open(output_file, "w")
        output_xyz.write(f"{number_of_atoms}\n")
        output_xyz.write("\n")
        for i in range(number_of_atoms):
            pos = coordinates[i]
            output_xyz.write("{:6s} {:10.7f} {:10.7f} {:10.7f}\n".format(labels[i], *(pos)))
