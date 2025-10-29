import subprocess as sp
import os
import sys
import shutil

from ..orca_tools.extract_from_orca_output import extract_from_orca_output
from ..logger import formatting
from ..orca_tools.get_values import get_values


class get_orbital_energies(extract_from_orca_output):
    """
    Extract the input commands
    """

    def read(self):
        self.result = {}
        self.result["gbw_file"] = self.parse_single_line(self.search_string("GBWName"))
        self.result["gbw_file"] = self.result["gbw_file"][-1]

        self.result["nel"] = round(
            int(get_values(self.lines, "Number of Electrons    NEL").result[-1].split()[-1])
        )

        hf_type = get_values(self.lines, "Hartree-Fock type      HFTyp").result[-1].split()[-1]
        if hf_type == "RHF":
            self.unrestricted = False
        else:
            self.unrestricted = True

        try:
            self.result["nalpha"] = round(
                float(get_values(self.lines, "N(Alpha)").result[-1].split()[-2])
            )
            self.result["nbeta"] = round(
                float(get_values(self.lines, "N(Beta)").result[-1].split()[-2])
            )
        except:
            if not self.unrestricted:
                self.result["nalpha"] = int(self.result["nel"] / 2)
                self.result["nbeta"] = self.result["nalpha"]
            else:
                # mult = 2S+1 or N+1
                res = get_values(self.lines, "Multiplicity           Mult", text="Multiplicity")
                multiplicity = int(res.result[-1].split()[-1])
                n_unpaired = int(multiplicity) - 1
                ntmp = self.result["nel"] - n_unpaired
                self.result["nalpha"] = int(ntmp / 2 + n_unpaired)
                self.result["nbeta"] = int(ntmp / 2)

        # check if it is a restricted or unrestricted calculation
        if self.unrestricted:
            alpha = self.search_string("SPIN UP ORBITALS")
            beta = self.search_string("SPIN DOWN ORBITALS")

            string = "SPIN UP ORBITALS"
            escape_string = ["*Only", ""]
            alpha = self.parse_blocks(string, escape_string, offset=4)

            string = "SPIN DOWN ORBITALS"
            escape_string = ["*Only", ""]
            beta = self.parse_blocks(string, escape_string, offset=4)
            self.result["orbitals"] = [alpha, beta]

        else:
            string = "ORBITAL ENERGIES"
            escape_string = ["*Only", ""]
            alpha = self.parse_blocks(string, escape_string, offset=4)
            self.result["orbitals"] = [alpha]

    def write(self, **kwargs):
        self.logger.info(formatting().dashes)
        nalpha = self.result["nalpha"]
        # nbeta = self.result["nbeta"]
        if len(self.result) == 0:
            self.logger.error("No orbital energies found!")
            sys.exit(1)

        if not self.unrestricted:
            orbitals = self.result["orbitals"][0][0]
            self.logger.info("Restricted calculation")
            self.logger.info("Only Alpha orbitals are available")
            homo = nalpha - 1
            lumo = homo + 1
            homo_ev = orbitals[homo].split()[-1]
            lumo_ev = orbitals[lumo].split()[-1]
            for i in range(-3, 0):
                erg = orbitals[homo + i].split()[-1]
                self.logger.verbose(self.fmt.format(f"HOMO{i} energy (eV)", erg))
            self.logger.info(self.fmt.format("HOMO energy (eV)", homo_ev))
            self.logger.info(self.fmt.format("LUMO energy (eV)", lumo_ev))
            for i in range(1, 4):
                erg = orbitals[lumo + i].split()[-1]
                self.logger.verbose(self.fmt.format(f"LUMO+{i} energy (eV)", erg))

        elif self.unrestricted:
            self.logger.info("Unrestricted calculation")
            self.logger.info("Alpha and Beta orbitals are available")
            ele = ["alpha", "beta"]

            for i in range(2):
                orbitals = self.result["orbitals"][i][0]
                homo = nalpha - 1
                lumo = homo + 1
                homo_ev = orbitals[homo].split()[-1]
                lumo_ev = orbitals[lumo].split()[-1]
                self.logger.debug(f"{ele[i]} electrons - HOMO={homo} | LUMO={lumo}")
                self.logger.info(self.fmt.format(f"{ele[i]}-HOMO energy (eV)", homo_ev))
                self.logger.info(self.fmt.format(f"{ele[i]}-LUMO energy (eV)", lumo_ev))

        else:
            self.logger.error("MOLECULAR ORBOTALS, SOMETHING WENT WRONG")
            sys.exit(1)

    def compute_mo(self, **kwargs):
        """
        1 - Enter type of plot
        2 - Enter no of orbital to plot
        3 - Enter operator of orbital (0=alpha,1=beta)
        4 - Enter number of grid intervals
        5 - Select output file format
        6 - Plot CIS/TD-DFT difference densities
        7 - Plot CIS/TD-DFT transition densities
        8 - Set AO(=1) vs MO(=0) to plot
        9 - List all available densities
        10 - Perform Density Algebraic Operations

        11 - Generate the plot
        12 - exit this program
        """
        nalpha = self.result["nalpha"]
        # nbeta = self.result["nbeta"]

        try:
            orcadir = os.path.dirname(shutil.which("orca_plot"))
        except:
            self.logger.error("Cannot find orca_plot executable")
            return

        gridvalue = 50
        mo_operator = 0

        filename = self.result["gbw_file"].split()[-1]

        homo = nalpha - 1
        lumo = nalpha

        mo_list = []
        if "range" in kwargs:
            norb = int(kwargs["range"])
            mo_list.extend([x for x in range(nalpha - norb, nalpha + norb)])

        if "orbitals" in kwargs:
            mo_list.extend(kwargs["orbitals"].split(","))

        if len(mo_list) == 0:
            mo_list = ["homo-1", "homo", "lumo", "lumo+1"]

        if not self.unrestricted:
            homo = nalpha - 1
            lumo = nalpha

            for mo in mo_list:
                if type(mo) is int:
                    mo_number = mo
                else:
                    if "homo" in mo:
                        mo_number = homo
                        try:
                            mo_number = homo + int(mo.replace("homo", ""))
                        except:
                            pass

                    elif "lumo" in mo:
                        mo_number = lumo
                        try:
                            mo_number = lumo + int(mo.replace("lumo", ""))
                        except:
                            pass

                    else:
                        self.logger.error(f"Cannot understand orbital {mo}")
                        sys.error(1)

                if mo_number < homo:
                    cube_file = f"HOMO{mo_number - homo}"
                elif mo_number == homo:
                    cube_file = "HOMO"
                elif mo_number == lumo:
                    cube_file = "LUMO"
                elif mo_number > lumo:
                    cube_file = f"LUMO+{mo_number - lumo}"

                self.logger.info(f"Generating cube file for obital number {mo} ({cube_file})")

                orca_input = f"5\n7\n4\n{gridvalue}\n3\n{mo_operator}\n2\n{mo_number}\n11\n12\n\n"
                cmd = [orcadir + "/orca_plot", filename, "-i"]
                output = sp.run(cmd, stdout=sp.PIPE, input=orca_input, encoding="ascii")

                if output.returncode != 0:
                    self.logger.error("Failed to create CUBE files for the selected MO")
                    sys.exit(1)
                else:
                    sp.run(f"mv orca.mo{mo}a.cube {cube_file}.cube", shell=True)

        elif self.unrestricted:
            return

        return
