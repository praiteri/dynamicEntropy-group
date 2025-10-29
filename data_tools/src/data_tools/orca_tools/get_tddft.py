import subprocess as sp
import os
import sys
import shutil

from ..orca_tools.extract_from_orca_output import extract_from_orca_output
from ..logger import formatting
from ..orca_tools.utils import centered_printout


class get_tddft(extract_from_orca_output):
    """
    Extract the TD-DFT results
    """

    def read(self):
        idx = self.search_string("ORCA TD-DFT/TDA CALCULATION")
        # Not a geometry optimisation run
        if idx is None:
            self.result = None
            return

        self.result = {}
        self.result_strings = {
            "abs": "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS",
            "cd": "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS",
            # "abs_vel" : "CD SPECTRUM",
            # "cd_vel" : "CD SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS",
        }

        for k, v in self.result_strings.items():
            self.result[k] = self.parse_blocks(v, [""], 5)[0]

        chunk = self.parse_blocks(
            "TD-DFT/TDA EXCITED STATES (SINGLETS)", ["TD-DFT/TDA-EXCITATION SPECTRA"], 5
        )
        self.result["characters"] = self.parse_blocks("STATE", [""], lines=chunk[0])

        # Natual Transition orbitals
        idx = self.search_string("Natural Transition Orbitals were saved in")
        if idx is None:
            return

        self.result["nto"] = self.parse_blocks(
            "NATURAL TRANSITION ORBITALS FOR STATE", [""], 9
        )
        self.result["nto_files"] = self.parse_single_line(idx)
        return

    def write(self, **kwargs):
        if self.result is None:
            return

        na = int(round(float(kwargs["nalpha"].split()[2]), 0))
        nb = int(round(float(kwargs["nbeta"].split()[2]), 0))

        if na != nb:
            self.logger.error("Alpha and Beta electrons are different!")
            self.logger.error("Parser not yet implemented for this case!")
            sys.exit(1)

        self.logger.info(formatting().dashes)
        self.logger.info("TD-DFT/TDA results")
        self.logger.info(
            "  # | Wavelength | Wavenumber | energy (eV)| Int (a.u) | Main character"
        )

        all_fosc = [float(x.split()[6]) for x in self.result["abs"]]
        fscale = max(all_fosc) / 100.0
        for idx in range(len(self.result["abs"])):
            x = self.result["abs"][idx].split()
            erg_ev = x[3]
            erg_cm = x[4]
            erg_nm = x[5]

            # Normalisation
            fosc = float(x[6]) / fscale

            homo = na - 1
            lumo = na

            char = self.result["characters"][idx]
            trans = {}
            for i in range(1, len(char)):
                x = char[i].split(":")

                t = x[0].replace("a", "").split("->")
                n1 = int(t[0]) - homo
                n2 = int(t[1]) - lumo
                if n1 == 0:
                    s1 = "HOMO"
                else:
                    s1 = f"HOMO{n1}"

                if n2 == 0:
                    s2 = "LUMO"
                else:
                    s2 = f"LUMO+{n2}"

                s = centered_printout(f"{s1} -> {s2}")
                trans[s] = x[1].split()[0]

            # Order the transitions based on the oscillator strength
            trans = sorted(trans.items(), key=lambda x: x[1], reverse=True)

            self.logger.info(
                "{:>3d} | {:>10s} | {:>10s} | {:>10s} | {:>9.1f} | {:s}  ({:5.1f}% )".format(
                    idx + 1,
                    erg_nm,
                    erg_cm,
                    erg_ev,
                    fosc,
                    trans[0][0],
                    float(trans[0][1]) * 100,
                )
            )
            for j in range(1, len(trans)):
                self.logger.verbose(
                    "{:>56s} {:s}  ({:5.1f}% )".format(
                        "|", trans[j][0], float(trans[j][1]) * 100
                    )
                )
        return

    def compute_nto(self, **kwargs):
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

        fosc_threshold = 0
        if "fosc" in kwargs:
            fosc_threshold = float(kwargs["fosc"])
        if fosc_threshold < 1:
            fosc_threshold *= 100
        self.logger.info(
            self.fmt.format("Minimum Oscillator Strength (% of max)", fosc_threshold)
        )

        gridvalue = 100
        if "grid" in kwargs:
            gridvalue = kwargs["grid"]
        self.logger.info(self.fmt.format("Grid dimension", gridvalue))

        mo_operator = 0
        if "spin" in kwargs:
            mo_operator = kwargs["spin"]
        self.logger.info(self.fmt.format("Orbitals for spin", mo_operator))

        try:
            orcadir = os.path.dirname(shutil.which("orca_plot"))
            self.logger.info(self.fmt.format("ORCA directory", orcadir))
        except:
            self.logger.error("Cannot find ORCA executable")
            sys.exit(1)

        all_fosc = [float(x.split()[6]) for x in self.result["abs"]]
        fscale = max(all_fosc) / 100

        for idx in range(len(self.result["nto_files"])):
            filename = self.result["nto_files"][idx].split()[-1]
            fosc = float(self.result["abs"][idx].split()[6]) / fscale
            if fosc < fosc_threshold:
                continue
            self.logger.info(formatting().dashes_short)
            self.logger.info(self.fmt.format("Processing NTO file", filename))

            for jdx in self.result["nto"][idx][1:]:
                x = jdx.split()
                if float(x[-1]) < 0.1:
                    continue
                tt = jdx.split(":")
                self.logger.info(
                    self.fmt.format(
                        "Generating CUBE files for transition", f"{tt[0]} ({tt[1]})"
                    )
                )

                for i in [0, 2]:
                    state_file = filename.replace(".nto", f".mo{jdx.split()[i]}.vmd")
                    cube_file = filename.replace(".nto", f".mo{jdx.split()[i]}.cube")
                    self.write_vmd_state_file(statefile=state_file, cubefile=cube_file)

                for mo_number in [x[0].replace("a", ""), x[2].replace("a", "")]:
                    cmd = f"5\n7\n4\n{gridvalue}\n3\n{mo_operator}\n2\n{mo_number}\n11\n12\n\n"
                    output = sp.run(
                        [orcadir + "/orca_plot", filename, "-i"],
                        stdout=sp.PIPE,
                        input=cmd,
                        encoding="ascii",
                    )

                if output.returncode != 0:
                    self.logger.error("Failed to create CUBE files for the selected MO")
                    sys.exit(1)

        return

    def write_vmd_state_file(self, **kwargs):
        """
        Write VMD state file with molecular orbital visualization settings.

        Args:
            **kwargs: Keyword arguments including:
                - filename (str): Output filename for VMD state file (default: "state.vmd")
                - homo (str): HOMO cube file path (default: "homo.cube")
                - lumo (str): LUMO cube file path (default: "lumo.cube")
        """
        script = """
color change rgb  0 0.10 0.20 0.70 ;# blue
color change rgb  1 0.70 0.20 0.10 ;# red
color change rgb  2 0.49 0.49 0.49 ;# gray
color change rgb  3 0.93 0.50 0.15 ;# orange
color change rgb  4 1.00 1.00 0.00 ;# yellow
color change rgb  5 0.50 0.50 0.20 ;# tan
color change rgb  6 0.60 0.60 0.60 ;# silver
color change rgb  7 0.10 0.70 0.20 ;# green
color change rgb  8 1.00 1.00 1.00 ;# white
color change rgb  9 1.00 0.60 0.60 ;# pink
color change rgb 10 0.10 0.70 0.80 ;# cyan
color change rgb 11 0.50 0.16 0.50 ;# purple
color change rgb 12 0.39 0.86 0.12 ;# lime
color change rgb 13 0.90 0.40 0.70 ;# mauve
color change rgb 14 0.50 0.30 0.00 ;# ochre
color change rgb 15 0.50 0.50 0.75 ;# iceblue
color change rgb 16 0.60 0.10 0.60 ;# black
color change rgb 17 0.82 0.88 0.00 ;# yellow2
color change rgb 18 0.55 0.90 0.02 ;# yellow3
color change rgb 19 0.10 0.47 0.22 ;# green2
color change rgb 20 0.00 0.90 0.50 ;# green3
color change rgb 21 0.00 0.88 1.00 ;# cyan2
color change rgb 22 0.00 0.76 1.00 ;# cyan3
color change rgb 23 0.02 0.38 0.67 ;# blue2
color change rgb 24 0.01 0.04 0.93 ;# blue3
color change rgb 25 0.27 0.00 0.98 ;# violet
color change rgb 26 0.45 0.00 0.90 ;# violet2
color change rgb 27 0.90 0.00 0.90 ;# magenta
color change rgb 28 1.00 0.00 0.66 ;# magenta2
color change rgb 29 0.98 0.00 0.23 ;# red2
color change rgb 30 0.81 0.00 0.00 ;# red3
color change rgb 31 0.87 0.67 0.20 ;# orange2
color change rgb 32 0.96 0.72 0.00 ;# orange3

color Display Background white
color Labels Atoms blue
color Labels Bonds blue
color Labels Angles blue
color Labels Dihedrals blue

after idle {{
    color Element C  gray
    color Element N  blue
    color Element O  red
    color Element F  yellow2
    color Element P  orange
    color Element S  orange2
    color Element Si yellow
    color Element Cl lime
    color Element Co blue
    color Element W  pink
    color Element Ir purple
    color Element Lu green2
}}

label textsize 1.0

light 0 on
light 1 on
light 2 off
light 3 off

# Display settings
display eyesep       0.065000
display focallength  2.000000
display height       6.000000
display distance     -2.000000
display projection   Orthographic
display nearclip set 0.001000
display farclip  set 10.000000
display depthcue   off
display cuestart   0.500000
display cueend     10.000000
display cuestart   0.500000
display cueend     10.000000
display cuedensity 0.320000
display cuemode    Exp2
display shadows off
display ambientocclusion off
display aoambient 0.800000
display aodirect 0.300000
display dof off
display dof_fnumber 64.000000
display dof_focaldist 0.700000
display resize 1200 1200

mol new {cube_file} type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation CPK 1.000000 0.300000 32.000000 32.000000
mol color Element
mol addrep top
mol representation Isosurface 0.030000 0 0 0 1 1
mol color ColorID 15
mol addrep top
mol representation Isosurface -0.030000 0 0 0 1 1
mol color ColorID 21
mol addrep top
mol rename top {cube_file}

axes location off

# Simple list approach for your specific atoms
set metal_atoms {{Ir Re Lu Ln}}
foreach atom $metal_atoms {{
    set sel [atomselect top "name $atom"]
    $sel set radius 2.5
    $sel delete
}}
mol bondsrecalc [molinfo top]

mol bondsrecalc top
topo retypebonds

menu main move 1500 196
menu graphics move 1500 455
menu display move 1000 100

lappend viewplist [molinfo top]
set topmol [molinfo top]
# done with molecule 0
foreach v $viewplist {{
    molinfo $v set {{center_matrix rotate_matrix scale_matrix global_matrix}} $viewpoints($v)
}}

#MAKE_IMAGE display samples 4
#MAKE_IMAGE render TachyonInternal high_quality.tga
#MAKE_IMAGE exec sips -s format tiff high_quality.tga --out high_quality.tiff
#MAKE_IMAGE file delete high_quality.tga
#MAKE_IMAGE quit
        """

        filename = kwargs.get("statefile", "state.vmd")
        self.logger.info(self.fmt.format("Writing VMD state file", filename))

        cube_file = kwargs.get("cubefile", None)

        with open(filename, "w") as file:
            file.write(script.format(cube_file=cube_file))
