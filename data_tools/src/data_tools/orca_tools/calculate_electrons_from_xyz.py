def calculate_total_electrons(xyz_file_path):
    """
    Read an XYZ file and calculate the total number of electrons in the system.

    Parameters:
    xyz_file_path (str): Path to the XYZ file

    Returns:
    int: Total number of electrons in the system
    """
    # Dictionary containing elements and their electron counts

    periodic_table = {
        "H": {"name": "Hydrogen", "atomic_number": 1, "atomic_mass": 1.008},
        "He": {"name": "Helium", "atomic_number": 2, "atomic_mass": 4.0026},
        "Li": {"name": "Lithium", "atomic_number": 3, "atomic_mass": 6.94},
        "Be": {"name": "Beryllium", "atomic_number": 4, "atomic_mass": 9.0122},
        "B": {"name": "Boron", "atomic_number": 5, "atomic_mass": 10.81},
        "C": {"name": "Carbon", "atomic_number": 6, "atomic_mass": 12.011},
        "N": {"name": "Nitrogen", "atomic_number": 7, "atomic_mass": 14.007},
        "O": {"name": "Oxygen", "atomic_number": 8, "atomic_mass": 15.999},
        "F": {"name": "Fluorine", "atomic_number": 9, "atomic_mass": 18.998},
        "Ne": {"name": "Neon", "atomic_number": 10, "atomic_mass": 20.180},
        "Na": {"name": "Sodium", "atomic_number": 11, "atomic_mass": 22.990},
        "Mg": {"name": "Magnesium", "atomic_number": 12, "atomic_mass": 24.305},
        "Al": {"name": "Aluminum", "atomic_number": 13, "atomic_mass": 26.982},
        "Si": {"name": "Silicon", "atomic_number": 14, "atomic_mass": 28.085},
        "P": {"name": "Phosphorus", "atomic_number": 15, "atomic_mass": 30.974},
        "S": {"name": "Sulfur", "atomic_number": 16, "atomic_mass": 32.06},
        "Cl": {"name": "Chlorine", "atomic_number": 17, "atomic_mass": 35.45},
        "Ar": {"name": "Argon", "atomic_number": 18, "atomic_mass": 39.948},
        "K": {"name": "Potassium", "atomic_number": 19, "atomic_mass": 39.098},
        "Ca": {"name": "Calcium", "atomic_number": 20, "atomic_mass": 40.078},
        "Sc": {"name": "Scandium", "atomic_number": 21, "atomic_mass": 44.956},
        "Ti": {"name": "Titanium", "atomic_number": 22, "atomic_mass": 47.867},
        "V": {"name": "Vanadium", "atomic_number": 23, "atomic_mass": 50.942},
        "Cr": {"name": "Chromium", "atomic_number": 24, "atomic_mass": 51.996},
        "Mn": {"name": "Manganese", "atomic_number": 25, "atomic_mass": 54.938},
        "Fe": {"name": "Iron", "atomic_number": 26, "atomic_mass": 55.845},
        "Co": {"name": "Cobalt", "atomic_number": 27, "atomic_mass": 58.933},
        "Ni": {"name": "Nickel", "atomic_number": 28, "atomic_mass": 58.693},
        "Cu": {"name": "Copper", "atomic_number": 29, "atomic_mass": 63.546},
        "Zn": {"name": "Zinc", "atomic_number": 30, "atomic_mass": 65.38},
        "Ga": {"name": "Gallium", "atomic_number": 31, "atomic_mass": 69.723},
        "Ge": {"name": "Germanium", "atomic_number": 32, "atomic_mass": 72.630},
        "As": {"name": "Arsenic", "atomic_number": 33, "atomic_mass": 74.922},
        "Se": {"name": "Selenium", "atomic_number": 34, "atomic_mass": 78.971},
        "Br": {"name": "Bromine", "atomic_number": 35, "atomic_mass": 79.904},
        "Kr": {"name": "Krypton", "atomic_number": 36, "atomic_mass": 83.798},
        "Rb": {"name": "Rubidium", "atomic_number": 37, "atomic_mass": 85.468},
        "Sr": {"name": "Strontium", "atomic_number": 38, "atomic_mass": 87.62},
        "Y": {"name": "Yttrium", "atomic_number": 39, "atomic_mass": 88.906},
        "Zr": {"name": "Zirconium", "atomic_number": 40, "atomic_mass": 91.224},
        "Nb": {"name": "Niobium", "atomic_number": 41, "atomic_mass": 92.906},
        "Mo": {"name": "Molybdenum", "atomic_number": 42, "atomic_mass": 95.95},
        "Tc": {"name": "Technetium", "atomic_number": 43, "atomic_mass": 98.0},
        "Ru": {"name": "Ruthenium", "atomic_number": 44, "atomic_mass": 101.07},
        "Rh": {"name": "Rhodium", "atomic_number": 45, "atomic_mass": 102.91},
        "Pd": {"name": "Palladium", "atomic_number": 46, "atomic_mass": 106.42},
        "Ag": {"name": "Silver", "atomic_number": 47, "atomic_mass": 107.87},
        "Cd": {"name": "Cadmium", "atomic_number": 48, "atomic_mass": 112.41},
        "In": {"name": "Indium", "atomic_number": 49, "atomic_mass": 114.82},
        "Sn": {"name": "Tin", "atomic_number": 50, "atomic_mass": 118.71},
        "Sb": {"name": "Antimony", "atomic_number": 51, "atomic_mass": 121.76},
        "Te": {"name": "Tellurium", "atomic_number": 52, "atomic_mass": 127.60},
        "I": {"name": "Iodine", "atomic_number": 53, "atomic_mass": 126.90},
        "Xe": {"name": "Xenon", "atomic_number": 54, "atomic_mass": 131.29},
        "Cs": {"name": "Cesium", "atomic_number": 55, "atomic_mass": 132.91},
        "Ba": {"name": "Barium", "atomic_number": 56, "atomic_mass": 137.33},
        "La": {"name": "Lanthanum", "atomic_number": 57, "atomic_mass": 138.91},
        "Ce": {"name": "Cerium", "atomic_number": 58, "atomic_mass": 140.12},
        "U": {"name": "Uranium", "atomic_number": 92, "atomic_mass": 238.03},
    }

    try:
        with open(xyz_file_path, "r") as file:
            lines = file.readlines()

        # First line is the number of atoms
        num_atoms = int(lines[0].strip())

        elements = [lines[i].strip().split()[0] for i in range(2, 2 + num_atoms)]
        total_electrons = sum([periodic_table[element]["atomic_number"] for element in elements])

        return total_electrons

    except FileNotFoundError:
        print(f"Error: File {xyz_file_path} not found")
        return None
    except Exception as e:
        print(f"Error reading XYZ file: {str(e)}")
        return None
