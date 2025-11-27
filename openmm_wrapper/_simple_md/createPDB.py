#!/usr/bin/env python3

import random

def generate_pdb_file(num_atoms, cell, output_filename='argon_atoms.pdb'):
    # PDB header
    pdb_lines = ["HEADER    Argon Atoms\n"]
    pdb_lines = [f"CRYST1   {cell:6.3f}   {cell:6.3f}   {cell:6.3f}  90.00  90.00  90.00\n"]

    # Generate ATOM records for Argon atoms
    for i in range(1, num_atoms + 1):
        x = random.uniform(0, cell)
        y = random.uniform(0, cell)
        z = random.uniform(0, cell)
        atom_line = f"ATOM  {i:4}  Ar   ARG  {i:4}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          Ar\n"
        pdb_lines.append(atom_line)

    # End of file record
    pdb_lines.append("END\n")

    # Write to the output file
    with open(output_filename, 'w') as f:
        f.writelines(pdb_lines)

if __name__ == "__main__":
    num_atoms = 200  # Change this to the desired number of Argon atoms
    cell = 50 # in Angstrom
    output_filename='coord.pdb'
    generate_pdb_file(num_atoms,cell,output_filename=output_filename)
    print("PDB file {} created with {} Argon atoms.".format(output_filename,num_atoms))

