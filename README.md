# Density-of-state
This script help you to obtain the density of state (DOS) of your structure from VASP results. Currently it can deal with 5 types of atom at most. If you have more than 5 types of atoms in your poscar. we can slightly change the code to obtain the DOS.
1. Putting the DOSCAR and POSCAR into this folder.
2. Running the script and it will generate "dos_total.txt" which contain 2 columns (1st column is energy, 2nd column is total dos of the compound); It will also generate a serise of "AtomName_dos.txt" which contain 5 columns (1st column is energy, 2nd column is total atomic dos, 3rd column is total s dos, 4th column is total p dos, 5th column is total d dos).
