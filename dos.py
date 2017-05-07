##Edited by Jin Zhang 07-13-2016
##This script can obtain DOS from VASP calcualtion. This code can deal with 5 types of atom at most.
##If you have more than 5 types of atoms in your poscar. we can slightly change the code to obtain the DOS.
## Spin=1 in INCAR. From line which including 10 columns in DOSCAR. These 10 columns represent energy, s, py, pz, px, dxy, dyz, dz2, dxz, dx2.
## The export files including the "dos_total.txt" which contain 2 columns (1st column is energy, 2nd column is total dos of the compound);
## "AtomName_dos.txt" which contain 5 columns (1st column is energy, 2nd column is total atomic dos, 3rd column is total s dos, 4th column is total p dos, 5th column is total d dos).
import numpy as np
import matplotlib.pyplot as plt


data = []
n=[]
with open ("DOSCAR") as doscar:
    data_0 = doscar.readlines()
    data_1 = data_0 [5]
    data_2 = data_0 [0]
    k_num = int (data_1.split()[2])
    e_fermi = float (data_1.split()[3])
    atom_num = int (data_2.split()[1])

with open ("DOSCAR") as doscar:
    content = [x.rstrip() for x in doscar]
    data_new = [x[4:] for x in content [6:]]
    for i in range(k_num):
        m = data_new[i]
        data.extend(m.split())
data = np.asarray(data, dtype = np.float64)
data = data.reshape(k_num, 3)
energy = data[:, 0]
dos_total = data[:, 1]
dos_total = np.asarray (dos_total, dtype = np.float64)
dos_total = dos_total.reshape(k_num, 1)
energy = [i-e_fermi for i in energy] # list
energy = np.asarray(energy, dtype = np.float64)
energy = energy.reshape(k_num, 1)
E_dos_total = np.hstack((energy, dos_total))
np.savetxt("dos_total.txt", E_dos_total)
print (k_num)

with open ("DOSCAR") as doscar:
    content = [x.rstrip() for x in doscar]
    data_new = [x[4:] for x in content [k_num+7:]]
    colume_num = len(data_new[0].split())
    for i in range(len(data_new)):
        if len(data_new[i].split()) > int(5):
            n.extend(data_new[i].split())
n = np.asarray(n, dtype = np.float64)
l = int(len(n)/colume_num)
num_dos = int(l/k_num)
n = n.reshape(l, colume_num)
dos_data = n [:, 1:]
u = dos_data.shape
dos_split = dos_data.reshape(num_dos, k_num, u[1])  # split dos and put them into same array

###########################
with open ("POSCAR") as poscar:
    poscar_0 = poscar.readlines()
    name_atom = poscar_0 [5]
    poscar_atom = poscar_0 [6]
    atom_type_num = len(poscar_atom.split())
    atom_type = poscar_atom.split()
    first_atom_num = int (poscar_atom.split()[0])
    second_atom_num = int (poscar_atom.split()[1])
##################################
dos_atom_1 = np.zeros((k_num, u[1]))
dos_atom_2 = np.zeros((k_num, u[1]))
dos_atom_3 = np.zeros((k_num, u[1]))
dos_atom_4 = np.zeros((k_num, u[1]))
dos_atom_5 = np.zeros((k_num, u[1]))

if atom_type_num==1:             # have 1 type of atom
    name_atom_1 = name_atom.split()[0]
    atom_type_1 = int(atom_type[0])
    for i in range(atom_type_1):
         dos_atom_11 = dos_split[i]
         dos_atom_1 = dos_atom_11 + dos_atom_1
    atom1_s =  (dos_atom_1[:, 0]).reshape(k_num, 1)
    atom1_p = (dos_atom_1[:, 1]+dos_atom_1[:, 2]+dos_atom_1[:, 3]).reshape(k_num, 1)
    atom1_d = (dos_atom_1[:, 4]+dos_atom_1[:, 5]+dos_atom_1[:, 6]+dos_atom_1[:, 7]+dos_atom_1[:, 8]).reshape(k_num, 1)
    atom1_total = atom1_s + atom1_p + atom1_d
    atom1_dos = np.hstack((energy, atom1_total, atom1_s, atom1_p, atom1_d)) # Energy, atom_total_dos, atom_s, atom_p, atom_d.
    np.savetxt( name_atom_1+"_dos.txt", atom1_dos)
    
################################    
if atom_type_num==2:           # have 2 types of atom
    name_atom_1 = name_atom.split()[0] 
    name_atom_2 = name_atom.split()[1]
    atom_type_1 = int(atom_type[0])
    atom_type_2 = int(atom_type[1])
    for i in range(atom_type_1):
        dos_atom_11 = dos_split[i]
        dos_atom_1 = dos_atom_11 + dos_atom_1
    atom1_s =  (dos_atom_1[:, 0]).reshape(k_num, 1)
    atom1_p = (dos_atom_1[:, 1]+dos_atom_1[:, 2]+dos_atom_1[:, 3]).reshape(k_num, 1)
    atom1_d = (dos_atom_1[:, 4]+dos_atom_1[:, 5]+dos_atom_1[:, 6]+dos_atom_1[:, 7]+dos_atom_1[:, 8]).reshape(k_num, 1)
    atom1_total = atom1_s + atom1_p + atom1_d
    atom1_dos = np.hstack((energy, atom1_total, atom1_s, atom1_p, atom1_d)) # Energy, atom_total_dos, atom_s, atom_p, atom_d.
    np.savetxt( name_atom_1+"_dos.txt", atom1_dos)
    for i in range(atom_type_1, atom_type_1+atom_type_2):
        dos_atom_22 = dos_split[i]
        dos_atom_2 = dos_atom_22 + dos_atom_2
    atom2_s = (dos_atom_2[:, 0]).reshape(k_num, 1)
    atom2_p = (dos_atom_2[:, 1]+dos_atom_2[:, 2]+dos_atom_2[:, 3]).reshape(k_num, 1)
    atom2_d = (dos_atom_2[:, 4]+dos_atom_2[:, 5]+dos_atom_2[:, 6]+dos_atom_2[:, 7]+dos_atom_2[:, 8]).reshape(k_num, 1)
    atom2_total = atom2_s + atom2_p + atom2_d
    atom2_dos = np.hstack((energy, atom2_total, atom2_s, atom2_p, atom2_d)) # Energy, atom_total_dos, atom_s, atom_p, atom_d.
    np.savetxt( name_atom_2+"_dos.txt", dos_atom_2)
    print (dos_atom_2.shape)
    
#################################
if atom_type_num==3:          # have 3 types of atom
    name_atom_1 = name_atom.split()[0]
    name_atom_2 = name_atom.split()[1]
    name_atom_3 = name_atom.split()[2]
    atom_type_1 = int(atom_type[0])
    atom_type_2 = int(atom_type[1])
    atom_type_3 = int(atom_type[2])
    for i in range(atom_type_1):
         dos_atom_11 = dos_split[i]
         dos_atom_1 = dos_atom_11 + dos_atom_1
    atom1_s =  (dos_atom_1[:, 0]).reshape(k_num, 1)
    atom1_p = (dos_atom_1[:, 1]+dos_atom_1[:, 2]+dos_atom_1[:, 3]).reshape(k_num, 1)
    atom1_d = (dos_atom_1[:, 4]+dos_atom_1[:, 5]+dos_atom_1[:, 6]+dos_atom_1[:, 7]+dos_atom_1[:, 8]).reshape(k_num, 1)
    atom1_total = atom1_s + atom1_p + atom1_d
    atom1_dos = np.hstack((energy, atom1_total, atom1_s, atom1_p, atom1_d)) # Energy, atom_total_dos, atom_s, atom_p, atom_d.
    np.savetxt( name_atom_1+"_dos.txt", atom1_dos)
    for i in range(atom_type_1, atom_type_1+atom_type_2):
         dos_atom_22 = dos_split[i]
         dos_atom_2 = dos_atom_22 + dos_atom_2
    atom2_s = (dos_atom_2[:, 0]).reshape(k_num, 1)
    atom2_p = (dos_atom_2[:, 1]+dos_atom_2[:, 2]+dos_atom_2[:, 3]).reshape(k_num, 1)
    atom2_d = (dos_atom_2[:, 4]+dos_atom_2[:, 5]+dos_atom_2[:, 6]+dos_atom_2[:, 7]+dos_atom_2[:, 8]).reshape(k_num, 1)
    atom2_total = atom2_s + atom2_p + atom2_d
    atom2_dos = np.hstack((energy, atom2_total, atom2_s, atom2_p, atom2_d)) # Energy, atom_total_dos, atom_s, atom_p, atom_d.
    np.savetxt( name_atom_2+"_dos.txt", dos_atom_2)
    for i in range(atom_type_1+atom_type_2, atom_type_1+atom_type_2+atom_type_3):
         dos_atom_33 = dos_split[i]
         dos_atom_3 = dos_atom_33 + dos_atom_3
    atom3_s = (dos_atom_3[:, 0]).reshape(k_num, 1)
    atom3_p = (dos_atom_3[:, 1]+dos_atom_3[:, 2]+dos_atom_3[:, 3]).reshape(k_num, 1)
    atom3_d = (dos_atom_3[:, 4]+dos_atom_3[:, 5]+dos_atom_3[:, 6]+dos_atom_3[:, 7]+dos_atom_3[:, 8]).reshape(k_num, 1)
    atom3_total = atom3_s + atom3_p + atom3_d
    atom3_dos = np.hstack((energy, atom3_total, atom3_s, atom3_p, atom3_d)) # Energy, atom_total_dos, atom_s, atom_p, atom_d.
    np.savetxt( name_atom_3+"_dos.txt", dos_atom_3)
    
################################
if atom_type_num==4:                # have 4 types of atom
    name_atom_1 = name_atom.split()[0]
    name_atom_2 = name_atom.split()[1]
    name_atom_3 = name_atom.split()[2]
    name_atom_4 = name_atom.split()[3]
    atom_type_1 = int(atom_type[0])
    atom_type_2 = int(atom_type[1])
    atom_type_3 = int(atom_type[2])
    atom_type_4 = int(atom_type[3])
    for i in range(atom_type_1):
        dos_atom_11 = dos_split[i]
        dos_atom_1 = dos_atom_11 + dos_atom_1
    atom1_s =  (dos_atom_1[:, 0]).reshape(k_num, 1)
    atom1_p = (dos_atom_1[:, 1]+dos_atom_1[:, 2]+dos_atom_1[:, 3]).reshape(k_num, 1)
    atom1_d = (dos_atom_1[:, 4]+dos_atom_1[:, 5]+dos_atom_1[:, 6]+dos_atom_1[:, 7]+dos_atom_1[:, 8]).reshape(k_num, 1)
    atom1_total = atom1_s + atom1_p + atom1_d
    atom1_dos = np.hstack((energy, atom1_total, atom1_s, atom1_p, atom1_d)) # Energy, atom_total_dos, atom_s, atom_p, atom_d.
    np.savetxt( name_atom_1+"_dos.txt", atom1_dos)
    for i in range(atom_type_1, atom_type_1+atom_type_2):
        dos_atom_22 = dos_split[i]
        dos_atom_2 = dos_atom_22 + dos_atom_2
    atom2_s = (dos_atom_2[:, 0]).reshape(k_num, 1)
    atom2_p = (dos_atom_2[:, 1]+dos_atom_2[:, 2]+dos_atom_2[:, 3]).reshape(k_num, 1)
    atom2_d = (dos_atom_2[:, 4]+dos_atom_2[:, 5]+dos_atom_2[:, 6]+dos_atom_2[:, 7]+dos_atom_2[:, 8]).reshape(k_num, 1)
    atom2_total = atom2_s + atom2_p + atom2_d
    atom2_dos = np.hstack((energy, atom2_total, atom2_s, atom2_p, atom2_d)) # Energy, atom_total_dos, atom_s, atom_p, atom_d.
    np.savetxt( name_atom_2+"_dos.txt", dos_atom_2)
    for i in range(atom_type_1+atom_type_2, atom_type_1+atom_type_2+atom_type_3):
        dos_atom_33 = dos_split[i]
        dos_atom_3 = dos_atom_33 + dos_atom_3
    atom3_s = (dos_atom_3[:, 0]).reshape(k_num, 1)
    atom3_p = (dos_atom_3[:, 1]+dos_atom_3[:, 2]+dos_atom_3[:, 3]).reshape(k_num, 1)
    atom3_d = (dos_atom_3[:, 4]+dos_atom_3[:, 5]+dos_atom_3[:, 6]+dos_atom_3[:, 7]+dos_atom_3[:, 8]).reshape(k_num, 1)
    atom3_total = atom3_s + atom3_p + atom3_d
    atom3_dos = np.hstack((energy, atom3_total, atom3_s, atom3_p, atom3_d)) # Energy, atom_total_dos, atom_s, atom_p, atom_d.
    np.savetxt( name_atom_3+"_dos.txt", dos_atom_3)
    for i in range(atom_type_1+atom_type_2+atom_type_3, atom_type_1+atom_type_2+atom_type_3+atom_type_4):
        dos_atom_44 = dos_split[i]
        dos_atom_4 = dos_atom_44 + dos_atom_4
    atom4_s = (dos_atom_4[:, 0]).reshape(k_num, 1)
    atom4_p = (dos_atom_4[:, 1]+dos_atom_4[:, 2]+dos_atom_4[:, 3]).reshape(k_num, 1)
    atom4_d = (dos_atom_4[:, 4]+dos_atom_4[:, 5]+dos_atom_4[:, 6]+dos_atom_4[:, 7]+dos_atom_4[:, 8]).reshape(k_num, 1)
    atom4_total = atom4_s + atom4_p + atom4_d
    atom4_dos = np.hstack((energy, atom4_total, atom4_s, atom4_p, atom4_d)) # Energy, atom_total_dos, atom_s, atom_p, atom_d.
    np.savetxt( name_atom_4+"_dos.txt", dos_atom_4)
    
######################################
if atom_type_num==5:              # have 5 types of atom
    name_atom_1 = name_atom.split()[0]
    name_atom_2 = name_atom.split()[1]
    name_atom_3 = name_atom.split()[2]
    name_atom_4 = name_atom.split()[3]
    name_atom_5 = name_atom.split()[4]
    atom_type_1 = int(atom_type[0])
    atom_type_2 = int(atom_type[1])
    atom_type_3 = int(atom_type[2])
    atom_type_4 = int(atom_type[3])
    atom_type_5 = int(atom_type[4])
    for i in range(atom_type_1):
        dos_atom_11 = dos_split[i]
        dos_atom_1 = dos_atom_11 + dos_atom_1
    atom1_s =  (dos_atom_1[:, 0]).reshape(k_num, 1)
    atom1_p = (dos_atom_1[:, 1]+dos_atom_1[:, 2]+dos_atom_1[:, 3]).reshape(k_num, 1)
    atom1_d = (dos_atom_1[:, 4]+dos_atom_1[:, 5]+dos_atom_1[:, 6]+dos_atom_1[:, 7]+dos_atom_1[:, 8]).reshape(k_num, 1)
    atom1_total = atom1_s + atom1_p + atom1_d
    atom1_dos = np.hstack((energy, atom1_total, atom1_s, atom1_p, atom1_d)) # Energy, atom_total_dos, atom_s, atom_p, atom_d.
    np.savetxt( name_atom_1+"_dos.txt", atom1_dos)
    for i in range(atom_type_1, atom_type_1+atom_type_2):
        dos_atom_22 = dos_split[i]
        dos_atom_2 = dos_atom_22 + dos_atom_2
    atom2_s = (dos_atom_2[:, 0]).reshape(k_num, 1)
    atom2_p = (dos_atom_2[:, 1]+dos_atom_2[:, 2]+dos_atom_2[:, 3]).reshape(k_num, 1)
    atom2_d = (dos_atom_2[:, 4]+dos_atom_2[:, 5]+dos_atom_2[:, 6]+dos_atom_2[:, 7]+dos_atom_2[:, 8]).reshape(k_num, 1)
    atom2_total = atom2_s + atom2_p + atom2_d
    atom2_dos = np.hstack((energy, atom2_total, atom2_s, atom2_p, atom2_d)) # Energy, atom_total_dos, atom_s, atom_p, atom_d.
    np.savetxt( name_atom_2+"_dos.txt", dos_atom_2)
    for i in range(atom_type_1+atom_type_2, atom_type_1+atom_type_2+atom_type_3):
        dos_atom_33 = dos_split[i]
        dos_atom_3 = dos_atom_33 + dos_atom_3
    atom3_s = (dos_atom_3[:, 0]).reshape(k_num, 1)
    atom3_p = (dos_atom_3[:, 1]+dos_atom_3[:, 2]+dos_atom_3[:, 3]).reshape(k_num, 1)
    atom3_d = (dos_atom_3[:, 4]+dos_atom_3[:, 5]+dos_atom_3[:, 6]+dos_atom_3[:, 7]+dos_atom_3[:, 8]).reshape(k_num, 1)
    atom3_total = atom3_s + atom3_p + atom3_d
    atom3_dos = np.hstack((energy, atom3_total, atom3_s, atom3_p, atom3_d)) # Energy, atom_total_dos, atom_s, atom_p, atom_d.
    np.savetxt( name_atom_3+"_dos.txt", dos_atom_3)
    for i in range(atom_type_1+atom_type_2+atom_type_3, atom_type_1+atom_type_2+atom_type_3+atom_type_4):
        dos_atom_44 = dos_split[i]
        dos_atom_4 = dos_atom_44 + dos_atom_4
    atom4_s = (dos_atom_4[:, 0]).reshape(k_num, 1)
    atom4_p = (dos_atom_4[:, 1]+dos_atom_4[:, 2]+dos_atom_4[:, 3]).reshape(k_num, 1)
    atom4_d = (dos_atom_4[:, 4]+dos_atom_4[:, 5]+dos_atom_4[:, 6]+dos_atom_4[:, 7]+dos_atom_4[:, 8]).reshape(k_num, 1)
    atom4_total = atom4_s + atom4_p + atom4_d
    atom4_dos = np.hstack((energy, atom4_total, atom4_s, atom4_p, atom4_d)) # Energy, atom_total_dos, atom_s, atom_p, atom_d.
    np.savetxt( name_atom_4+"_dos.txt", dos_atom_4)
    for i in range(atom_type_1+atom_type_2+atom_type_3+atom_type_4, atom_type_1+atom_type_2+atom_type_3+atom_type_4+atom_type_5):
        dos_atom_55 = dos_split[i]
        dos_atom_5 = dos_atom_55 + dos_atom_5
    atom5_s = (dos_atom_5[:, 0]).reshape(k_num, 1)
    atom5_p = (dos_atom_5[:, 1]+dos_atom_5[:, 2]+dos_atom_5[:, 3]).reshape(k_num, 1)
    atom5_d = (dos_atom_5[:, 4]+dos_atom_5[:, 5]+dos_atom_5[:, 6]+dos_atom_5[:, 7]+dos_atom_5[:, 8]).reshape(k_num, 1)
    atom5_total = atom5_s + atom5_p + atom5_d
    atom5_dos = np.hstack((energy, atom5_total, atom5_s, atom5_p, atom5_d)) # Energy, atom_total_dos, atom_s, atom_p, atom_d.
    np.savetxt( name_atom_5+"_dos.txt", dos_atom_5)

   
