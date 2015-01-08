#!/usr/bin/env python
import argparse

#___ADD_PARSER_ARGUMENTS_______________________________________________________________________________________
parser = argparse.ArgumentParser(description="Simulate molecular dynamics of a structure")
parser.add_argument("--iF", type=str, help="Path to the input file.")
parser.add_argument("--kB", type=float, default = 40000.0, help="Same as kb.")
parser.add_argument("--kN", type=float, default = 400.0, help="Same as kn.")

parser.add_argument("--nbCutoff", type=float, default = 0.50, help="Distance within which atoms should be considered as having nonbonded interactions, if they are not a bonded pair")
parser.add_argument("--m", type=float, default = 12.0, help="Atom mass, a constant to be applied to all atoms.")
parser.add_argument("--dt", type=float, default = 0.001, help="Length of time step.")
parser.add_argument("--n", type=int, default = 1000, help="Number of time steps to iterate, useful for debugging.")
parser.add_argument("--out", type=str, help="Prefix of the filename to output.")

args = parser.parse_args()
if not args.iF:
	exit()
if not args.out:
	args.out = args.iF.split('.')[0] + "_out"


#___FUNCTION: Initialize__________________________________
#Reads in data from input files and parses the header values
#Creates the coordinate, velocity, and bond connection matrices (that are passed in as parameters)
#Opens the output files and writes the headers
#INPUT: empty coordinate, velocity, and bonded-atoms matrix
#OUTPUT: the erg and tvc output files
def initialize(cm,vm,bm):
	out_erg = open(args.out + ".erg",'w+')
	out_rvc = open(args.out + ".rvc",'w+')
	out_erg.write("#\tstep\tE_k\tE_b\tE_nB\tE_tot\n")
	lines = file(args.iF).readlines()
	out_rvc.write(lines[0])
	for line in lines[1:]:
		out_rvc.write(line)
		arr = line.split()
		floats = [float(n) for n in arr[1:7]]
		cm.append(floats[:3])
		vm.append(floats[3:])
		bm.append([int(n)-1 for n in arr[7:]])
	return [out_erg,out_rvc]

#___FUNCTION: INITIALiZE_MATRIX__________________________________
#Creates a two dimensional matrix of the specified value
#INPUT: number of rows, columns, and the values to fill the matrix with
#OUTPUT: desired matrix
def initialize_matrix(num_rows,num_cols,value):
	matrix = []
	for i in range(num_rows):
		row = []
		for j in range(num_cols):
			row.append(value)
		matrix.append(row)
	return matrix

#___FUNCTION: RESET_MATRIX__________________________________
#Resets all values of a two dimensional matrix to the specified value
#INPUT: the matrix to reset, the number of rows, number of columns, and the value to reset the matrix with
#OUTPUT: none
def reset_matrix(matrix,x,y,value):
	for i in range(x):
		for j in range(y):
			matrix[i][j] = value

#___FUNCTION: EUCLIDEAN_____________________________________
#INPUT: two lists floats that are of the same length
#OUTPUT: euclidean distance between the two lists
def euclidean(x,y):
	if len(x) != len(y):
		return -1
	sumSq = 0.0
	for i in range(len(x)):
		sumSq+=(x[i]-y[i])**2
	return (sumSq**0.5)

#___FUNCTION: PDIST_____________________________________
#INPUT: Coordinate matrix of atoms
#OUTPUT: matrix of (euclidean) distances between every pair of atoms (pdist[i,j] = pdist[j,i])
def pdist(cm):
	pdist = initialize_matrix(len(cm),len(cm),0.0)
	for i in range(len(cm)):
		for j in range(len(cm)):
			if i < j:
				pdist[i][j] = euclidean(cm[i],cm[j])
			elif j == i:
				pdist[i][j] = 0.0
			else:
				pdist[i][j] = pdist[j][i]
	return pdist

#___FUNCTION: FIND_NB_ATOMS_____________________________________
#INPUT: reference (initial) distnace matrix and the nested list of bonded atoms
#OUTPUT: nested list of all bonded atoms(nbm) using the reference distance matrix (rdist)
#Each index of the returned list contains a list of all atoms that are within the nbCutoff distance from the indexed atom
def find_nb_atoms(rdist,bm):
	nbm = []
	for i in range(len(rdist)):
		row = []
		for j in range(len(rdist)):
			if (rdist[i][j] < args.nbCutoff) and (i != j) and (j not in bm[i]):
				row.append(j)
		nbm.append(row)
	return nbm

#___FUNCTION: UPDATE_CM_____________________________________
#Updates coordinate matrix of atoms(cm), using the half time velocity matrix(vm)
#INPUT: the (half-time)velocity and coordinate matrices of the atoms
#OUTPUT: none
def update_cm(vm,cm):
	for i in range(len(cm)):
		for k in range(3):
			cm[i][k] += vm[i][k]*args.dt

#___FUNCTION: UPDATE_FM_____________________________________
#Updates force matrix of atoms(fm), using the reference distance matrix(rdist), current distance matrix(dist), nested list of bonded atoms(bm), nested list of non-bonded atoms(nbm), and current coordinate matrix(cm)
#In addition, it calculates and returns the bonded and non-bonded potential energy values in a list
#INPUT: the reference distance matrix, bonded atoms list, non-bonded atoms list, (updated)coordinate matrix, and current force matrix 
#OUTPUT: the accumulated bonded potential energy andd non-bonded potential energy of the forces
def update_fm(rdist,bm,nbm,cm,fm):
	bEnergy = 0.0
	nbEnergy = 0.0
	for i in range(len(cm)):
		for atomIndex in bm[i]: 
			if atomIndex > i:
				b = euclidean(cm[i],cm[atomIndex])
				diff = b-rdist[i][atomIndex]
				bEnergy += args.kB*(diff**2)/2
				F = args.kB*diff #magnitude of force
				for k in range(3):
					f = F*(cm[atomIndex][k]-cm[i][k])/b
					fm[i][k] += f
					fm[atomIndex][k] -= f
		#add all non-bonding forces
		for atomIndex in nbm[i]:
			if atomIndex > i:
				b = euclidean(cm[i],cm[atomIndex])
				diff = b-rdist[i][atomIndex]
				nbEnergy += args.kN*(diff**2)/2
				F = args.kN*diff #magnitude of force
				for k in range(3):
					f = F*(cm[atomIndex][k]-cm[i][k])/b
					fm[i][k] += f
					fm[atomIndex][k] -= f
	return [bEnergy,nbEnergy]

#___FUNCTION: UPDATE_VM_____________________________________
#Updates velocity matrix of atoms(vm), using the force matrix(fm)
#Note that this function is used to find both half-time and regular velocities: the differentiating aspect is the force matrix passed in
#In addition, it calculates and returns the kinetic energy value
#INPUT: the updated force matrix and (half-time)velocity matrix
#OUTPUT: the calculated sum of kinetic energies of each atom
def update_vm(fm,vm):
	kEnergy = 0.0
	for i in range(len(vm)):
		for k in range(3):
			a = fm[i][k]/args.m
			vm[i][k] += a*args.dt/2
			kEnergy += args.m*(vm[i][k])**2/2
	return kEnergy

#___FUNCTION: WRITE_TO_ERG______________________________________
#Write the energy information from the time step to the erg file
#INPUT: the step number, energy values, and out file
#OUTPUT: none
def write_to_erg(step,kEnergy,bEnergy,nbEnergy,tEnergy,out):
	out.write(str(step)+'\t'+str(round(kEnergy,1))+'\t'+str(round(bEnergy,1))+'\t'+str(round(nbEnergy,1))+'\t'+str(round(tEnergy,1))+'\n')

#___FUNCTION: WRITE_TO_RVC______________________________________
#Write the atom informations from the time step to the rvc file
#INPUT: the step number, coordinate and velocity matrix, total energy, bonded atoms list, and out file
#OUTPUT: none
def write_to_rvc(step,cm,vm,tEnergy,bm,out):
	out.write("#At time step "+str(step)+",energy = "+str(round(tEnergy,3))+"kJ\n")
	for i in range(len(cm)):
		b_atoms = ""
		for atom in bm[i]:
			b_atoms += '\t' + str(atom+1)
		out.write(str(i+1)+'\t'+str(round(cm[i][0],4))+'\t'+str(round(cm[i][1],4))+'\t'+str(round(cm[i][2],4))+'\t'+str(round(vm[i][0],4))+'\t'+str(round(vm[i][1],4))+'\t'+str(round(vm[i][2],4))+b_atoms+'\n')

#___FUNCTION: WRITE_TO_EUC______________________________________
#Write the distance between the binding sites (manually specified) from this time step to the euc file
#INPUT: step number, coordinate matrix, outfile
#OUTPUT: none
def write_to_euc(step,cm,out):
	line = str(step)
	line += '\t' + str(round(euclidean(cm[17],cm[96]),4))
	out.write(line + '\n')

#___FUNCTION: MAIN_____________________________________________________________________________________________
#Performs the main function of calculating position and velocity of the input atoms as well as their kinetic, potential, and total energies at each time step.
#It also writes this information to two files (_out.rvc and _out.erg) according to the required specifications.
#Note: the commented out functions are optional for writing information to an euc file
def main():
	cm = []
	vm = []
	bm = []
	[out_erg,out_rvc] = initialize(cm,vm,bm)
	rdist = pdist(cm) #initial matrix of pair-wise distances, used for REFERENCE
	nbm = find_nb_atoms(rdist,bm) #non-bond matrix --> atom IDS for nonbonding atoms within the nbcutoff
	fm = initialize_matrix(len(cm),3,0.0)
	dist = initialize_matrix(len(cm),len(cm),0.0)
	energy0 = update_vm(fm,vm) #calculates initial energy (with zero matrix)
	# out_euc = open(args.out + ".euc",'w+')
	# write_to_euc(0,cm,out_euc)

	for n in range(args.n):
		update_vm(fm,vm) #Calculates half time velocities in vm
		update_cm(vm,cm) #Updates cm according to half time velocities
		reset_matrix(fm,len(fm),3,0.0)
		[bEnergy,nbEnergy] = update_fm(rdist,bm,nbm,cm,fm) #Updates force matrix
		kEnergy = update_vm(fm,vm) #calculates velocities based on updated force and half time velocities
		tEnergy = bEnergy + nbEnergy + kEnergy
		if tEnergy > 10*energy0:
			print "Energy exceeded 10x initial energy on step " + str(n+1) + ": Molecule Unstable"
			break
		if (n+1)%10 == 0:
			write_to_erg(n+1,kEnergy,bEnergy,nbEnergy,tEnergy,out_erg)
			write_to_rvc(n+1,cm,vm,tEnergy,bm,out_rvc)
			# write_to_euc(n+1,cm,out_euc)

main()
