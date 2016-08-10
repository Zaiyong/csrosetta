from pymol import cmd
import pymol

def rename_atoms(a1, residues, atom_names):
	aa1 = a1
	if a1 == "QD2" or a1 == "HD2":
		aa1 = "HD21"
	if a1 == "QQD" :
		aa1 = "HD11"
	if a1 == "QQG" :
		aa1 = "HG11"
	if a1 == "QE2" :
		aa1 = "HE21"
	if a1 == "CEN" :
		aa1 = "CB"
	aa1 = aa1.replace("Q", "C")
	if (aa1 == "HG13" or aa1 == "3HG1") and residues[-1] == "ILE" and not aa1 in atom_names:
		aa1 = "1HG1"
	if ( (aa1 == "HB3" or aa1 == "3HB") and not residues[-1] == "ALA" and not aa1 in atom_names):
		aa1 = "HB1"
	if ( (aa1 == "QE2" or aa1 == "HE2") and not residues[-1] == "GLN"):
		aa1 = "1HE2"
	return aa1


def show_restraints( restraint_file , sele = "all", show_sidechains = True, print_label = False ):

	tp = [0,0,0]
	fp = [0,0,0]

	restraint_count = 0
	restraint_count_shown = 0


	for line in open(restraint_file,'r'):

		if line.strip() == "":
			continue

		res = line.split()

		restraint_count += 1
		label = "_dist_" + str(restraint_count)

		myspace = {'atom_names_1': [],'atom_names_2': [], 'resn1': [], 'resn2' : []}
		cmd.iterate( sele + " and resi " + res[2], 'atom_names_1.append(name)', space=myspace)
		cmd.iterate( sele + " and resi " + res[4], 'atom_names_2.append(name)', space=myspace)
		cmd.iterate( sele + " and resi " + res[2], 'resn1.append(resn)', space=myspace)
		cmd.iterate( sele + " and resi " + res[4], 'resn2.append(resn)', space=myspace)

		aa1 = rename_atoms( res[1], myspace[ 'resn1' ], myspace[ 'atom_names_1' ] )
		aa2 = rename_atoms( res[3], myspace[ 'resn2' ], myspace[ 'atom_names_2' ] )


		if not aa1 in myspace["atom_names_1"]:
			temp = aa1[-1]+aa1[:-1]
			if temp in  myspace["atom_names_1"]:
				aa1 = temp
			else:
				print "Problem finding %s/%s" %(res[2],aa1)
				continue

		if not aa2 in  myspace["atom_names_2"]:
			temp = aa2[-1]+aa2[:-1]
			if temp in  myspace["atom_names_2"]:
				aa2 = temp
			else:
				print "Problem finding %s/%s" %(res[4],aa2)
				continue

		dist = float(cmd.distance(label,"(/%s///%s/%s)"%(sele, res[2], aa1), "(/%s///%s/%s)"%(sele, res[4], aa2)))

		if show_sidechains == True:
			# get colors of sidechains and color accordingly
			pymol.color_list = []
			cmd.iterate("/"+sele+"///"+res[2]+"/"+aa1, 'pymol.color_list.append(color)')
			cmd.show("lines", "%s and resi %s" %(sele,res[2]) )
			cmd.color(pymol.color_list[0], "%s and resi %s" %(sele,res[2]) )
			pymol.color_list = []
			cmd.iterate("/" + sele+"///"+res[4]+"/"+aa2, 'pymol.color_list.append(color)')
			cmd.show("lines", "%s and resi %s" %(sele,res[4]) )
			cmd.color(pymol.color_list[0], "%s and resi %s" %(sele,res[4]) )

		#color restraints based on actual distance
		lb = 0
		ub = 100
		if res[5].lower() == "sigmoid":
			ub = float(res[6])
			lb = 0.0
		elif res[5].lower() == "bounded":
			ub = float(res[7])
			lb = float(res[6])
		else:
			print "Function not recognized. Setting lower and upperbound to 0 and 100 ... "


		ff = ""

		if dist <= ub and dist >= lb:
			color = "blue"
			tp[0] += 1
			tp[1] += 1
			tp[2] += 1
			ff = "SATISFIED"
		elif dist <= ub + 0.5 and dist >= lb - 0.5:
			color = "purpleblue"
			fp[0] += 1
			tp[1] += 1
			tp[2] += 1
		elif dist <= ub + 0.5 and dist >= lb - 0.5:
			color = "magenta"
			tp[2] += 1
			fp[0] += 1
			fp[1] += 1
		else:
			color = "red"
			fp[0] += 1
			fp[1] += 1
			fp[2] += 1

		print "%s/%s %s/%s with distance %.2f   Restraint: %.2f - %.2f     %s"%(res[2], aa1, res[4], aa2, dist, lb, ub, ff)

		restraint_count_shown += 1
		cmd.color(color, label)
		if not print_label:
			cmd.hide("label",label)
		cmd.set('dash_gap','0')
	print "Total Restraints: %.0f"%restraint_count
	print "Shown Restraints: %.0f"%restraint_count_shown
	print "Padding: [0, 0.5, 1]"
	print "Satisfied Restraints: " + str(tp)
	print "Unsatisfied Restraints: " + str(fp)

cmd.extend("showRestraints",show_restraints)
