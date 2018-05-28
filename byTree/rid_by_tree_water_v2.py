#!/usr/bin/python2

__author__ = "Jean-David Grattepanche"
__version__ = "6, April 21, 2016"
__email__ = "jeandavid.grattepanche@gmail.com"


import re,os
from Bio import SeqIO
from Bio import Phylo
toKeep = []
toRemove = []
seqDict = {}

def rid(RIDin, inputfile, IDmin):
	infile=SeqIO.parse(open(inputfile,'r'),'fasta')
	for seq in infile:
		try:
			name = seq.id #.replace('.','_')
#			print name, "dict"
		except:
			name = seq.id #.replace('.','_')	
		seqDict[name] = seq.seq
	


#build a mafft guide tree.
		
	os.system('mafft --auto --reorder --anysymbol --thread 4 ' +   inputfile +'  >  '+ inputfile.split('.')[0] + '_mafftout') # --adjustdirection can be added if not sure of the direction of the sequences
	os.system('FastTree -nt ' + inputfile.split('.')[0] + '_mafftout >  '+ inputfile.split('.')[0] + '_FastTree.tree') # model can be added (-gtr -gamma )
	treefile = inputfile.split('.')[0] + '_FastTree.tree'
#	print treefile
	tree = Phylo.read(treefile,'newick')
	checkClade(tree, RIDin,inputfile,IDmin)

#take a leaf, compare it to its sister
def get_parent(tree, child_clade):
		#############################################################
		#############################################################
	node_path = tree.get_path(child_clade)
	try:
		return node_path[-2]				
	except:
		return 'None' 	
	
def checkWater(tree, parent, child, keepone, seen, RIDin, IDmin):
	answer = []
	f = str(tree.name)
	for leaf in parent.get_terminals():
		#newname = 
		out = open('seq1.fasta','w')
#		print leaf.name, "treename"
		if leaf.name[0:3] == "_R_":
			leafrename = leaf.name[3:] #.split('_')[1]+'_'+leaf.name.split('_')[2]+'.'+leaf.name.split('_')[3]
		else:
			try:
				leafrename = leaf.name #.split('_')[1]+'.'+leaf.name.split('_')[2]
			except:
				leafrename = leaf.name #.split('_')[1]
#		print leafrename, "rename"
		try:
			out.write('>' + leafrename  + '\n' + str(seqDict[leafrename]) + '\n')
		except:
			out.write('>' + leafrename  + '\n' + str(seqDict[leafrename]) + '\n')
                
		out.close()

	for j in range(len(parent.get_terminals())):
		
		if parent.get_terminals()[j] != leaf:
			out2 = open('seq2.fasta','w')			
			if parent.get_terminals()[j].name[0:3] == "_R_":
				parentrename = parent.get_terminals()[j].name[3:] #.split('_')[1]+'_'+parent.get_terminals()[j].name.split('_')[2]+'.'+parent.get_terminals()[j].name.split('_')[3]
			else:
				try:
					parentrename = parent.get_terminals()[j].name #.split('_')[1]+'.'+parent.get_terminals()[j].name.split('_')[2]
				except:
					parentrename = parent.get_terminals()[j].name #.split('_')[1]
#			print parentrename
			try:
				out2.write('>' + parentrename  + '\n' + str(seqDict[parentrename]) + '\n')
			except:
				out2.write('>' + parentrename  + '\n' + str(seqDict[parentrename]) + '\n')
				
			out2.close()
#	
			cline = 'water -outfile=' + f + 'water.txt -asequence=seq1.fasta -bsequence=seq2.fasta -gapopen=10 -gapextend=0.5' #WaterCommandline(asequence="../Temp/seq1.fasta", bsequence="../Temp/seq2.fasta", gapopen=10, gapextend=0.5, outfile=('../Temp/' + f + "water.txt"))		
			os.system(cline)
		
			os.system('cat ' + f + 'water.txt >> allwaterout.txt' )
		
			import WaterRids_v2
			x = WaterRids_v2.test(f, RIDin, IDmin)			
#			print 'x = ' + str(x)
			answer.append(x)
			
			if x == True:
				break
					
	#if all(x in cladeList for x in TestList)	
	if all(ans == False for ans in answer): #all are within cut off, o up a node and try again
		newparent = get_parent(tree, parent)
		if newparent != 'None':
			checkWater(tree, newparent, parent, keepone, seen, RIDin, IDmin)
		else:
			keepone.append(child) #the parent has a node outside the cut off, but the child does not.  Keep one from child, 'rid' the rest.
			for leaf in child.get_terminals():
				seen.append(leaf)
	else:
		keepone.append(child) #the parent has a node outside the cut off, but the child does not.  Keep one from child, 'rid' the rest.
		for leaf in child.get_terminals():
			seen.append(leaf)
	return seen
		
			
def checkClade(tree, RIDin,inputfile, IDmin):
	keepone = []
	seen = []
	written = []
	for clade in tree.get_terminals():
#		print clade.name, "checkclade"
		if clade not in seen:
			#seen.append(clade)
			parent = get_parent(tree,clade)

			if parent == 'None':
				seqlist =  clade.name #.split('_')
				seq1name = ""	
				for i in range(len(seqlist)-1):
					if i != 0:
						seq1name = seq1name + '_' + seqlist[i]
				toKeep.append(seq1name)	
				
			else:			

				seen = checkWater(tree, parent, clade, keepone, seen,RIDin, IDmin)
				
	for clade in keepone:
		#print clade
		out = open(inputfile.split(".")[0] +'_RiddedSeqs.fasta','a')
		out2 = open(inputfile.split(".")[0] +'_WhoRiddedWho.txt','a')
		out2.write('\n' + str(clade.get_terminals()[0].name) + ':\n')
		if clade.get_terminals()[0].name not in written:
			written.append(clade.get_terminals()[0].name)
			if clade.get_terminals()[0].name[0:3] == "_R_":
				claderename = clade.get_terminals()[0].name[3:] #.split('_')[1]+'_'+clade.get_terminals()[0].name.split('_')[2]+'.'+clade.get_terminals()[0].name.split('_')[3]
			else:
				try:
					claderename = clade.get_terminals()[0].name #.split('_')[1]+'.'+clade.get_terminals()[0].name.split('_')[2]
				except:
					claderename = clade.get_terminals()[0].name #.split('_')[1]
#			print claderename
			out.write('>' + claderename  + '\n' + str(seqDict[claderename]) + '\n')
		out.close()
		
		
		for seq in clade.get_terminals():
			if seq != clade.get_terminals()[0]:
				out2.write('>' + str(seq) + ',')
		out2.close()
		
		
