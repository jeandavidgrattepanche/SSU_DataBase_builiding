#!/usr/bin/python3
#from Bio import Entrez
#Entrez.email = 'jgrattepanche@smith.edu'
import string
import re, sys, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio import Align
from Bio.Align import MultipleSeqAlignment

def maskalignment(arg, percent,filetype):
	maskFileName = arg + '_masked_' + str(percent) + '.fas'
	outFile = open(maskFileName,'a')
	alignment = AlignIO.read(arg, filetype)
	trimAlign = MultipleSeqAlignment([])
	numRows = len(alignment)
	x = float(percent) * float(numRows) / 100.0
	numGap = numRows - int(x)
	numCol = alignment.get_alignment_length()

	print("Total number of rows: %i" % numRows)
	print("Number of gapped sequences allowed at a given site: %i" % numGap)
	print("Total number of columns: %i" % numCol)
	my_array = {}
	colToKeep=[]
	for i in range(numCol):
		#print (i)
		lineName = "line_" + str(i)
		my_array[lineName] = alignment[:,i]
		if my_array[lineName].count('-') > numGap:
			print ("get rid of column %i" % i)
		else:
			colToKeep.append(i)
	
	for record in alignment:
		newseq = ""
		for i in colToKeep:
			newseq= newseq + (record[i])
			
		newRecord = SeqRecord(Seq(newseq), id=record.id)
		trimAlign.append(newRecord)
		outFile.write('>' + record.id + '\n' + newseq + '\n')
		
	print("Total number of columns remaining: %i" % trimAlign.get_alignment_length())



def writelog(input):
	outfile = open('log.txt','a')
	outfile.write(input)
	outfile.close


def main():
	print ("*************************************************************************************************")
	print ("This script will take an alignment and return an alignment with the gapped columns removed." )
	print ("Useage is 'python mask_gaps.py <inputAlignment>'")
	print ("*************************************************************************************************\n\n")
	
	y = input('What type of file is this? (fasta, nexus, phyllip) ')
	if y[0] == 'n':
		filetype = 'nexus'
	elif y[0] == 'f':
		filetype = 'fasta'	
	elif y[0] == 'p':
		filetype = 'phylip'
	else:
		print('that is not a valid file type.  try again')
		main()
	
	print ('At what percent do you want to mask? (i.e. at 75% a column will be removed if fewer than 3/4 the sequences have a character at that position) ')
	x = input('Hit enter for the default of 75%) ')
	try:
		num = float(x) + 1
	except TypeError:
		print ('Your input must be a number.  Try again. ')
		main()
	except ValueError:
		x = ""
	if x == "":	
		percent = 75.0
	else:
		percent = float(x)

	for arg in sys.argv[1:]:
		maskalignment(arg, percent,filetype)


main()
