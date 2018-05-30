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
	maskedcolumn =  open(arg.split('.')[0] + '_mask_' + str(percent) + '.txt', 'w+')
	outFile = open(arg + '_masked_' + str(percent) + '.fas','w+')
	checkgap = open(arg.split('.')[0] + '_missingcharacter.txt', 'w+')
	alignment = AlignIO.read(arg, filetype)
	trimAlign = MultipleSeqAlignment([])
	numRows = len(alignment)
	x = float(percent) * float(numRows) / 100.0
	numGap = numRows - float(x)
	numCol = alignment.get_alignment_length()

	print("Total number of rows: %i" % numRows)
	print("Number of gapped sequences allowed at a given site: %i" % numGap)
	print("Total number of columns: %i" % numCol)
	checkgap.write("Total number of rows: \t"+ str(numRows) + '\nNumber of gapped sequences allowed at a given site: \t'+ str(numGap) +'\n Total number of columns: \t'+ str(numCol)+ '\n\n cutoff : \t'+ str(x)+ '\n\n\n')
	checkgap.write("Position \t Missing Characters \t Characters \n")
	my_array = {}
	colToKeep=[]
	for i in range(numCol):
		#print i
		lineName = "line_" + str(i)
		my_array[lineName] = alignment[:,i]
		chapre = int(numRows)-int(my_array[lineName].count('-'))
		checkgap.write(str(i)+'\t'+ str(my_array[lineName].count('-'))+'\t'+ str(chapre)+'\n')
		if my_array[lineName].count('-') > numGap:
			print("get rid of column %i" % i)
			maskedcolumn.write(str(i) +'\n')
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
	print("*************************************************************************************************")
	print("This script will take an alignment and return an alignment with the gapped columns removed." )
	print("Useage is 'python mask_gaps.py <inputAlignment>'")
	print("*************************************************************************************************\n\n")
	
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
	
	print('At what percent do you want to mask? (i.e. at 75% a column will be removed if fewer than 1/4 the sequences have a character at that position i.e. 75% of missing data) ')
	x = input('Hit enter for the default of 75%) ')
	try:
		num = float(x) + 1
	except TypeError:
		print('Your input must be a number.  Try again. ')
		main()
	except ValueError:
		x = ""
	if x == "":	
		percent = 25.0
	else:
		percent = 100.00-float(x)

	for arg in sys.argv[1:]:
		maskalignment(arg, percent,filetype)


main()
