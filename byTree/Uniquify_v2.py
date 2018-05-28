#!/usr/bin/python2

__author__ = "Jean-David Grattepanche"
__version__ = "7, August 30, 2016"
__email__ = "jeandavid.grattepanche@gmail.com"

'''
This script will remove closely related sequences based on the RID percent you choose.

'''


from Bio import Entrez
Entrez.email = 'jgrattepanche@smith.edu'
import string
import re
import sys
import os
import os.path
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from datetime import date

recordlist= []
chars_to_remove = ["(", ")" , ":" , ";" , "[", "]", "." , "," ,]

def Rename(inputfile):
	print "start Adding Name"
	for record in SeqIO.parse(open(inputfile,'r'),'fasta'): 
		ID2  = record.id.split('|')[3] #.split('>')[1]
		print(ID2)
		handle2 = Entrez.efetch(db="nucleotide", id=ID2, rettype="gb")
		record2 = SeqIO.read(handle2,"gb")
		print(record2.annotations["taxonomy"],record2.annotations["organism"]) #and float(alignment.hsps[0].identities)/float(alignment.hsps[0].align_length) > 90: ??
		try:
			newname = record2.annotations["taxonomy"][1][0:3] + "_" +record2.annotations["taxonomy"][2][0:2] + "_" +record2.annotations["taxonomy"][3][0:2] + "_" +  record2.annotations["organism"].replace(' ','_') + '_' + ID2 # need to be modified related to the taxonomy list =>record2.annotations["taxonomy"][14][0:3]
		except:
			try: 
				newname = record2.annotations["taxonomy"][1][0:3] + "_" +record2.annotations["taxonomy"][2][0:2]  + "__" +  record2.annotations["organism"].replace(' ','_') + '_' + ID2 # need to be modified related to the taxonomy list =>record2.annotations["taxonomy"][14][0:3]
			except:
				try: 
					newname = record2.annotations["taxonomy"][1][0:3] + "___" +  record2.annotations["organism"].replace(' ','_') + '_' + ID2 # need to be modified related to the taxonomy list =>record2.annotations["taxonomy"][14][0:3]
				except:
					newname = "SAR___" + record2.annotations["organism"].replace(' ','_') + '_' + ID2
		print(newname)
		rename2 = newname.translate(None, ''.join(chars_to_remove))
		print record.description, " newname => " , rename2
		outrename = open(inputfile.split('.')[0] + '_rename.fas','a')
		outrename.write('>' + rename2 + '\n' + str(record.seq) + '\n')
		outrename.close()
		
def Addnumber(inputfile):
	print "start Adding Number"
	for record in SeqIO.parse(open(inputfile.split(".")[0] +'_RiddedSeqs.fasta','r'),'fasta'):
		ID2  = record.id
		whor = open(inputfile.split(".")[0] +'_WhoRiddedWho.txt','r')
		riddedp = ""
		for ridded in whor:
			if record.id == riddedp.split(":")[0]:
#				print ridded
				numberridded = ridded.count('>')
				name_num = riddedp.split(":")[0] + "_rid_" + str(numberridded)
				print name_num
				outnamnum = open(inputfile.split(".")[0] +"_uniquified.fas","a")
				outnamnum.write(">" + name_num + "\n" + str(record.seq) + '\n' )
				outnamnum.close()
			riddedp = ridded


def main():
	print "This script will remove closely related sequences based on the RID percent you choose. \n"
	inputfile = raw_input('What is the file of gene sequences?  ')
	i = raw_input('What percent would you like to RID at? (hit return for default of 99%) ')
	try:
		num = float(i) + 1
	except TypeError:
		print 'Your input must be a number.  Try again. '
		main()
	except ValueError:
		i = ""	
	if i == "":
		RIDin = 99.00
	else:
		RIDin = float(i)
	i = raw_input('What is the minimum identity? (hit return for default of 60%) ')
	try:
		num = float(i) + 1
	except TypeError:
		print 'Your input must be a number.  Try again. '
		main()
	except ValueError:
		i = ""	
	if i == "":
		IDmin = 60.00
	else:
		IDmin = float(i)
	z = raw_input('Hit return when you are ready to continue. ')
	
#	Rename(inputfile)
	
	import rid_by_tree_water_v2
	rid_by_tree_water_v2.rid(RIDin, inputfile, IDmin)
	
	Addnumber(inputfile)	
	
#	os.system('mafft --auto --reorder --anysymbol --thread 4 ' + inputfile.split(".")[0] +'_uniquified.fas  >  ' + inputfile.split(".")[0] +'_uniquified_aligned.fas')
#	os.system('FastTree -nt ' + inputfile.split(".")[0] +'_uniquified_aligned.fas  >  ' + inputfile.split(".")[0] + '_uniquified.tree')
#	os.system('raxmlHPC-PTHREADS-SSE3 -m GTRGAMMAI -s ' + inputfile.split(".")[0] +'_uniquified_aligned.fas -n uniquified_tree_bootstrap -T 4 -p 12345 -x 12345 -# autoMRE -f a')
	sys.stderr.close()
main()