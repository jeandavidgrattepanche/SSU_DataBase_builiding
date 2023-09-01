#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "1, September 1, 2023"
__email__ = "jeandavid.grattepanche@gmail.com"



import sys
import os
import re
import urllib.request
import time
import string
import os.path
from Bio import SeqIO
from sys import argv


def sort_cluster(HTSfile,idmin ):
	print("\n\n\n\n\nsort sequence by length\n\n")
	fastalist = []
	fastadict= {}
	for record in SeqIO.parse(HTSfile,'fasta'): 
		ID  = record.id
		IDL  = record.id, int(len(record.seq))
		fastalist.append(IDL)
		fastadict[ID] = record.seq
	print(len(fastalist))
	out = open(HTSfile.split('.fa')[0]+'_sorted.fasta','w+')
	for length in sorted(fastalist, reverse=True, key=lambda x: x[1]):#.sort():
		out.write('>'+length[0] + '\n'+ str(fastadict[length[0]])+'\n')
	out.close()	
# 	print("vsearch --input "+HTSfile.split('.fa')[0]+'_sorted.fasta'+" --uc results"+HTSfile.split('.fa')[0]+".uc --id "+str(idmin))
# 	os.system("vsearch --cluster_smallmem "+HTSfile.split('.fa')[0]+'_sorted.fasta'+" --uc results"+HTSfile.split('.fa')[0]+".uc --id "+str(idmin))

	masterseq = [];mastersequnique =[]; masterseqrided = {}
	input1 = open('results'+HTSfile.split('.fa')[0]+'.uc','r')
	for row in input1:
		if row.split('\t')[0] == 'H':
			cluster = row.split('\t')[1]
			if row.split('\t')[8].split('_')[4] == row.split('\t')[9].split('\n')[0].split('_')[4]:
				print(cluster, "Idem at ", idmin*100 , " : ", row.split('\t')[8], row.split('\t')[9].split('\n')[0])
			else:
				if os.path.is_file('Cluster/cluster_'+cluster+'.fasta'):
					input2 = open('Cluster/cluster_'+cluster+'.fasta','a+')
					input2.write('>'+row.split('\t')[9]+'\n'+str(fastadict[row.split('\t')[9]]) +'\n')
					input2.close()
				else:
					input2 = open('Cluster/cluster_'+cluster+'.fasta','w+')
					input2.write('>'+row.split('\t')[8]+'\n'+str(fastadict[row.split('\t')[8]]) +'\n')	
					input2.write('>'+row.split('\t')[9]+'\n'+str(fastadict[row.split('\t')[9]]) +'\n')	
					input2.close()
	os.system('rm '+ 	HTSfile.split('.fa')[0]+'_sorted.fasta')
def main():
	print("\n\n*******************************************************************")
	print("This script:  \n\t Create cluster by similarity \n")
	print("\n\npython3 Contamination_removal_clustering.py Sequencefile.fas ID% \n\n")
	print("\nFor example\npython3 Contamination_removal_clustering.py EF1_Alpha_MostLKH_Data.fasta 98 \n\n")
	print("*******************************************************************\n\n")
	script,  HTSfile, idminy = argv
	idmin = int(idminy)/100
	sort_cluster(HTSfile,idmin )
main()