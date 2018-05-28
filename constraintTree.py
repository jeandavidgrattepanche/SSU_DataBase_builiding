#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "1, February 7, 2017"
__email__ = "jeandavid.grattepanche@gmail.com"

import string
import re
import sys
import os
from sys import argv
from Bio import SeqIO
from Bio.SeqUtils import GC

seqlist=[];seqlen=[]; seqGC=[]; genuslist= []

def main():
	script, seqfile = argv
	out = open(seqfile.split('.fasta')[0]+ "_Amtaxalist.txt",'w+');out2 = open(seqfile.split('.fasta')[0]+ "_outtaxalist.txt",'w+'); outseq = open(seqfile.split('.fasta')[0]+ "_clean.fasta",'w+')
	am = []; other = []
	for seq in SeqIO.parse(open(seqfile,'r'),'fasta'):
		if seq.seq in seqlist:
			print(seq.id)
		else:
			seqlist.append(seq.seq)
#			print (seq.seq)
			outseq.write(">" + seq.id + '\n' + str(seq.seq) + '\n')
			if "Am_" in seq.id:
				out.write(seq.id.replace(';size=','_')+',')
				amtaxa = seq.id.replace(';size=','_').replace('/','').replace(':','') #+':0'
				am.append(amtaxa)
#				print("okay")
			
			else:
				print("out", seq.id)
				out2.write(seq.id.replace(';size=','_')+',')
				othertaxa = seq.id.replace(';size=','_').replace('/','').replace(':','') #+':0'
				other.append(othertaxa)
	outfinaltree = open('constrainttree.tre','w+')
	outfinaltree.write('('+str(am).replace('[','').replace(']','').replace("'","")+',('+str(other).replace('[','').replace(']','').replace("'","")+'))')
main()