#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "2, September 14, 2016"
__email__ = "jeandavid.grattepanche@gmail.com"



import sys
import os
import re
import urllib.request
import time
import string
import os.path
from Bio import SeqIO
from Bio import Entrez
Entrez.email = 'jgrattepanche@smith.edu'
from sys import argv
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW


def blastonline(NGSfile, idmin, Emin):
	if not os.path.exists('temp/'):
		os.makedirs('temp/')
	if not os.path.exists('BLAST_results/'):
		os.makedirs('BLAST_results/')
	print("start BLAST online eukaryote morphospecies")
	fastalist = []
	fastadict= {}
	for record in SeqIO.parse(NGSfile,'fasta'): 
		ID  = record.id
		fastalist.append(ID)
		fastadict[ID] = record.seq
	for j in range(0, len(fastalist), 100):
		fastalist2 = fastalist[j:j+100]
		print("Create batch :", j, " - ", j+100, "on", len(fastalist), "for BLASTing against morphospecies online.", len(fastalist2), "sequences per batch" )
		for ID2 in fastalist2:
			tempseq = open("temp/seqsbatch"+str(j)+".fas",'a')
			tempseq.write('>' + ID2 + '\n' + str(fastadict[ID2]) + '\n')
			tempseq.close()
		entrezQuery = '("eukaryotae"[organism] AND (all [filter] NOT(environmental samples[filter] OR metagenomes[orgn])))'
		fastafile = open("temp/seqsbatch"+str(j)+".fas").read()
		#for Seq in SeqIO.parse(NGSfile + "_noBlast_result.txt",'fasta'):
		#	print(Seq.id, " added to the blast online")
		#	print("Start Blast online")
		result_handle = NCBIWWW.qblast("blastn", "nr", fastafile, ncbi_gi=False, format_type="XML", entrez_query= entrezQuery, expect = Emin, perc_ident=idmin, hitlist_size=1)
		blast_records = NCBIXML.parse(result_handle)
		for blast_record in blast_records:
			if blast_record.descriptions:
				for i in range(1):
					try:
#						print (i)
						ident =  blast_record.alignments[i].hsps[0].identities
						match = blast_record.alignments[i].hsps[0].match
						Matchlist = []
						cnt = 0
						for Match in match:
							Matchlist.append('x')
						cnt = Matchlist.count('x') 
						evalue = blast_record.alignments[i].hsps[0].expect
						Sim = round(float(ident) / float(cnt) * 100)
						ID = blast_record.alignments[i].title.split('|')[4]
						gb = blast_record.alignments[i].title.split('|')[3]
						print(blast_record.query, 'blasted with', ID , " at ", Sim, "%", "E-value:", evalue)
						out1 = open("BLAST_results/Blast_result.txt",'a')
						out1.write(blast_record.query + '\t' + ID.replace('_','\t')  +'\t' + gb + '\t'+str(Sim) + '\t'+ str(evalue) +'\n')
						out1.close()
					except:
						print("No more blast for", blast_record.query)
			else:
				print("NO BLAST for ", blast_record.query)
				for record in SeqIO.parse(NGSfile,'fasta'): 
					if record.id == blast_record.query:
						outnb2 = open("BLAST_results/noBlast_online.txt",'a')
						outnb2.write('>' + blast_record.query + '\n' + str(record.seq) +'\n')
						outnb2.close()



def Rename(inputfile):
	print("start Adding Name")
	OTUIDlist= []
	OTUIDlist2= []
	IDlist = []
	chars_to_remove = ["(", ")" , ":" , ";" , "[", "]", "," , "/"]
	blastrecord = {}
	for record in open("BLAST_results/Blast_result.txt",'r'): 
		gbID  = record.split('\t')[-3].split('.')[0]
		OTUID = record.split('\t')[0]
		blastrecord[OTUID] = [gbID, record.split('\t')[-2]+"_"+record.split('\n')[0].split('\t')[-1]]
		OTUIDlist.append(OTUID)
		if gbID not in IDlist:
			IDlist.append(gbID)
	for k in range(0, len(IDlist), 200):
		list = IDlist[k:k+200]
		print(k, " - ", k+200, "on", len(IDlist), "efetch batch of:", len(list), "sequences" )
		list2 = str(list).replace('[','').replace(']','').replace("'","").replace(" ","")
		print("efetch:", list2)
		try:
			time.sleep(10)
			handle = Entrez.efetch(db="nucleotide", id=list2, rettype="gb" )
			records = SeqIO.parse(handle,"gb")
		except:
			time.sleep(120)
			try:
				print("try 2")
				handle = Entrez.efetch(db="nucleotide", id=list2, rettype="gb" )
				records = SeqIO.parse(handle,"gb")
			except urllib.request.HTTPError as err:
				if err.code == 502:
					time.sleep(60)
					print("retry efetch")
					raise
		for record2 in records:
#			print(record2.annotations["taxonomy"],record2.annotations["organism"]) #and float(alignment.hsps[0].identities)/float(alignment.hsps[0].align_length) > 90: ??
			ID2 = record2.annotations["accessions"][0] #.replace("_","")
#			print(ID2)
			try:
				newname = record2.annotations["taxonomy"][1] + "\t" +record2.annotations["taxonomy"][2] + "\t" +record2.annotations["taxonomy"][3] + "\t" +record2.annotations["taxonomy"][4] +  "\t" +  record2.annotations["organism"]  # need to be modified related to the taxonomy list =>record2.annotations["taxonomy"][14][0:3]
			except:
				try: 
					newname = record2.annotations["taxonomy"][1] + "\t" +record2.annotations["taxonomy"][2] + "\t" +record2.annotations["taxonomy"][3] + "\t\t" +  record2.annotations["organism"] # need to be modified related to the taxonomy list =>record2.annotations["taxonomy"][14][0:3]
				except:
					try: 
						newname = record2.annotations["taxonomy"][1]+ "\t" +record2.annotations["taxonomy"][2]+ "\t\t\t" +  record2.annotations["organism"] # need to be modified related to the taxonomy list =>record2.annotations["taxonomy"][14][0:3]
					except:
						try: 
							newname = record2.annotations["taxonomy"][1]+ "\t\t\t\t" +  record2.annotations["organism"] # need to be modified related to the taxonomy list =>record2.annotations["taxonomy"][14][0:3]
						except:
							newname = "\t\t\t\t\t" + record2.annotations["organism"].replace(' ',' ')
#			print(newname)
			rename2 = newname.translate(''.join(chars_to_remove)) 
			OTUlist= []
			for key, value in blastrecord.items():
				if value[0] == ID2:
					OTUlist.append(key)
					OTUIDlist2.append(key)
#					print(key, value) 
					outrename = open("BLAST_results/Blast_result_renamed.txt",'a')
					outrename.write(key + '\t' + rename2 + '\t'+ ID2 +  '\t'+ value[1].replace('_','\t')+'\n')
					outrename.close()
			print(len(OTUlist), " OTUs match with => " , rename2)
			log = open("log.txt",'a')
			log.write(str(len(OTUlist))+ "\t OTUs match with  \t" + rename2 + '\n')
			log.close()
			
	for OTU in OTUIDlist:
		if OTU not in OTUIDlist2:
			print(OTU, blastrecord[OTU]) 
			
def main():
	script,  NGSfile, idminy, Eminz = argv
	idmin = float(idminy)
	Emin = float(Eminz)
	blastonline(NGSfile,idmin, Emin)
	Rename(NGSfile)
main()