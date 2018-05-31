import os, re, string
import xml.etree.ElementTree

def main():
	search = input('Have you changed the name, password and key in the script (look at file lines 10,11 and 12, replace by your codes)? Y or N \n')
	if search[0] == 'N':
		main()
	if search[0] == "Y":
		URL="https://cipresrest.sdsc.edu/cipresrest/v1"
		CRA_USER="your user name"
		PASSWORD="your password"
		search = input('Do you want to build a new tree? Y or N \n')
	if search[0] == 'Y':
		search = input('Have you changed the name of the alignment in the script (look at file line 23, replace test_tree_pipeline.fasta by your sequence file)? Y or N \n')
		if search[0] == 'N':
			main()
		if search[0] == "Y":
#start new tree:
			file = 'test_tree_pipeline.fasta'
# for protein add '-F vparam.datatype_=protein' and replace '-F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I' by '-F vparam.prot_sub_model_=PROTGAMMA -F vparam.prot_matrix_spec_=LG'
# see http://www.phylo.org/index.php/rest/raxmlhpc8_rest_xsede.html for more details
# with 10 bootstraps
#			print('curl -u $CRA_USER:$PASSWORD -H cipres-appkey:$KEY $URL/job/$CRA_USER -F tool=RAXMLHPC8_REST_XSEDE -F vparam.select_analysis_=fa -F vparam.choose_bootstrap_=x -F vparam.choose_bootstop_=specify -F vparam.bootstrap_value_=10 -F input.infile_=@./' +file+' -F vparam.runtime_=168 -F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I  -F vparam.provide_parsimony_seed_=1 -F vparam.parsimony_seed_val_=12345 -F vparam.outsuffix_=' + file.split('.')[0] + '_outrax.tree -F metadata.statusEmail=true > temp.xml')
#			os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL+"/job/"+CRA_USER+" -F tool=RAXMLHPC8_REST_XSEDE -F vparam.select_analysis_=fa -F vparam.choose_bootstrap_=x -F vparam.choose_bootstop_=specify -F vparam.bootstrap_value_=10 -F input.infile_=@./" +file+" -F vparam.runtime_=168 -F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I -F vparam.provide_parsimony_seed_=1 -F vparam.parsimony_seed_val_=12345 -F vparam.outsuffix_=" + file.split('.')[0] + "_outrax.tree -F metadata.statusEmail=true > temp.xml") 
# without bootstrap
			print('curl -u $CRA_USER:$PASSWORD -H cipres-appkey:$KEY $URL/job/$CRA_USER -F tool=RAXMLHPC8_REST_XSEDE -F vparam.select_analysis_=fd -F input.infile_=@./' +file+' -F vparam.runtime_=48 -F vparam.runtime_=168 -F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I -F vparam.parsimony_seed_val_=12345 -F vparam.outsuffix_=' + file.split('.')[0] + '_outrax.tree -F metadata.statusEmail=true > temp.xml')
			os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL+"/job/"+CRA_USER+" -F tool=RAXMLHPC8_REST_XSEDE -F vparam.select_analysis_=fd -F input.infile_=@./" +file+" -F vparam.runtime_=168 -F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I -F vparam.parsimony_seed_val_=12345 -F vparam.outsuffix_=" + file.split('.')[0] + "_outrax.tree -F metadata.statusEmail=true > temp.xml") 
# with OTU and reference tree
#			refTree = 'referenceTree.tre'
#			print('curl -u $CRA_USER:$PASSWORD -H cipres-appkey:$KEY $URL/job/$CRA_USER -F tool=RAXMLHPC8_REST_XSEDE -F vparam.select_analysis_=fv -F input.infile_=@./' +file+' -F vparam.runtime_=48 -F vparam.runtime_=168 -F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I -F input.treetop_=@./'+refTree+' -F vparam.outsuffix_=' + file.split('.')[0] + '_outrax.tree -F metadata.statusEmail=true > temp.xml')
#			os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL+"/job/"+CRA_USER+" -F tool=RAXMLHPC8_REST_XSEDE -F vparam.select_analysis_=fv -F input.infile_=@./" +file+" -F vparam.runtime_=168 -F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I -F input.treetop_=@./"+refTree+" -F vparam.outsuffix_=" + file.split('.')[0] + "_outrax.tree -F metadata.statusEmail=true > temp.xml") 

# temp.xml store the commandline of the last job you sent to CIPRES REST API
	
			subxml = xml.etree.ElementTree.parse('temp.xml').getroot()
			for field in subxml.iter('title'):
				if 'RAXML' in field.text:
					print(field.text)
					out= open('List_of_Submitted_Tree.txt','a')
					out.write(file+ '\t'+ str(field.text) +'\n')
					out.close()
	if search[0] == 'N':
		search2 = input('downloading a tree? Y or N \n')
		if search2[0] == 'Y':
#get back the tree:
			listfiles = open('List_of_Submitted_Tree.txt','r')
			for line in listfiles:
				os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL+"/job/"+CRA_USER+"/"+line.split('\n')[0].split('\t')[1]+'/output > outputlist.xml')
		
			outxml = xml.etree.ElementTree.parse('output.xml').getroot()	
			for results in outxml.iter('jobfile'):
				for resultname in results.iter('filename'):
					print(resultname.text)
#					if '.' in resultname.text: # for downloading everything
					if 'BestTree' in resultname.text:
						for resultcode in results.iter('outputDocumentId'):
							print(resultcode.text)
							os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL+"/job/"+CRA_USER+"/"+line.split('\n')[0].split('\t')[1]+'/output/'+str(resultcode.text)+' > '+str(resultname.text)) 
		if search2[0] == 'N':
			print ('Please answer yes or no. Try again. ')
main()
