import os, re, string
import xml.etree.ElementTree

def main():
	search0 = input('Have you changed the name, password and key in the script (look at file lines 10,11 and 12, replace by your codes)? Y or N \n')
	if search0[0] == 'N':
		main()
	if search0[0] == "Y":
		URL="https://cipresrest.sdsc.edu/cipresrest/v1"
		CRA_USER="CRA_USER"
		PASSWORD="PASSWORD"
		KEY="KEY"
  search = input('Do you want to send a new job? Y or N \n')
		if search[0] == 'Y':
			search1 = input('Have you changed the name of the file/alignment in the script (look at file line 20, replace test_tree_pipeline.fasta by your sequence file), and the reference l 30 or 34 if needed ? Y or N \n')
			if search1[0] == 'N':
				main()
			if search1[0] == "Y":
			#start new tree:
				file = 'Am_ref_28June18_clean1_mafft_masked_90.fas'
# for protein add '-F vparam.datatype_=protein' and replace '-F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I' by '-F vparam.prot_sub_model_=PROTGAMMA -F vparam.prot_matrix_spec_=LG'
# see http://www.phylo.org/index.php/rest/raxmlhpc8_rest_xsede.html for more details and https://www.phylo.org/restusers/docs/cipresXml
# MAFFT auto for quick alignment
# 				print("curl -u CRA_USER:PASSWORD -H cipres-appkey:KEY URL/job/CRA_USER -F metadata.statusEmail=true -F tool=MAFFT_XSEDE -F input.infile_=@./" +file+" -F vparam.auto_analysis_='1'  -F vparam.adjust_direction_='1' -F vparam.anysymbol_='1' -F vparam.datatype_='dna' -F vparam.dnaMatrix_='200' -F vparam.noScore_='1' -F vparam.opPenaltyGroupToGroup_='1.53' -F vparam.outputOrder_='--reorder' -F vparam.which_mafft_='7394' -F vparam.runtime_='168' -F vparam.more_memory_='1'")
# 				os.sytem("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL+"/job/"+CRA_USER+" -F metadata.statusEmail=true -F tool=MAFFT_XSEDE -F input.infile_=@./" +file+" -F vparam.auto_analysis_='1' -F vparam.adjust_direction_='1' -F vparam.anysymbol_='1' -F vparam.datatype_='dna' -F vparam.dnaMatrix_='200' -F vparam.noScore_='1' -F vparam.opPenaltyGroupToGroup_='1.53' -F vparam.outputOrder_='--reorder' -F vparam.which_mafft_='7394' -F vparam.runtime_='168' -F vparam.more_memory_='1'")
# 
# MAFFT Linsi for accurate alignment
# 				print("curl -u CRA_USER:PASSWORD -H cipres-appkey:KEY URL/job/CRA_USER -F metadata.statusEmail=true -F tool=MAFFT_XSEDE -F input.infile_=@./" +file+" -F vparam.accurate_executable_='linsi68' -F vparam.adjust_direction_='1' -F vparam.analysis_type_='accurate' -F vparam.anysymbol_='1' -F vparam.datatype_='dna' -F vparam.dnaMatrix_='200' -F vparam.noScore_='1' -F vparam.opPenaltyGroupToGroup_='1.53' -F vparam.outputOrder_='--reorder' -F vparam.which_mafft_='7394' -F vparam.runtime_='168' -F vparam.more_memory_='1'")
# 				os.sytem("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL+"/job/"+CRA_USER+" -F metadata.statusEmail=true -F tool=MAFFT_XSEDE -F input.infile_=@./" +file+" -F vparam.accurate_executable_='linsi68' -F vparam.adjust_direction_='1' -F vparam.analysis_type_='accurate' -F vparam.anysymbol_='1' -F vparam.datatype_='dna' -F vparam.dnaMatrix_='200' -F vparam.noScore_='1' -F vparam.opPenaltyGroupToGroup_='1.53' -F vparam.outputOrder_='--reorder' -F vparam.which_mafft_='7394' -F vparam.runtime_='168' -F vparam.more_memory_='1'")
# MAFTT --addfragment
#  				refalignment = 'Am_ref_28June18_clean1_mafft_masked_90.fas'
# 				print("curl -u CRA_USER:PASSWORD -H cipres-appkey:KEY URL/job/CRA_USER -F metadata.statusEmail=true -F tool=MAFFT_XSEDE -F input.infile_=@./" +file+" -F input.ref_alignment1_=@./" +refalignment+" -F vparam.anysymbol_='1' -F vparam.auto_analysis_='1' -F vparam.noScore_='1' -F vparam.opPenaltyGroupToGroup_='1.53' -F vparam.outputOrder_='--reorder' -F vparam.reorder_add_='1' -F vparam.use_addfrag_='1' -F vparam.use_keeplength_='1' -F vparam.use_mapout_='1' -F vparam.which_mafft_='7394' -F vparam.runtime_='168' -F vparam.more_memory_='1'")
# 				os.sytem("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL+"/job/"+CRA_USER+" -F metadata.statusEmail=true -F tool=MAFFT_XSEDE -F input.infile_=@./" +file+" -F input.ref_alignment1_=@./" +refalignment+" -F vparam.anysymbol_='1' -F vparam.auto_analysis_='1' -F vparam.noScore_='1' -F vparam.opPenaltyGroupToGroup_='1.53' -F vparam.outputOrder_='--reorder' -F vparam.reorder_add_='1' -F vparam.use_addfrag_='1' -F vparam.use_keeplength_='1' -F vparam.use_mapout_='1' -F vparam.which_mafft_='7394' -F vparam.runtime_='168' -F vparam.more_memory_='1'")
# with 10 bootstraps
#				print('curl -u $CRA_USER:$PASSWORD -H cipres-appkey:$KEY $URL/job/$CRA_USER -F tool=RAXMLHPC8_REST_XSEDE -F vparam.select_analysis_=fa -F vparam.choose_bootstrap_=x -F vparam.choose_bootstop_=specify -F vparam.bootstrap_value_=10 -F input.infile_=@./' +file+' -F vparam.runtime_=168 -F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I  -F vparam.provide_parsimony_seed_=1 -F vparam.parsimony_seed_val_=12345 -F vparam.outsuffix_=' + file.split('.')[0] + '_outrax.tree -F metadata.statusEmail=true > temp.xml')
#				os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL+"/job/"+CRA_USER+" -F tool=RAXMLHPC8_REST_XSEDE -F vparam.select_analysis_=fa -F vparam.choose_bootstrap_=x -F vparam.choose_bootstop_=specify -F vparam.bootstrap_value_=10 -F input.infile_=@./" +file+" -F vparam.runtime_=168 -F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I -F vparam.provide_parsimony_seed_=1 -F vparam.parsimony_seed_val_=12345 -F vparam.outsuffix_=" + file.split('.')[0] + "_outrax.tree -F metadata.statusEmail=true > temp.xml") 
# without bootstrap
# 				print('curl -u $CRA_USER:$PASSWORD -H cipres-appkey:$KEY $URL/job/$CRA_USER -F tool=RAXMLHPC8_REST_XSEDE -F vparam.select_analysis_=fd -F input.infile_=@./' +file+' -F vparam.runtime_=48 -F vparam.runtime_=168 -F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I -F vparam.parsimony_seed_val_=12345 -F vparam.outsuffix_=' + file.split('.')[0] + '_outrax.tree -F metadata.statusEmail=true > temp.xml')
# 				os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL+"/job/"+CRA_USER+" -F tool=RAXMLHPC8_REST_XSEDE -F vparam.select_analysis_=fd -F input.infile_=@./" +file+" -F vparam.runtime_=168 -F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I -F vparam.parsimony_seed_val_=12345 -F vparam.outsuffix_=" + file.split('.')[0] + "_outrax.tree -F metadata.statusEmail=true > temp.xml") 
# constraint tree
# 				constraintTree = 'constrainttree.tre'
# 				print('curl -u $CRA_USER:$PASSWORD -H cipres-appkey:$KEY $URL/job/$CRA_USER -F tool=RAXMLHPC8_REST_XSEDE -F vparam.select_analysis_=fd -F input.infile_=@./' +file+' -F vparam.runtime_=48 -F vparam.runtime_=168 -F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I -F input.constraint_=@./'+constraintTree+' -F vparam.outsuffix_=' + file.split('.')[0] + '_constT.tree -F metadata.statusEmail=true > temp.xml')
# 				os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL+"/job/"+CRA_USER+" -F tool=RAXMLHPC8_REST_XSEDE -F vparam.select_analysis_=fd -F input.infile_=@./" +file+" -F vparam.runtime_=168 -F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I -F input.constraint_=@./"+constraintTree+" -F vparam.outsuffix_=" + file.split('.')[0] + "_constT.tree -F metadata.statusEmail=true > temp.xml") 
# with OTU and reference tree
# 				refTree = 'RAxML_bestTree.SSU_SAR_only_v14_for_Euk_SAR_Rhi_Fo_Outgroups_mafft_outrax_85.tree'
# 				print('curl -u $CRA_USER:$PASSWORD -H cipres-appkey:$KEY $URL/job/$CRA_USER -F tool=RAXMLHPC8_REST_XSEDE -F vparam.select_analysis_=fv -F input.infile_=@./' +file+' -F vparam.runtime_=48 -F vparam.runtime_=168 -F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I -F input.treetop_=@./'+refTree+' -F vparam.outsuffix_=' + file.split('.')[0] + '_outrax.tree -F metadata.statusEmail=true > temp.xml')
# 				os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL+"/job/"+CRA_USER+" -F tool=RAXMLHPC8_REST_XSEDE -F vparam.select_analysis_=fv -F input.infile_=@./" +file+" -F vparam.runtime_=168 -F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I -F input.treetop_=@./"+refTree+" -F vparam.outsuffix_=" + file.split('.')[0] + "_outrax.tree -F metadata.statusEmail=true > temp.xml") 

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
		
					outxml = xml.etree.ElementTree.parse('outputlist.xml').getroot()	
					for results in outxml.iter('jobfile'):
						for resultname in results.iter('filename'):
							print(resultname.text)
# 							if '.txt' in resultname.text: # for downloading everything
							if 'mafft' in resultname.text: # for downloading everything
# 							if 'bestTree' in resultname.text:
# 							if 'labelledTree' in resultname.text:
								for resultcode in results.iter('outputDocumentId'):
									print(resultcode.text)
									os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" "+URL+"/job/"+CRA_USER+"/"+line.split('\n')[0].split('\t')[1]+'/output/'+str(resultcode.text)+' > '+str(line.split('\n')[0].split('\t')[0]+'_'+resultname.text)) 
				print("When you have done, please remember to delete your file on the CIPRES server. To do it copy the job you want to delete in a file named List_of_Tree_to_be_deleted.txt. Thanks!")
			if search2[0] == 'N':
				search3 = input('Deleting old Tree? Y or N \n')
				if search3[0] == 'Y':
					listfiles = open('List_of_Tree_to_be_deleted.txt','r')
					for line in listfiles:
						os.system("curl -u "+CRA_USER+":"+PASSWORD+" -H cipres-appkey:"+KEY+" -X DELETE "+URL+"/job/"+CRA_USER+"/"+line.split('\n')[0].split('\t')[1]+' > deleted.xml')
				if search3[0] == 'N':
					print ('Please answer yes or no. Try again. ')
					main()
main()
