# SSU_DataBase_builiding
Here are a set of python scripts I generated to built a curated SSU database, including removing redundant sequences and removing column with x% of missing characters

The sequences for the taget lineage(s) are downloaded from GenBank (or other public database).


The first set of scripts located in byTree folder use a tree approach to remove too similar sequences.
Required: MAFFT, FastTree, Water (EMBOSS) 

The second script located in byClustering folder uses a clustering approache to remove too simialr sequences.
Required: Vsearch

The scripts require the software packages to be added to your path.

constraintTree.py allows to create a constraint tree to use in RaxML. In this purpose, the sequence names are coded such as:

Al_ci_in_sp_Strombidium_rassoulzadegani_AY257125 : Alveolata; Ciliophora; Intramacronucleata; Spirotrichea; and the species name and GenBank accession number
Am_My_My_Didymium_nigripes_AF239230: Amoebozoa; Mycetozoa; Myxogastria; and the species name and GenBank accession number
The script uses this code to constraint for outgroup.
This script generates a file named constrainttree.txt

mask_gaps.py allows to remove column with missing character at x % from an alignment. 
Required: BioPython and aligned sequences file. I use MAFFT to do my alignments ($ MAFFT inputsequencefile > outputalignedsequencefile.fas)
to run the script: $ python mask_gaps.py outputalignedsequencefile.fas x%
generate a file named: outputalignedsequencefile_maskedx.fas


These database can be use within the MiSeq pipeline (see other repository). For this, you need to create a folder with:
- the sequence file before alignment (inputsequencefile)
- the alignment outputalignedsequencefile.fas
- the masked alignment outputalignedsequencefile_maskedx.fas
- a tree for reference in the pipeline, I use RAxML to do my tree ($ raxmlHPC-PTHREADS-AVX2 -s outputalignedsequencefile_maskedx.fas -m GTRGAMMAI -g constrainttree.txt -n ContraintTree.tre -p 3498756 (a random number)-T -1 (number of threads you have -1)
- BLAST database ($ makeblastdb -in inputsequencefile -dbtype nucl)
