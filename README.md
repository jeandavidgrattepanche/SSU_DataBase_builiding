# SSU_DataBase_builiding
Here are a set of python scripts I generated to built a curated SSU database, including removing redundant sequences and removing column with x% of missing characters

The sequences for the taget lineage(s) are download from GenBank (or other public database).


The first set of scripts located in byTree folder use a tree approach to remove too similar sequences.
Required: MAFFT, FastTree, Water (EMBOSS) 

The second script located in byClustering folder uses a clustering approache to remove too simialr sequences.
Required: Vsearch

The scripts require the software packages to be added to your path.

mask_gaps.py allows to remove column with missing character at x %. 
Required: BioPython and aligned sequences file.
