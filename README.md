
# Phytocytokine Scripts

# Aim of the project

The aim of this project was the identification of orthologs of phytocytokines in *Zea mays*.

# Mining of cross-species phytocytokines

Therefore a BLAST database containing 79 predicted plant proteoms retrieved from the ensemble plants database (https://plants.ensembl.org/info/data/ftp/index.html) was created using the NCBI BLAST+ 2.2.18 application was created. 

Additionally, 18 previously published phytocytokines amino acid sequences were selected as queries for the BLAST searches.

psiBLAST searches of the queries in the database were performed with five iterations and an e-value threshold of 0.05 to detect distant relationships between proteins. A parallel BLASTp with the same settings was performed.

The following steps were performed with the script "psiblast_ensemble_propeptides_analysis.py" if not noted differentely.

After the BLASTs, the results for all queried propeptides were combined, and unique identifiers were assigned to each hit linking it to the respective query. 

The psiBLAST of all queried propeptides included many redundant hits (a hit found in rounds 1 to 5 was in the results five times) and could not be compared to the BLASTp results. Therefore, these redundant hits were removed after sorting all hits by their bit score.
Additionally hits with low quality were removed (below 10% protein identities, query coverage below 25%) by the "process_data" function.

Therefore the script "ensemble_accession_organism.py" was used to link each id with the source proteom. 
The created table was loaded and linked to the hits with the "create_map" function in the script.

The remaining hits varied strongly in their amino acid sequence length. Therefore, all hits whose amino acid sequence length was 1.5-times longer or 1.5-times shorter than the longest or shortest hit within the group of hits with more than 75% identity, respectively, were removed (function "rel_length", length filtering).

In the last step, hits were filtered according to their similarity to the respective peptide hormone of each query.
Initial multiple sequence alignments (MSA) with ClustalOmega 1.2.2 (prepared with function "create_clustal_cmd") of the hits indicated that there were still some false-positive hits without an alignment in the region of the peptide hormone, which might disturb the MSA and result in large gaps within the MSA.
Therefore, all hits for each queried propeptide were aligned with the ClustalOmega application, and the region containing the peptide hormone was analyzed and scored (functions "read_clustal", "clustal_analysis", and "filter_data").
The similarity of each alignment region containing the peptide hormone was scored based on a BLOSUM62 matrix and low gap costs (gap: -8, extension: -1) using the biopython package:

Later motif diagrams were created with the logomaker 0.8 package (function "create_fasta" and "create_weblogo").

# Figures 

# Phylogenetic tree

The phylogenetic tree of the subset of organisms was created as follows: 
140 complete proteins present in all proteoms of the selected organisms either as one or multiple copies were selected using BUSCO 5.1.2 in its standard settings. 
If multiple copies were present, only the best hit was used to create the phylogenetic tree. 
The orthologs of each BUSCO protein were aligned via ClustalOmega 1.2.2. 
The alignments of all proteins were concatenated. 
A phylogenetic tree of the concatenated multiple sequence alignments was calculated via RAxML 8.2.12 using the standard settings. 
These procedures were performed with the python script ”psiblast_ensemble_propeptides_analysis.py”.

# Venn diagram

The venn diagram was created with "venn_diagram.R".

# qPCR graphs

Graphs displaying the fold change of gene expression were created with "3h+24h_foldchange.R".

# Motif comparison AtPEP1 and ZmPEP1

Protein sequences of hits for AtPEP1 and ZmPEP1 were clustered together with CLUSTAL Omega. 
Afterwards, results were seperated based on the query. 
Per query the peptide sequence was extracted and a weblogo motif was created with the script "weblogo.py".