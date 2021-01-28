Written by Janet Harwood June 22nd 2018.
----------------------------------------

Script gene_sets_entrez_to_ensembl_long.py takes long format gene sets with entrez gene ids and converts them to ENSEMBL gene ids.

Script gene_sets_entrez_to_ensembl_wide.py takes gene sets with one gene set per line with entrez gene ids and converts them to ENSEMBL gene ids.


This script will run from any directory.

The input files must be in the same directory as the script

ENSG.genes.txt
ncbi_human_entrez_to_ensembl_all_18_6_2018.txt

the gene sets must be in a folder called ‘gene_sets’

the output goes to a folder called ‘gene_sets_ensembl’


To run this script on the command line:

python gene_sets_entrez_to_ensembl_long.py

This script is dependent on the python packages: ‘re’ and ‘ ‘default dict’


The ENSG.genes.txt file was used to do the id conversion from NCBI genes to ENSEMBL. It was derived by Dr. Kyoko Watanabe, from Prof. Danielle Posthuma's group (https://ctg.cncr.nl/research_teams/psychiatric_and_statistical_genomics/).

This is her code:

library(biomaRt)
version = 92
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37, version=version)
chr <- c(1:22, "X")
attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name","start_position", "end_position", "strand", "gene_biotype", "hgnc_symbol", "mim_gene_accession")
ENSG.genes <- getBM(attributes = attributes, filters="chromosome_name", values=chr, mart=ensembl)
ENSG.genes <- ENSG.genes[with(ENSG.genes, order(chromosome_name, start_position)),]
ENSG.genes <- ENSG.genes[order(ENSG.genes$chromosome_name),]

The NCBI mapping file ncbi_human_entrez_to_ensembl_all_18_6_2018.txt was derived from the NCBI gene info file by Dr. Janet Harwood at Cardiff University.

It contains all human genes: entrez gene id: symbol: ENSEMBL gene id: gene type


The converted gene sets are output with the file name: *.ensembl_magma.txt

Unmapped ids are output in the file *.ensembl_unmapped.txt


the  *.ensembl_unmapped.txt contains mapping information from both annotation files:

format of *.ensembl_unmapped.txt file:

entrez gene id: ensembl gene id: type of gene: gene symbol: annotation source

ensembl_GRCh37 = ENSG.genes.txt file

ncbi = ncbi_human_entrez_to_ensembl_all_18_6_2018.txt
