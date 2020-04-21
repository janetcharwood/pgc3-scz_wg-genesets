# pgc3-scz_wg-genesets

1. GO gene sets

Details of the construction of the GO gene sets are in ANNOTATION_AUTO/README/AUTOMATED_ANNOTATION_README.rtf section 6. GO: constructing gene sets from GO annotations.

GO gene sets were derived using the annotation files downloaded on 04_07_2018 and using the scripts in the automated annotation process:ANNOTATION_AUTO

The file GO_STRICT_PC_10-2000_MAGMA.txt which contains sets of between 10 and 2000 protein coding genes derived from GO 
annotations was used in further analysis. Entrez gene ids were converted to ensembl gene ids using the tools in the folder: Gene_sets_to_ensembl_220618: The final gene set file used in the analysis was GO_STRICT_PC_10-2000_MAGMA_ensembl_magma.txt


2. Developmental gene sets.

Shank3_Dlg4_Dlgap1_dev_gene_sets.txt contains the human gene sets derived from each of the baits, Shank3, Dlg4 and Dlgap1 and each of the developmental stages: e14,p7,p14 and adult.

Shank3_Dlg4_Dlgap1_dev_gene_UpSet_subsets.txt (Supplementary Figure 11_ PSD Upset_v5) contains the subsets of non-overlapping gene sets derived from the UpSet plots.  Set notation is used in the gene set names. 
