
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Publicly available gene sets used in PGC3-SCZ

This repository contains code that can be used to automatically extract
biological annotations from the [GO](http://geneontology.org/) and
[MGI](http://www.informatics.jax.org) databases. It also contains a
cleaned and Ensembl-translated version of GO suitable to be used with
[MAGMA](https://ctg.cncr.nl/software/magma).

## Gene Ontology gene sets

Details of the construction of the GO gene sets are given in section 6
of
[AUTOMATED\_ANNOTATION\_README.rtf](ANNOTATION_AUTO/README/AUTOMATED_ANNOTATION_README.rtf):
*“GO: constructing gene sets from GO annotations”*.

The file
[GO\_STRICT\_PC\_10-2000\_MAGMA\_ensembl.txt](GO_STRICT_PC_10-2000_MAGMA_ensembl.txt)
contains sets of between 10 and 2000 protein coding genes derived from
GO annotations (release 09/11/2020), as described in the
[PGC3-SCZ](https://doi.org/10.1101/2020.09.12.20192922) paper. Entrez
gene ids were converted to Ensembl gene ids using the tools in the
folder [Gene\_sets\_to\_ensembl\_220618](Gene_sets_to_ensembl_220618/).
Obsolete gene sets that remained in the raw and processed GO files have
been removed from this file, and the entries are given in the file
[OBSOLETE\_GO\_TERMS.txt](OBSOLETE_GO_TERMS.txt)

## Citation

If any of these scripts or data files is helpful for your work, please
reference the main
[PGC3-SCZ](https://doi.org/10.1101/2020.09.12.20192922) paper.

## Contact

Please submit suggestions and bug-reports at
<https://github.com/janetcharwood/pgc3-scz_wg-genesets/issues>.
