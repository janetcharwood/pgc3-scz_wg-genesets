#!/bin/bash
#$ -V
#$ -q all.q
#$ -l h_vmem=40G
#$ -cwd
#$ -o annotation_log
#$ -e annotation_err
#$ -S /bin/bash

export PATH=/share/apps/anaconda2/bin:$PATH

# activate numpy pandas environment

source activate py2-numpy-pandas

# # create directories.
root_dir=$(pwd)

mkdir processed
mkdir logs
mkdir downloads

cd processed

# create sub directories.

sub_dir_list=(GO_ANNOTATIONS GO_ATTRIBUTES GO_ATTRIBUTES_UPDATED \
GO_EVIDENCE GO_PATHS GO_TREES HOMOLOGENE GENE_SETS GO_GENE_CHECK MP_ANNOTATIONS \
MP_EVIDENCE MP_ID_MAPPING MP_TREES)

for subdir in ${sub_dir_list[*]}
do
    mkdir $subdir
done

cd GENE_SETS

mkdir GO
mkdir MAMMALIAN_PHENOTYPE

cd GO
mkdir GO_MAGMA_FORMAT

cd ..

cd MAMMALIAN_PHENOTYPE
mkdir MOUSE

cd MOUSE
mkdir MOUSE_MAGMA_FORMAT

###################################################################

cd $root_dir/scripts

python ./annotation_downloads_hbp.py
################################################################

# check that all of the downloaded files exist before proceeding.

cd $root_dir/downloads

files=(MGI_EntrezGene.rpt  gene_info.gz MGI_PhenoGenoMP.rpt gene2go.gz homologene.data MPheno_OBO.ontology go-basic.obo)

for file in ${files[*]}
do
  if [ -s $file ]
  then
    echo "$0: File '${file}' has downloaded."
  else
    echo "$0: File '${file}' does not exist or is empty."
    exit 1
  fi
done

################################################################################
################################################################################

## process mouse annotations
# cd to the directory containing all the scripts to process the data.
cd $root_dir/scripts

# MP ID mapping

python ./MGI_Marker_ID_to_entrez_hbp.py

cd $root_dir/processed/MP_ID_MAPPING

# check that the output files exist
ID_mapping_files=(MGI_markerID_to_entrezID_ALL.txt	MGI_markerID_to_entrezID_pc.txt)


for file in ${ID_mapping_files[*]}
do
  if [ -s $file ]
  then
    echo "$0: File '${file}' exists and is not empty."
  else
    echo "$0: File '${file}' does not exist or is empty."
    exit 1
  fi
done

###############################################################################
## cd to the directory containing all the scripts to process the data.
cd $root_dir/scripts

#construct the evidence files:

python ./Mouse_pheno_JH_hbp.py

# check that the evidence files exist
cd $root_dir/processed/MP_EVIDENCE

evidence_files=(MGI_PhenoGeno_single_gene_ALL.txt	MGI_PhenoGeno_single_protein_coding_gene.txt)


for file in ${evidence_files[*]}
do
  if [ -s $file ]
  then
    echo "$0: File '${file}' exists and is not empty."
  else
    echo "$0: File '${file}' does not exist or is empty."
    exit 1
  fi
done

################################################################################
## Extracting mouse phenotype ontology to tree

cd $root_dir/scripts

python ./obo_to_tree_mouse_hbp.py

# check that the tree files exist
cd $root_dir/processed/MP_TREES

mouse_trees_files=(MP_tree.txt MP_attr.txt)

for file in ${mouse_trees_files[*]}
do
  if [ -s $file ]
  then
    echo "$0: File '${file}' exists and is not empty."
  else
    echo "$0: File '${file}' does not exist or is empty."
    exit 1
  fi
done

###############################################################################
# Process mouse trees files to create paths and update the attributes files.
cd $root_dir/scripts

python ./obo_tree_to_paths_hbp.py

# check that the path and updated attributes files exist
cd $root_dir/processed/MP_TREES

mouse_paths_files=(MP_paths.txt MP_attr_level.txt)

for file in ${mouse_paths_files[*]}
do
  if [ -s $file ]
  then
    echo "$0: File '${file}' exists and is not empty."
  else
    echo "$0: File '${file}' does not exist or is empty."
    exit 1
  fi
done

##############################################################################
## expand the mouse annotations.
cd $root_dir/scripts

python ./expand_pheno_mouse_hbp.py

# check that the final output files 'expand.txt' exist
cd $root_dir/processed/MP_ANNOTATIONS

mouse_expand_files=(MGI_single_gene_Pheno_protein_coding_annotation.txt)

for file in ${mouse_expand_files[*]}
do
  if [ -s $mouse_expand_files ]
  then
    echo "$0: File '${mouse_expand_files}' exists and is not empty."
  else
    echo "$0: File '${mouse_expand_files}' does not exist or is empty."
    exit 1
  fi
done

##############################################################################
echo "MP_PHENO Processing complete at $(date)"

########################################################################
########################################################################
# Process the homologene data
cd $root_dir/downloads

# extract mouse and human data as two separate files:

awk '$2 == "9606" {print $0  }' homologene.data > $root_dir/processed/HOMOLOGENE/human_homologene_jh.txt

awk '$2 == "10090" {print $0  }' homologene.data > $root_dir/processed/HOMOLOGENE/mouse_homologene_jh.txt

# check that the homologene_cut files exist
cd $root_dir/processed/HOMOLOGENE

cut_files=(human_homologene_jh.txt mouse_homologene_jh.txt )

for file in ${cut_files[*]}
do
  if [ -r $file ]
  then
    echo "$0: File '${file}' exists."
  else
    echo "$0: File '${file}' does not exist."
    exit 1
  fi
done

######################

# extract human and mouse gene information from the NCBI gene file and make homologene gene sets

cd $root_dir/scripts

python ./homologene_gene_set_extractor_hbp.py

# check that the human and mouse NCBI gene files exist

cd $root_dir/processed/HOMOLOGENE

homologene_gene_sets=(homologene_human_protein_coding.txt homologene_mouse_protein_coding.txt homologene_human_all_minus_pseudo.txt homologene_mouse_all_minus_pseudo.txt)

for file in ${homologene_gene_sets[*]}
do
  if [ -s $file ]
  then
    echo "$0: File '${file}' exists and is not empty."
  else
    echo "$0: File '${file}' does not exist or is empty."
    exit 1
  fi
done

###########################################################################
## extract the one-to-one, one-to many etc relationships between the human and mouse homologene gene sets.

cd $root_dir/scripts

python ./homologene_merge_hbp.py

# check that the homologene_merge files exist
cd $root_dir/processed/HOMOLOGENE

homologene_merge_files=(hm_one_to_one_homol_ALL.txt hm_one_to_many_homol_ALL.txt hm_many_to_many_homol_ALL.txt  hm_many_to_one_homol_ALL.txt hm_one_to_one_homol_PC.txt hm_one_to_many_homol_PC.txt  hm_many_to_one_homol_PC.txt hm_many_to_one_homol_PC.txt hm_many_to_many_homol_PC.txt)

for file in ${homologene_merge_files[*]}
do
  if [ -s $file ]
  then
    echo "$0: File '${file}' exists and is not empty."
  else
    echo "$0: File '${file}' does not exist or is empty."
    exit 1
  fi
done

############################################################################
# Extract human homologues of mouse genes from phenotypes with single gene manipulations.

cd $root_dir/scripts

python ./Mouse_pheno_to_human_PC_gene_hbp.py

# check that the mouse_to_human files exist

cd $root_dir/processed/GENE_SETS/MAMMALIAN_PHENOTYPE/MOUSE

mouse_pheno_to_human_gene_files=(MGI_single_gene_Pheno_to_human_protein_coding_gene.txt)

for file in ${mouse_pheno_to_human_gene_files[*]}
do
  if [ -s $file ]
  then
    echo "$0: File '${file}' exists and is not empty."
  else
    echo "$0: File '${file}' does not exist or is empty or is empty."
    exit 1
  fi
done

echo "homologene processing complete at $(date)"

###############################################################################
###############################################################################
## Process GO ontology.
# cd to the directory containing all the scripts to process the data.
cd $root_dir/scripts

# construct the evidence files:

python ./human_gene2go_extract_hbp.py

# check that the evidence files exist
cd $root_dir/processed/GO_EVIDENCE

evidence_files=(Gene_to_GO_ALL_ev_ALL_genes_evidence.txt Gene_to_GO_ALL_ev_PC_genes_evidence.txt Gene_to_GO_STRICT_ev_ALL_genes_evidence.txt Gene_to_GO_STRICT_ev_PC_genes_evidence.txt )

for file in ${evidence_files[*]}
do
  if [ -s $file ]
  then
    echo "$0: File '${file}' exists and is not empty."
  else
    echo "$0: File '${file}' does not exist or is empty or is empty."
    exit 1
  fi
done

###############################################################################
# Process the GO obo file to generate three ontology trees and their attribute files
cd $root_dir/scripts

python ./GO_obo_to_tree_hbp.py

# check that the trees files exist
cd $root_dir/processed/GO_TREES

trees_files=(BP_tree.txt CC_tree.txt MF_tree.txt)

for file in ${trees_files[*]}
do
  if [ -s $file ]
  then
    echo "$0: File '${file}' exists and is not empty."
  else
    echo "$0: File '${file}' does not exist or is empty or is empty."
    exit 1
  fi
done

#check that the attributes files exist
cd $root_dir/processed/GO_ATTRIBUTES

attr_files=(BP_attr.txt	CC_attr.txt	MF_attr.txt)

for file in ${attr_files[*]}
do
  if [ -s $file ]
  then
    echo "$0: File '${file}' exists and is not empty."
  else
    echo "$0: File '${file}' does not exist or is empty or is empty."
    exit 1
  fi
done

###############################################################################
# Process GO trees files to create paths and update the attributes files.
cd $root_dir/scripts

python ./GO_trees_to_paths_hbp.py

#check that the paths files exist
cd $root_dir/processed/GO_PATHS

paths_files=(BP_paths.txt CC_paths.txt MF_paths.txt)

for file in ${paths_files[*]}
do
  if [ -s $file ]
  then
    echo "$0: File '${file}' exists and is not empty."
  else
    echo "$0: File '${file}' does not exist or is empty or is empty."
    exit 1
  fi
done

# check that the updated attributes files exist

cd $root_dir/processed/GO_ATTRIBUTES_UPDATED

attr_level_files=(BP_attr_level.txt CC_attr_level.txt MF_attr_level.txt)

for file in ${attr_level_files[*]}
do
  if [ -s $file ]
  then
    echo "$0: File '${file}' exists and is not empty."
  else
    echo "$0: File '${file}' does not exist or is empty or is empty."
    exit 1
  fi
done

##############################################################################
# expand the GO term to gene annotations.
cd $root_dir/scripts

python ./expand_GO2_hbp.py

# check that the final output files 'expand.txt' exist
cd $root_dir/processed/GO_ANNOTATIONS

GO_expand_files=(BP_ALL_PC.txt CC_ALL_PC.txt MF_ALL_PC.txt BP_STRICT_PC.txt CC_STRICT_PC.txt MF_STRICT_PC.txt)

for file in ${GO_expand_files[*]}
do
  if [ -s $file ]
  then
    echo "$0: File '${file}' exists and is not empty."
  else
    echo "$0: File '${file}' does not exist or is empty or is empty."
    exit 1
  fi
done

# concatenate the GO ontology files together to make gene sets

cat BP_ALL_PC.txt CC_ALL_PC.txt MF_ALL_PC.txt > GO_ALL_PC.txt
cat BP_STRICT_PC.txt CC_STRICT_PC.txt MF_STRICT_PC.txt > GO_STRICT_PC.txt

mv GO_ALL_PC.txt $root_dir/processed/GENE_SETS/GO
mv GO_STRICT_PC.txt $root_dir/processed/GENE_SETS/GO

# ##############################################################################
#
# # Make id:gene_set: parents file for each ontology
cd $root_dir/scripts

python ./GO_gene_checker_hbp.py

# check that the final output files 'expand.txt' exist
cd $root_dir/processed/GO_GENE_CHECK

GO_id_genes_parents_files=(BP_ALL_id_genes_parents.txt CC_ALL_id_genes_parents.txt MF_ALL_id_genes_parents.txt)

for file in ${GO_id_genes_parents_files[*]}
do
  if [ -s $file ]
  then
    echo "$0: File '${file}' exists and is not empty."
  else
    echo "$0: File '${file}' does not exist or is empty."
    exit 1
  fi
done

##################################################################################
## make gene set files in magma format for GO and MGI files.

cd $root_dir/scripts

python ./MGI_gene_sets_to_Magma_hbp.py

python ./GO_gene_sets_to_Magma_hbp.py

python ./Mouse_subsets_to_Magma_hbp.py

################################################################################
# copy the processed data to the shared location.

cd $root_dir

NOW=$(date +"%d_%m_%Y-%H_%M")

# make a date stamped folder
DIR="ANNOTATION_$NOW"

# copy the processed data to the date-stamped folder

cp -a processed $DIR

cp -a logs $DIR

cp -a downloads $DIR

cp -a README $DIR

cp -a MOUSE_SUBSETS_CNS $DIR

cp annotation_err $DIR/logs

cp annotation_log $DIR/logs

# copy the datestamped folder to the shared area

cp -a $DIR /mnt/databank/ANNOTATION_AUTO_ROCKS
cp -a $DIR /home/sbijch


# remove the unwanted data at the end of the process
rm -r logs
rm -r downloads
rm -r processed
rm annotation_log
rm annotation_err
rm -r $DIR

################################################################################
