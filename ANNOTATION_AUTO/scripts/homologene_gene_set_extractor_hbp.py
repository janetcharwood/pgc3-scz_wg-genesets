#!/usr/bin/env python

# homologene.data is a tab delimited file containing the following columns:

# 1) HID (HomoloGene group id)
# 2) Taxonomy ID
# 3) Gene ID
# 4) Gene Symbol
# 5) Protein gi
# 6) Protein accession

import sys
import os
import re
import time
import gzip

###########################################################
#set directory paths

# get the curent working directory
cwd = os.getcwd()
# print 'current working directory is', cwd

# get the root directory
root_dir = os.path.split(cwd)[0]
# print 'root directory is', root_dir

NCBI_gene_info_fname = root_dir + '/downloads/gene_info.gz'

mouse_homologene_fname = root_dir + '/processed/HOMOLOGENE/mouse_homologene_jh.txt'
human_homologene_fname = root_dir + '/processed/HOMOLOGENE/human_homologene_jh.txt'

# outfile paths

log_fname = root_dir + '/logs/homologene_gene_set_extractor_log.txt'

human_out_all_fname = root_dir + '/processed/HOMOLOGENE/homologene_human_all_minus_pseudo.txt'
mouse_out_all_fname = root_dir + '/processed/HOMOLOGENE/homologene_mouse_all_minus_pseudo.txt'

human_out_PC_fname = root_dir + '/processed/HOMOLOGENE/homologene_human_protein_coding.txt'
mouse_out_PC_fname = root_dir + '/processed/HOMOLOGENE/homologene_mouse_protein_coding.txt'

#######################################
# open output files for writing

log_file = open(log_fname, 'w')
out_file1 = open(human_out_all_fname, 'w')
out_file2 = open(human_out_PC_fname , 'w')
out_file3 = open(mouse_out_all_fname, 'w')
out_file4 = open(mouse_out_PC_fname, 'w')

########################################

# make empty lists and dictionaries.

human_pc_genes = []
mouse_pc_genes = []

human_all_genes = []
mouse_all_genes = []

human_homol_dict ={}
mouse_homol_dict ={}

########################################
# define regular expressions

# define regex.
mouse_tax_ID = re.compile(r"^10090$")
human_tax_ID = re.compile(r"^9606$")

########################################
#set counters.

count1=0
count2=0
count3=0
count4=0
count5=0
count6=0
count7=0
count8=0

##########################################
# functions

# print a dictionary out to the screen
def print_dictionary(d):
    for k, v in d.items():
        print k, "=>", v

# make a dictionary of homologene ID GeneID, and gene symbol from the human and mouse homologene files
def make_homolgene_dict(homologene_fname, homol_dict):
    for records in open(homologene_fname, 'r'):
        record = records.strip().split('\t')
        homol_ID = record[0]# homologene ID
        gene_ID =record[2]  # entrez gene ID
        gene_symb =record[3] # entrez gene symbol
        # make a dictionary of  GeneID, homologene ID and gene symbol
        homol_dict[gene_ID]= [homol_ID,gene_symb]

def make_homologene_gene_sets(homol_dict, ncbi_gene_set, out_file):
    for gene_ID in homol_dict.keys():
        if gene_ID in ncbi_gene_set:
            output = homol_dict.get(gene_ID)
            #print '',output[0], gene_ID, output[1]
            out_file.write(output[0] + '\t' + gene_ID + '\t' + output[1]  + '\n')

##########################################

# Main program

print 'process began at :', time.strftime("%Y-%m-%d %H:%M")

log_file.write('process began at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

# make mouse and human homologene dictionaries.

make_homolgene_dict(human_homologene_fname, human_homol_dict)

#print_dictionary(human_homol_dict)

make_homolgene_dict(mouse_homologene_fname, mouse_homol_dict)

#print_dictionary(mouse_homol_dict)

##########################################

# open the NCBI 'gene_info' zip file and read it in.
# extract human and mouse gene sets
infile = gzip.open(NCBI_gene_info_fname, 'r')
records = infile.readlines()
for record in records:
    #skip the header line
    if not record.startswith("#"):
        data = record.strip().split('\t')
        #print data
        if re.match(mouse_tax_ID, data[0]):
            count1+=1
            if  data[9]=='protein-coding' :
                #print '\n', data[0],data[1],data[2],data[4],data[5],data[6],data[9]
                mouse_pc_genes.append(data[1])
                count2+=1

            if data[9]!='pseudo' :
                #print '\n', data[0],data[1],data[2],data[4],data[5],data[6],data[9]
                mouse_all_genes.append(data[1])
                count3+=1

            if data[9]=='pseudo' :
                count4+=1

        if re.match(human_tax_ID, data[0]):
            count5+=1
            if data[9]=='protein-coding' :
                #print '\n', data[0],data[1],data[2],data[4],data[5],data[6],data[9]
                human_pc_genes.append(data[1])
                count6+=1

            if data[9]!='pseudo' :
                #print '\n', data[0],data[1],data[2],data[4],data[5],data[6],data[9]
                human_all_genes.append(data[1])
                count7+=1

            if data[9]=='pseudo' :
                count8+=1


log_file.write('Total number of mouse genes =  '  + str(count1) + '\n')
log_file.write('Total number of mouse protein coding genes =  '  + str(count2) + '\n')
log_file.write('Total number of mouse all_minus_pseudo genes =  '  + str(count3) + '\n')
log_file.write('Total number of mouse pseudo genes =  '  + str(count4) + '\n')
log_file.write('Total number of human genes =  '  + str(count5) + '\n')
log_file.write('Total number of human protein coding genes =  '  + str(count6) + '\n')
log_file.write('Total number of human all_minus_pseudo genes =  '  + str(count7) + '\n')
log_file.write('Total number of human pseudo genes =  '  + str(count8) + '\n')


# make four output files using the function 'make_homologene_gene_sets':
# mouse all genes with their homologene IDs.
# human all genes with their homologene IDs.
# mouse pc genes with their homologene IDs.
# mouse pc genes with their homologene IDs.


make_homologene_gene_sets(human_homol_dict, human_all_genes, out_file1)

make_homologene_gene_sets(human_homol_dict, human_pc_genes,  out_file2)

make_homologene_gene_sets(mouse_homol_dict, mouse_all_genes, out_file3)

make_homologene_gene_sets(mouse_homol_dict, mouse_pc_genes, out_file4)


log_file.write('files: homologene_human_all_minus_pseudo.txt, homologene_human_protein_coding.txt,homologene_mouse_all_minus_pseudo.txt,homologene_mouse_protein_coding.txt created' + '\n')


log_file.write('process complete:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

print 'process complete at :', time.strftime("%Y-%m-%d %H:%M")

#############################################################################

# Used the output files from this script to run in a new version of homologene merge:
# New script: homologene_merge_protein_coding.py
