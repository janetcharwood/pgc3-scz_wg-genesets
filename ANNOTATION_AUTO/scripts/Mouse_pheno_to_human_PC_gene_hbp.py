#!/usr/bin/env python

import sys
import os
import time

#############################################################

# create file paths
# get the curent working directory
cwd = os.getcwd()
# print 'current working directory is', cwd

# get the root directory
root_dir = os.path.split(cwd)[0]
# print 'root directory is', root_dir

# create infile paths.

f_name1 = root_dir + '/processed/HOMOLOGENE/hm_one_to_one_homol_PC.txt'

f_name2 = root_dir + '/processed/MP_ANNOTATIONS/MGI_single_gene_Pheno_protein_coding_annotation.txt'

# define log file path
log_fname = root_dir + '/logs/mouse_to_human_genes.log'

# define outfilepath
out_fname = root_dir + '/processed/GENE_SETS/MAMMALIAN_PHENOTYPE/MOUSE/MGI_single_gene_Pheno_to_human_protein_coding_gene.txt'

# open output file
outfile = open(out_fname, 'w')

# open log file for writing
log_file = open(log_fname, 'w')

#############################################################
# make empty dictionary
M_HU_gene_map = {}

#############################################################
# functions
# print a dictionary out to the screen
def print_dictionary(d):
    for k, v in d.items():
        print k, "=>", v

#############################################################

# main process

print 'process began at :', time.strftime("%Y-%m-%d %H:%M")

log_file.write('process began at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

# make a dictionary of mouse to human IDs from the homologene file
# open hm_one_to_one_homol_PC.txt file

for record in open(f_name1,'r'):
    # make each line into a list
    records = record.strip().split('\t')
    #print records
    h_geneID=records[1]
    h_symb=records[2]
    m_geneID=records[3]
    m_symb=records[4]
    #make a dictionary using the mouse geneID as the key
    M_HU_gene_map[m_geneID] = [m_symb,h_geneID,h_symb]

# print the dictionary to check it.
# print_dictionary(M_HU_gene_map)


#read in MGI_single_gene_Pheno_protein_coding_annotation.txt

for line in open(f_name2,'r'):
    lines = line.strip().split('\t')
    # print lines
    # split the gene IDs and gene symbols up into a list:
    mouse_genes = lines[2].split('|')
    #print mouse_genes
    # iterate through the list of mouse genes and get the human gene and its symbol from the M_HU_gene_map, use lists to keep the order
    hu_gene_id_list =[]
    hu_gene_sym_list =[]
    for ID in mouse_genes:
        if ID in M_HU_gene_map.keys():
            homol_data = M_HU_gene_map.get(ID)
            # print '\n'
            # print ID
            # print homol_data
            # print homol_data[1] # human gene ID
            # print homol_data[2] # human gene symbol
            hu_gene_id_list.append(homol_data[1])
            hu_gene_sym_list.append(homol_data[2])

    # check that there is at least one human homologue by checking that the list is not empty, before making output file.
    if hu_gene_id_list :
        #convert the human gene id and gene symbol lists to strings
        hu_symbol_str = '|'.join(hu_gene_id_list)
        hu_gene_str  = '|'.join(hu_gene_sym_list)
        # check output format: human
        # print '\n',lines[0], lines[1],hu_symbol_str,hu_gene_str
        # Make output file , format  = MP ID: pheno name:  human genes: human symbols.
        outfile.write(lines[0] + '\t' + lines[1]  + '\t' + hu_symbol_str + '\t' +  hu_gene_str + '\n')


log_file.write('process complete at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')
print 'process complete at :', time.strftime("%Y-%m-%d %H:%M")
