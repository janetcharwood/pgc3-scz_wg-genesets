#!/usr/bin/env python

import sys
import os
import re
import time
###############################################

# get the curent working directory
cwd = os.getcwd()

# get the root directory
root_dir = os.path.split(cwd)[0]
print root_dir

# create file paths  processed/GENE_SETS/MAMMALIAN_PHENOTYPE/MOUSE

in_file = root_dir + '/processed/GENE_SETS/MAMMALIAN_PHENOTYPE/MOUSE/MGI_single_gene_Pheno_to_human_protein_coding_gene.txt'

gene_set_path = root_dir + '/MOUSE_SUBSETS_CNS'

out_dir_path = root_dir + '/processed/GENE_SETS/MAMMALIAN_PHENOTYPE/MOUSE/MOUSE_MAGMA_FORMAT'

outfile_type =  'MAGMA_subset.txt'

log_fname = root_dir + '/logs/mouse_subsets_log.txt'

log_file = open(log_fname, 'w')

##########################################
# functions

# print a dictionary out to the screen
def print_dictionary(d):
    for k, v in d.items():
        print k, "=>", v

##########################################

# Main process

# read in the mouse pheno to human pc gene file.

log_file.write('process began at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

mouse_pheno_dict = {}
for line in open(in_file, 'r'):
    records = line.strip().split('\t')
    #print records
    # make a  dictionary of mouse pheno id and gene list
    mouse_pheno_dict[records[0]]= [records[1],records[2]]

#print_dictionary(mouse_pheno_dict)

# process the subsets of mouse phenotypes.

for infilename in os.listdir(gene_set_path):
    subset = []
    count = 0
    print 'processing infile =', infilename
    in_fname = os.path.join(gene_set_path, infilename)
    # print '', in_fname
    file_id = re.split('\.', infilename)[0]
    print 'file_id =',file_id
    log_file.write('file_id =' + file_id + '\n')
    #make outfile names from the infile name.
    outfilename = "_".join(['MGI',file_id,outfile_type])
    print 'outfile =', outfilename
    log_file.write('outfile =' + outfilename + '\n')
    out_fname = os.path.join(out_dir_path, outfilename)
    #open outfile
    outfile = open(out_fname, 'w')
    for line in open(in_fname, 'r'):
        subset_name = line.strip()
        subset.append(subset_name)
    print 'number of phenotype IDs = ', str(len(subset))
    log_file.write('number of phenotype IDs = '+ str(len(subset)) + '\n')
    # extract the subsets from the mouse pheno dictionary
    for ID in subset:
        if ID in mouse_pheno_dict.keys():
            name_text = '_'.join(mouse_pheno_dict[ID][0].split())
            #print name_text
            name = '_'.join([ID,name_text])
            #print name
            genes = mouse_pheno_dict[ID][1].replace('|'," ")
            gene_list = genes.split()
            if len(gene_list) >=20:
                count += 1
            outfile.write(name + '\t' + genes + '\n' )
    print ' number of phenotypes processed:', count
    log_file.write('number of phenotypes processed:'+ str(count) + '\n')

log_file.write('process complete at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

print 'end of processing'
