#!/usr/bin/env python

import sys
import os
from collections import defaultdict
import time


# get the curent working directory
cwd = os.getcwd()
# print 'current working directory is', cwd

# get the root directory
root_dir = os.path.split(cwd)[0]
# print 'root directory is', root_dir

# initialise infile variables
path_fname = root_dir + '/processed/MP_TREES/MP_paths.txt'
attr_fname = root_dir + '/processed/MP_TREES/MP_attr.txt'
evidence_fname = root_dir + '/processed/MP_EVIDENCE/MGI_PhenoGeno_single_protein_coding_gene.txt'

# initialise outfile variables
expand_fname = root_dir + '/processed/MP_ANNOTATIONS/MGI_single_gene_Pheno_protein_coding_annotation.txt'

log_fname = root_dir + '/logs/MP_expand_log.txt'

#######################################################

# create output files
# open  log file for writing
log_file = open(log_fname, 'w')

# open  outfile for writing
MP_expand_out = open(expand_fname, 'w')

#############################################################
# functions

# print a dictionary out to the screen
def print_dictionary(d):
    for k, v in d.items():
        print k, "=>", v

#############################################################
# make empty dictionaries, lists and sets

id_set = set()

id_paths = []

pheno_check = set()

# id parent dictionary
id_parents = defaultdict(set)

# phenotype attributes dictionary
pheno_attr =defaultdict(list)

# make a pheno gene dictionary:
pheno_geno = defaultdict(set)

# make a gene_ID, gene symbol dictionary
gene_symbol_map = defaultdict(list)

##########################################
# Main process

print 'process began at : ', time.strftime("%Y-%m-%d %H:%M")

log_file.write('process began at: ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

# read in the attributes file, extract the MP_ID and name and make the MP_attr  dictionary

for line in open(attr_fname, 'r'):
    # skip the header line
    if not line.startswith('#'):
        # make each line into a list
        attr = line.strip().split('\t')
        MP_attr_ID = attr[0]
        MP_name = attr[1]
        # get a list of MP ids
        id_set.add(MP_attr_ID)
        # make a dictionary from MP ID and name
        pheno_attr[MP_attr_ID]=MP_name

    # print_dictionary(pheno_attr)

print "MP_attr dictionary constructed from ", os.path.basename(attr_fname), 'at', time.strftime("%Y-%m-%d %H:%M")

###########

# check number of unique IDs in the atrribute file
print 'number of unique IDs in the atrribute file =  ', len(id_set)
log_file.write('number of unique IDs in the atrribute file =  ' + str(len(id_set)) + '\n')

# process the paths file : make ID parents dictionary and get a list of unique IDs

for line in open(path_fname, 'r'):
    # make each line into a list
    data = line.strip().split('\t')
    id_paths.append(data)

# make the ID parents dictionary
for node in id_set:
    id_parents[node] = set([node])

# print_dictionary(id_parents)

# update the ID parents dictionary
for path in id_paths:
    id_parents[path[-1]].update(path)

# print_dictionary(id_parents)

print "id_parents dictionary constructed from ", os.path.basename(path_fname), 'at', time.strftime("%Y-%m-%d %H:%M")
log_file.write('id_parents dictionary constructed from ' + os.path.basename(path_fname) + '\n\n')
print 'size of id parents dictionary  =  ', len(id_parents)
log_file.write('size of id parents dictionary  =  ' + str(len(id_parents)) + '\n')

######################################################################################################################

# read in the evidence file
# make a gene_symbol_map dictionary

for line in open(evidence_fname, 'r'):
    evidence = line.strip().split('\t')
    # print evidence
    gene_id = evidence[1]
    gene_symbol = evidence[2]
    pheno_id = evidence[6]
    pheno_check.add(pheno_id)
    gene_symbol_map[gene_id] = gene_symbol
    # make a phenotype: gene_ID dictionary
    pheno_geno[pheno_id].add(gene_id)

####################################################################################
# print_dictionary(pheno_geno)

# print_dictionary(gene_symbol_map)

####################################################################################

# update the pheno_geno dictionary

for ID, parents in id_parents.items():
    gene_set = pheno_geno[ID]
    for parent in parents:
        pheno_geno[parent].update(gene_set)

########################################################################################################
# make the ouput file.

for MP_id, gene_set in pheno_geno.items():
    if gene_set:
        MP_attr_name = str(pheno_attr[MP_id])
        # print MP_attr_name
        symbol_list = []
        gene_output = []
        for gene in gene_set:
            symbol = gene_symbol_map[gene]
            symbol_list.append(symbol)
            gene_output.append(gene)
        # convert the gene and symbol lists to concatenated strings
        symbol_str = str('|'.join(symbol_list))
        gene_str = str('|'.join(gene_output))
        # write to output file.
        MP_expand_out.write(MP_id + '\t' + MP_attr_name + '\t' + gene_str + '\t' + symbol_str + '\n')


########################################################################################################

# checks

print 'number of unique IDs in the pheno atrribute dictionary  = ', len(pheno_attr)
log_file.write('number of unique IDs in the pheno atrribute dictionary  = ' + str(len(pheno_attr)) + '\n')

print 'number of unique genes in evidence file: ', len(gene_symbol_map)
log_file.write('number of unique genes in evidence file: = ' + str(len(gene_symbol_map)) + '\n')

print 'number of unique phenotypes in evidence file: ', len(pheno_check)
log_file.write('number of unique phenotypes in evidence file: = ' + str(len(pheno_check)) + '\n')

print 'length of pheno_geno dictionary =  ', len(pheno_geno)
log_file.write('length of pheno_geno dictionary = ' + str(len(pheno_geno)) + '\n')


#######################################################################################################
print 'end of processing at : ', time.strftime("%Y-%m-%d %H:%M")

log_file.write('end of processing at :' + time.strftime("%Y-%m-%d %H:%M") + '\n')

#  END OF PROCESSING.
log_file.close()
