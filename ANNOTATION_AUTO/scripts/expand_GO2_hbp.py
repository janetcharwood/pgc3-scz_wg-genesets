#!/usr/bin/env python

import sys
import os
from collections import defaultdict
import time
import re

# create file paths

# get the curent working directory
cwd = os.getcwd()

# get the root directory
root_dir = os.path.split(cwd)[0]

evidence_fname1 = root_dir + '/processed/GO_EVIDENCE/Gene_to_GO_ALL_ev_PC_genes_evidence.txt'
evidence_fname2 = root_dir + '/processed/GO_EVIDENCE/Gene_to_GO_STRICT_ev_PC_genes_evidence.txt'
evidence_fname3 = root_dir + '/processed/GO_EVIDENCE/Gene_to_GO_ALL_ev_ALL_genes_evidence.txt'
evidence_fname4 = root_dir + '/processed/GO_EVIDENCE/Gene_to_GO_STRICT_ev_ALL_genes_evidence.txt'

path_dir_path = root_dir + '/processed/GO_PATHS'
attr_dir_path = root_dir + '/processed/GO_ATTRIBUTES'
out_dir_path = root_dir + '/processed/GO_ANNOTATIONS'

########################################

# define infile type
infile_type1 = 'attr.txt'

# define outfile types
outfile_type1 = 'ALL_PC.txt'
outfile_type2 = 'STRICT_PC.txt'
outfile_type3 = 'ALL_ALL.txt'
outfile_type4 = 'STRICT_ALL.txt'

# define log file
log_fname = root_dir + '/logs/expand_GO.log'

# open log file for writing
log_file = open(log_fname, 'w')

#############################################################

# make empty dictionaries, lists and sets

gene_symbol_map_ALL_PC = defaultdict(list)
gene_symbol_map_STRICT_PC = defaultdict(list)
gene_symbol_map_ALL_ALL = defaultdict(list)
gene_symbol_map_STRICT_ALL = defaultdict(list)


gene_GO_ALL_PC = defaultdict(set)
gene_GO_STRICT_PC = defaultdict(set)
gene_GO_ALL_ALL = defaultdict(set)
gene_GO_STRICT_ALL = defaultdict(set)

GO_Master_list_ALL_PC = set()
GO_Master_list_STRICT_PC = set()
GO_Master_list_ALL_ALL = set()
GO_Master_list_STRICT_ALL = set()

##########################################
# functions

# print a dictionary out to the screen
def print_dictionary(d):
    for k, v in d.items():
        print k, "=>", v


# read in the evidence file and make a gene_symbol_map dictionary
# and a gene to GO dictionary
def read_evidence(evidence_file,GO_Master_list,gene_symbol_map,gene_GO):
    for line in open(evidence_file, 'r'):
        evidence = line.strip().split('\t')
        gene_id = evidence[0]
        gene_symbol = evidence[1]
        # split the GO ids into a set
        GO_set = set(evidence[2].split("|"))

        # add each GO_id to GO_Master_list to get set of unique GO_IDs
        # from the evidence file.
        for id in GO_set:
            GO_Master_list.add(id)

        # make a gene-symbol map dictionary
        gene_symbol_map[gene_id] = gene_symbol

        # make a gene_GO dictionary
        gene_GO[gene_id] = GO_set

#############################################################################
def process_ontologies (outfile_type,gene_GO, gene_symbol_map):
    for pathfilename in os.listdir(path_dir_path):
        if pathfilename.endswith("paths.txt"):

            # make empty id parent dictionary
            id_parents = defaultdict(set)

            # make empty GO_gene dictionary
            GO_gene = defaultdict(set)

            #GO attributes dictionary
            GO_attr = defaultdict(list)

            # make a  master ID list
            master_ID = set()

            GO_ID_attr_set = set()

            # read in the path files
            file_id = re.split('_|\.', pathfilename)[0]
            #print 'file_id =:',file_id

            # derive directory paths
            path_fname = os.path.join(path_dir_path, pathfilename)
            #print 'path infilename =', path_fname

            #get the corresponding attribute file
            attr_infilename = "_".join([file_id, infile_type1])

            # derive the attr file directory path
            attr_fname = os.path.join(attr_dir_path, attr_infilename)
            #print 'attributes file directory path = ', attr_fname

            # derive the outfile name and path.
            outfilename = "_".join([file_id, outfile_type])
            #print 'outfilename = ', outfilename
            outfile_fname = os.path.join(out_dir_path, outfilename)
            #print 'outfilepath = ', outfile_fname

            #open the output file
            GO_expand_out = open(outfile_fname, 'w')

            #print 'processing paths file =', pathfilename
            # log_file.write('path infile = ' + pathfilename + '\n')
            #print 'attributes infilename = ', attr_infilename
            # log_file.write('attributes infile = '+ attr_infilename + '\n')

            ###################################################################

            # read in the attributes file, extract the GO_ID and name
            # and make the GO_attr  dictionary

            id_set = set()

            for line in open(attr_fname, 'r'):
                #skip the header line
                if not line.startswith('#'):
                    #make each line into a list
                    attr = line.strip().split('\t')
                    GO_ID_attr = attr[0]
                    GO_ID_name = attr[1]
                    # get a list of GO ids
                    id_set.add(GO_ID_attr)
                    # make a dictionary from GO ID and name
                    GO_attr[GO_ID_attr] = GO_ID_name

            # print_dictionary(GO_attr)

            # print "GO_attr dictionary constructed from ", os.path.basename(attr_fname), 'at', time.strftime("%Y-%m-%d %H:%M")

            # #check number of unique IDs in the atrribute file
            # print 'number of unique IDs in the atrribute file =  ', len(id_set)
            # log_file.write('number of unique IDs in the atrribute file =  ' + str(len(id_set)) + '\n')

            #####################################################################
            # process the paths file : make ID parents dictionary and get a list of unique IDs
            id_paths = []
            for line in open(path_fname, 'r'):
                #make each line into a list
                data = line.strip().split('\t')
                id_paths.append(data)

            # make the ID parents dictionary
            for node in id_set: id_parents[node] = set([node])

            # print_dictionary(id_parents)

            for path in id_paths: id_parents[path[-1]].update(path)

            # print_dictionary(id_parents)

            ####################################################################

            # print "id_parents dictionary constructed from ", os.path.basename(path_fname), 'at', time.strftime("%Y-%m-%d %H:%M")
            # log_file.write('id_parents dictionary constructed from ' + os.path.basename(path_fname) + '\n\n')
            # print 'size of id parents dictionary  =  ', len(id_parents)
            # log_file.write('size of id parents dictionary  =  ' + str(len(id_parents)) + '\n')

            ####################################################################

            # make a GO_gene dictionary from the gene-GO dictionary

            # print 'constructing GO_gene dictionary'

            for GO_id in id_set:
                for gene_id, GO_set in gene_GO.items():
                    if GO_id in GO_set:
                    # make a GO_gene dictionary
                        GO_gene[GO_id].add(gene_id)

            # print_dictionary(GO_gene)

            # print ' GO_gene dictionary constructed'

            # print 'length of GO_gene dictionary BEFORE update =  ', len(GO_gene)
            # log_file.write('length of GO_gene dictionary BEFORE update =  ' + str(len(GO_gene)) + '\n')

            ################################################################
            # update the GO-gene dictionary

            for ID, parents in id_parents.items():
                gene_set = GO_gene[ID]
                for parent in parents:
                    GO_gene[parent].update(gene_set)

            # print 'checking GO_gene dictionary entry AFTER updating with parent-gene mappings.'

            # print 'GO_gene dictionary upated at :', time.strftime("%Y-%m-%d %H:%M")

            ###################################################################
            # line 257 needed to remove any ids without genes.
            line_count = 0
            for GO_id, gene_set in GO_gene.items():
                if gene_set:
                    GO_ID_name = GO_attr[GO_id]
                    symbol_list = []
                    gene_output = []
                    for gene in gene_set:
                        symbol = gene_symbol_map[gene]
                        symbol_list.append(symbol)
                        gene_output.append(gene)
                    # convert the gene and symbol lists to concatenated strings
                    symbol_str = '|'.join(symbol_list)
                    gene_str = '|'.join(gene_output)
                    # check the output with a print statement
                    #print '', GO_id,  GO_ID_name,  gene_str,  symbol_str
                    #write to output file.
                    GO_expand_out.write(GO_id + '\t' + GO_ID_name + '\t' + gene_str + '\t' + symbol_str + '\n')
                    line_count += 1

            # print 'processed  ', os.path.basename(path_fname)
            # print 'number of lines output to file =    ', line_count
            # log_file.write('output'  + str(line_count) + 'lines to ' + os.path.basename(path_fname) + '\n\n')
            # log_file.write('processed ' + os.path.basename(path_fname) + 'at' + time.strftime("%Y-%m-%d %H:%M") + '\n\n')

#############################################################################
#  check an entry in the GO_gene dictionary
def parent_checker():
    for x in parent_check:
        for GO_id, gene_set in GO_gene.items():
            if x == GO_id:
                print GO_id, "=>", gene_set

############################################################################
#main process

print 'process began at :', time.strftime("%Y-%m-%d %H:%M")

log_file.write('process began at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

######################################################################################################

# process evidence file 1

log_file.write('processing   ' + os.path.basename(evidence_fname1) + '\n')

read_evidence(evidence_fname1, GO_Master_list_ALL_PC, gene_symbol_map_ALL_PC, gene_GO_ALL_PC)

log_file.write('number of unique genes in evidence file: '+ str(len(gene_symbol_map_ALL_PC)) + '\n')

# check length of GO_Master_list
# print 'number of unique GO ids in ', os.path.basename(evidence_fname1) , '=', len(GO_Master_list_ALL_PC)

log_file.write('number of unique GO ids in' + os.path.basename(evidence_fname1) + '=' + str(len(GO_Master_list_ALL_PC))  + '\n')

process_ontologies(outfile_type1,gene_GO_ALL_PC, gene_symbol_map_ALL_PC)

######################################################################################################

# process evidence file 2

log_file.write('processing   ' + os.path.basename(evidence_fname2) + '\n')

read_evidence(evidence_fname2, GO_Master_list_STRICT_PC, gene_symbol_map_STRICT_PC, gene_GO_STRICT_PC)

log_file.write('number of unique genes in STRICT_PC evidence file:  ' + str(len(gene_symbol_map_STRICT_PC)) + '\n')

# check length of GO_Master_list

log_file.write('number of unique GO ids in   ' + os.path.basename(evidence_fname2) + '=' + str(len(GO_Master_list_STRICT_PC)) + '\n')

process_ontologies (outfile_type2,gene_GO_STRICT_PC, gene_symbol_map_STRICT_PC)

#########################################################################################################
# process evidence file 3

# log_file.write('processing   ' + os.path.basename(evidence_fname3) + '\n')
#
# read_evidence(evidence_fname3, GO_Master_list_ALL_ALL, gene_symbol_map_ALL_ALL, gene_GO_ALL_ALL)
#
# log_file.write('number of unique genes in ALL_ALL evidence file:  ' + str(len(gene_symbol_map_ALL_ALL)) + '\n')
#
# # check length of GO_Master_list
#
# log_file.write('number of unique GO ids in   ' + os.path.basename(evidence_fname3) + '=' + str(len(GO_Master_list_ALL_ALL)) + '\n')
#
# process_ontologies (outfile_type3,gene_GO_ALL_ALL, gene_symbol_map_ALL_ALL)

#########################################################################################################

# process evidence file 4

# log_file.write('processing   ' + os.path.basename(evidence_fname4) + '\n')
#
# read_evidence(evidence_fname4, GO_Master_list_STRICT_ALL, gene_symbol_map_STRICT_ALL, gene_GO_STRICT_ALL)
#
# log_file.write('number of unique genes in STRICT_ALL evidence file:  ' + str(len(gene_symbol_map_STRICT_ALL)) + '\n')
#
# # check length of GO_Master_list
#
# log_file.write('number of unique GO ids in  ' + os.path.basename(evidence_fname4) + '=' +  str(len(GO_Master_list_STRICT_ALL)) + '\n')
#
# process_ontologies(outfile_type4, gene_GO_STRICT_ALL, gene_symbol_map_STRICT_ALL)

#########################################################################################################

print 'end of processing at : ', time.strftime("%Y-%m-%d %H:%M")

log_file.write('end of processing at :   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

#END OF PROCESSING.
log_file.close()
