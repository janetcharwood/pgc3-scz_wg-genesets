#!/usr/bin/python


# this script extracts human data from the NCBI gene2go file.

# headers for NCBI 'gene_info' file
# tax_id	GeneID	Symbol	LocusTag	Synonyms	dbXrefs	chromosome	map_location
# description	type_of_gene	Symbol_from_nomenclature_authority	Full_name_from_nomenclature_authority
# Nomenclature_status	Other_designations	Modification_date

# Possible qualifiers are:
# set(['NOT colocalizes_with', 'contributes_to', 'NOT contributes_to', '-', 'colocalizes_with', 'NOT'])
# removed all negative qualifiers using a regex.
# NOTE: blank values in the input file are '-', where all the values in a list are '-', they have been replaced with 'None' in the output file.
# finished and checked.

import sys
import os
import re
from collections import defaultdict
import time
import gzip

# get the curent working directory
cwd = os.getcwd()
# print 'current working directory is', cwd

# get the root directory
root_dir = os.path.split(cwd)[0]
# print 'root directory is', root_dir

dir_path = root_dir + '/downloads'

# output file paths.

log_fname = root_dir + '/logs/human_gene_2go_extract_log.txt'

GO_evidence_ALL_ALL_fname = root_dir + '/processed/GO_EVIDENCE/Gene_to_GO_ALL_ev_ALL_genes_evidence.txt'
GO_evidence_ALL_PC_fname = root_dir + '/processed/GO_EVIDENCE/Gene_to_GO_ALL_ev_PC_genes_evidence.txt'

GO_evidence_STRICT_ALL_fname = root_dir + '/processed/GO_EVIDENCE/Gene_to_GO_STRICT_ev_ALL_genes_evidence.txt'
GO_evidence_STRICT_PC_fname = root_dir + '/processed/GO_EVIDENCE/Gene_to_GO_STRICT_ev_PC_genes_evidence.txt'

#############################

# open output files for writing
log_file = open(log_fname, 'w')

out_file1 = open(GO_evidence_ALL_ALL_fname, 'w')
out_file2 = open(GO_evidence_ALL_PC_fname, 'w')
out_file3 = open(GO_evidence_STRICT_ALL_fname, 'w')
out_file4 = open(GO_evidence_STRICT_PC_fname, 'w')

#############################

# define regex
# Remove the following qualifiers:['NOT', 'NOT colocalizes_with', 'NOT contributes_to'])
neg_qualifier = re.compile("^NOT\s|^NOT")

# regex to get human entries using the tax id.
tax_ID = re.compile(r"^9606$")

# make a list of ALL evidence codes
evidence_codes_all = ['EXP', 'IDA', 'IPI', 'IEP', 'IGI', 'IMP', 'ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'IBA', 'IKR', 'IBD', 'IRD', 'TAS', 'ND', 'IC', 'IEA', 'NAS', 'RCA']

# make lists of evidence codes excluding: IEA, NAS, RCA.
evidence_codes_strict = ['EXP', 'IDA', 'IPI', 'IEP', 'IGI', 'IMP', 'ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'IBA', 'IKR', 'IBD', 'IRD', 'TAS', 'ND', 'IC']

# make sets of genes.
genes_ALL = set()
genes_STRICT = set()

# make a list of lists for the records filtered by species and attributes.
records_ALL = []

# make a list of lists for filtered records .
records_STRICT = []

# make a set to of GO term gene types in the genes_to_go file.
GO_gene_types = set()
filtered_qualifiers = set()

# make gene-info dictionary for all genes.
gene_info_map = {}

# make gene-GO dictionaries
gene2go_dict_all = defaultdict(list)
gene2go_dict_strict = defaultdict(list)

###########################################################################
def print_dictionary(test_dict):
    for k, v in test_dict.items():
        print k, "=>", v

###########################################################################
# get the gene info for the selected genes from the gene_info map and make a new dictionary
def get_gene_info(gene_set, gene2go_dict):
    for gene_ID in gene_info_map.keys():
        if gene_ID in gene_set:
            gene_info = gene_info_map.get(gene_ID)
            # construct gene2go_dict
            gene2go_dict[gene_ID] = gene_info
############################################################################
# function to process data for output
def process_genes(gene_ID, record_subset):
    ID = set()
    symb = set()
    GO = []
    Evidence_codes = []
    Qual = []
    Pubmed = []
    for i in range(0, len(record_subset)):
        data = record_subset[i]
        if gene_ID == data[1]:
            ID.add(data[1])
            symb.add(gene_info[0])
            GO.append(data[2])
            Evidence_codes.append(data[3])
            Qual.append(data[4])
            Pubmed.append(data[6])

    ID_string = ''.join(ID)
    symb_string = ''.join(symb)
    GO_string = "|".join(GO)
    Evidence_string = "|".join(Evidence_codes)

    # replace lists without any values for the qualifiers with 'None',
    # otherwise join the lists into pipe delimited strings.
    if all(x == '-' for x in Qual):
        Qual_string = 'None'
    else:
        Qual_string = "|".join(Qual)

    # replace lists without any values for the Pubmed Ids with 'None',
    # otherwise join the lists into pipe delimited strings.
    if all(y == '-' for y in Pubmed):
        Pubmed_string = 'None'
    else:
        Pubmed_string = ";".join(Pubmed)

    str_output1 = str(ID_string + '\t' + symb_string + '\t' + GO_string + '\t' + Evidence_string + '\t' + Qual_string + '\t' + Pubmed_string + '\n')
    return str_output1
###########################################################################
# select evidence codes
def select_evidence(evidence_codes):
    for record in records_ALL:
        if record[3] in evidence_codes_strict:
            # print record
            # make another subset of records here for the evidence codes selected
            records_STRICT.append(record)
            # get a list of unique GENE_IDs in this subset
            genes_STRICT.add(record[1])

###########################################################################
# MAIN PROCESS:

print 'process began at :', time.strftime("%Y-%m-%d %H:%M")

log_file.write('process began at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

# change directory to the path above
os.chdir(dir_path)

# make a gene info dictionary
# the below used for testing
# infile = open ('gene_info', 'r')

# open the NCBI 'gene_info' zip file and read it in.
infile = gzip.open('gene_info.gz', 'r')
records = infile.readlines()
for record in records:
    # skip the header line
    if not record.startswith("#"):
        gene_data = record.strip().split('\t')
        # print gene_data
        if re.match(tax_ID, gene_data[0]):
            gene_ID = gene_data[1]
            gene_symbol = gene_data[2]
            gene_type = gene_data[9]
            # print '' , gene_data[0], gene_ID, gene_symbol, gene_type
            # make a gene-info map
            gene_info_map[gene_ID] = [gene_symbol, gene_type]

print 'made gene info dictionary'

# print the gene info dictionary to check it.
# print_dictionary(gene_info_map)

#############################################################################
# process the gene2go.gz file, get a unique list of genes from this file
# and make a lists of lists (records) of the data.

infile2 = gzip.open('gene2go.gz', 'r')
lines = infile2.readlines()
for line in lines:
    # skip the header line
    if not line.startswith("#"):
        gene_data = line.strip().split('\t')
        # print gene_data
        # get human records  without 'NOT' in the attributes line.
        if re.match(tax_ID, gene_data[0]) and not re.match(neg_qualifier, gene_data[4]):
            # print gene_data
            filtered_qualifiers.add(gene_data[4])
            # derive a subset of human GO records without 'NOT'
            # in the attributes line.
            records_ALL.append(gene_data)
            # get a list of unique GENE_IDs from this subset of records
            genes_ALL.add(gene_data[1])

# get the gene info for the 'ALL' evidence codes set of records
# and make 'gene2go_dict_all' dictionary.

get_gene_info(genes_ALL, gene2go_dict_all)

# print out the gene2go_dict_all dictionary to check it.
# print_dictionary(gene2go_dict_all)

print 'length of genes_ALL = ', len(genes_ALL)

log_file.write('length of genes_ALL = ' + str(len(genes_ALL)) + '\n')

print 'length of gene2go_dict_all = ', len(gene2go_dict_all)
log_file.write('length of gene2go_dict_all = ' + str(len(gene2go_dict_all)) + '\n')

# make output file for all evidence codes and all gene types

for gene_ID, gene_info in gene2go_dict_all.items():
    # print(process_genes(gene_ID, records_ALL))
    # call the process genes function and write output to outfile1.
    out_file1.write(process_genes(gene_ID, records_ALL))

log_file.write('constructed file' + os.path.basename(GO_evidence_ALL_ALL_fname) + '\n')


# make output file for all evidence codes and protein coding genes only.
for gene_ID, gene_info in gene2go_dict_all.items():
    if gene_info[1] == 'protein-coding':
        # print gene_ID, "=>", gene_info
        # print(process_genes(gene_ID, records_ALL))
        # call the process genes function and write output to outfile2.
        out_file2.write(process_genes(gene_ID, records_ALL))

log_file.write('constructed file' + os.path.basename(GO_evidence_ALL_PC_fname) + '\n')

#################################################################################################
# Derive  records for STRICT evidence codes using 'select_evidence' function

select_evidence(evidence_codes_strict)

print 'number of strict records = ', len(records_STRICT)
print 'number of strict genes = ', len(genes_STRICT)

# get the gene info for the 'STRICT' evidence codes set of records and make 'gene2go_dict_strict' dictionary.

get_gene_info(genes_STRICT, gene2go_dict_strict)

# print out the gene2go_dict_strict dictionary to check it.
# print_dictionary(gene2go_dict_strict)

print 'length of genes_STRICT = ', len(genes_STRICT)
print 'length of gene2go_dict_strict = ', len(gene2go_dict_strict)

# make output file for all evidence codes and all gene types

for gene_ID, gene_info in gene2go_dict_strict.items():
    # print(process_genes(gene_ID, records_ALL))
    # call the process genes function and write output to outfile3.
    out_file3.write(process_genes(gene_ID, records_STRICT))

log_file.write('constructed file' + os.path.basename(GO_evidence_STRICT_ALL_fname) + '\n')

# make output file for all evidence codes and protein coding genes only.
for gene_ID, gene_info in gene2go_dict_strict.items():
    if gene_info[1] == 'protein-coding':
        # print(process_genes(gene_ID, records_ALL))
        # call the process genes function and write output to outfile4.
        out_file4.write(process_genes(gene_ID, records_STRICT))

log_file.write('constructed file' + os.path.basename(GO_evidence_STRICT_PC_fname) + '\n')

##############################################################################################################

log_file.write('process complete at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

log_file.close()
print 'process complete at :', time.strftime("%Y-%m-%d %H:%M")
