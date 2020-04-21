#!/usr/bin/python

# Empty values in the lists are replaced with 'None'

import sys
import os
import re
from pandas import DataFrame, Series
import pandas as pd
import collections
import time

# get the curent working directory
cwd = os.getcwd()
#print 'current working directory is', cwd

# get the root directory
root_dir = os.path.split(cwd)[0]
#print 'root directory is', root_dir

in_file1_fname = root_dir + '/downloads/MGI_PhenoGenoMP.rpt'
in_file2_fname = root_dir + '/processed/MP_ID_MAPPING/MGI_markerID_to_entrezID_ALL.txt'

# open output file paths.
log_fname = root_dir + '/logs/MGI_evidence_log.txt'

MGI_PhenoGeno_ALL_fname = root_dir + '/processed/MP_EVIDENCE/MGI_PhenoGeno_single_gene_ALL.txt'
MGI_PhenoGeno_PC_fname = root_dir + '/processed/MP_EVIDENCE/MGI_PhenoGeno_single_protein_coding_gene.txt'

# open log file for writing

log_file = open(log_fname, 'w')

# open  outfiles for writing

outfile1 = open(MGI_PhenoGeno_ALL_fname, 'w')
outfile2 = open(MGI_PhenoGeno_PC_fname, 'w')

##########
# make empty lists
master_marker_list = []
single_gene = []

# set counters
count1 = 0
count2 = 0
single_gene_count = 0

#################################################################################################
# Main process

print 'process began at :', time.strftime("%Y-%m-%d %H:%M")

log_file.write('process began at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

# open the mouse phenotype infile for reading
for line in open(in_file1_fname, 'r'):
    # skip the header line
    if not line.startswith('#'):
        count1 += 1
        # make each line into a list
        data = line.strip().split('\t')
        # exclude transgenes
        if not data[0].startswith("Tg"):
            # Get single gene/phenotype relationships
            # get the MGI marker ID field , split data[5] by comma-delimiter and make into a list
            MGI_marker_list = data[5].split(',')
            if len(MGI_marker_list) == 1:
                single_gene.append(data)
                single_gene_count += 1
                # add these data to a list of lists, then convert to a dataframe.

# convert 'single_gene' list of lists to a data frame

single_gene_df = pd.DataFrame(single_gene,  columns=['Allelic_Composition', 'Allelic_Symbol', 'Genetic_Background', 'Pheno_ID', 'PubMed_ID', 'MGI_marker_ID'])
# print single_gene_df[:10]

# open the MGI_marker_ID_to_entrez_ID_ALL.txt file for reading

try:
    in_file2 = open(in_file2_fname, 'r')
except IOError, message:
    print >> sys.stderr, "infile2 could not be opened:"
    sys.exit(1)

print 'processing gene ID data'


# read in the MGI_marker_ID_to_entrez_ID_ALL.txt file as a dataframe
MGI_marker_ID_to_entrez_ID_df = pd.read_table(in_file2, sep='\t')

# merge the two dataframes by the MGI marker ID.
mouse_single_gene_pheno_to_entrez_gene_ID_df = pd.merge(MGI_marker_ID_to_entrez_ID_df, single_gene_df, on ='MGI_marker_ID')
# print mouse_single_gene_pheno_to_entrez_gene_ID_df

mini_df = mouse_single_gene_pheno_to_entrez_gene_ID_df.iloc[:, [0, 2, 3, 5, 6, 8, 9, 10]]
# print mini_df

# convert columns 2 and 9 to lists
# zip them
# add the new list to the data frame

############################
Geno = mini_df['GeneID'].tolist()
Pheno = mini_df['Pheno_ID'].tolist()

Ge_Phe = map(lambda Geno, Pheno: str(Geno) + str(Pheno), Geno, Pheno)
# print Ge_Phe

# make the list into a series
se = pd.Series(Ge_Phe)
# add the series to the dataframe
mini_df.insert(0, 'Geno_Pheno', se.values)

# get a list of unique 'Geno_Pheno' IDs
ID = list(set(Ge_Phe))
# print ID
#############################

print 'merging Gene and Phenotype data'

for y in ID:
    tmp = pd.DataFrame(mini_df[mini_df.Geno_Pheno == y])
    # print tmp
    output = []
    for i in range(1, 9):
        values = list(set(tmp.iloc[:, i]))
        # print values
        # replace empty values with 'None' (see the pubmed ID column, not all evidence has a pubmed ID)
        values = ['None' if not x else x for x in values]
        # convert data to string and join values by a pipe
        str_values = '|'.join(map(str, values))
        #print str_values
        output.append(str_values)
    # print output
    str_output = str('\t'.join(output))
    outfile1.write (str_output + '\n')

    if output[3] == 'protein-coding':
        count2 += 1
        pc_str_output = str('\t'.join(output))
        outfile2.write(pc_str_output + '\n')

print 'Gene and Phenotype data merged at :', time.strftime("%Y-%m-%d %H:%M")

print "total number of records =  ", str(count1)
log_file.write('number of records    ' + str(count1) + '\n')

print "number of single gene records excluding transgenes =  ", str(single_gene_count)
log_file.write('number of single gene-pheno records excluding transgenes   ' + str(single_gene_count) + '\n')

print "number of single protein coding gene records  =  ", str(count2)
log_file.write('number of single protein coding gene-pheno records   ' + str(count2) + '\n')

# close the out files
outfile1.close()
outfile2.close()

log_file.write('process complete at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

log_file.close()

print 'process complete at :', time.strftime("%Y-%m-%d %H:%M")
