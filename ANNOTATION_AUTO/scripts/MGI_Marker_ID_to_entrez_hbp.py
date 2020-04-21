#!/usr/bin/python

# This script requires pandas numpy environment.

# Purpose of script: to get the entrez IDs, symbol and  gene type for all the MGI marker IDs from the NCBI_gene_info file

import sys
import os
import re

from pandas import DataFrame, Series
import pandas as pd
import gzip
import collections
import time

################

# get the curent working directory
cwd = os.getcwd()
# print 'current working directory is', cwd

# get the root directory
root_dir = os.path.split(cwd)[0]
# print 'root directory is', root_dir

# infile paths.
MGI_EntrezGene_fname = root_dir + '/downloads/MGI_EntrezGene.rpt'
NCBI_gene_info_fname = root_dir + '/downloads/gene_info.gz'

# output file paths.
log_fname = root_dir + '/logs/MGI_markerID_to_entrezID_log'

markerID_to_entrezID_ALL_fname = root_dir + '/processed/MP_ID_MAPPING/MGI_markerID_to_entrezID_ALL.txt'
markerID_to_entrezID_pc_fname = root_dir + '/processed/MP_ID_MAPPING/MGI_markerID_to_entrezID_pc.txt'

#################

# open output files for writing
log_file = open(log_fname, 'w')

out_file1 = open(markerID_to_entrezID_ALL_fname, 'w')
out_file2 = open(markerID_to_entrezID_pc_fname, 'w')

#################

# define regex.
mm_tax_ID = re.compile(r"^10090$")

#################################################################################################
# Main process:

print 'process began at :', time.strftime("%Y-%m-%d %H:%M")

log_file.write('process began at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

# open the MGI_EntrezGene.rpt.txt for reading
try:
    in_file1 = open(MGI_EntrezGene_fname, 'r')
except IOError, message:
    print >> sys.stderr, "infile1 could not be opened:"
    sys.exit(1)

log_file.write('Name of infile 1:  ' + os.path.basename(MGI_EntrezGene_fname) + '\n')

print "Processing ", os.path.basename(MGI_EntrezGene_fname)

# read the data into a data frame
MGI_markers = pd.read_table(in_file1, sep = '\t', header = None, dtype = str)
# make sure data type =string, otherwise the entrez gene IDs get decimalised!

# extract MGI marker ID, Type and Entrez gene ID colums to a new data frame
MGI_mini = MGI_markers.iloc[:, [0, 6, 8]]

# add headers to the MGI_mini dataframe
MGI_headers = ['MGI_marker_ID', 'Type', 'GeneID']
MGI_mini.columns = MGI_headers

# filter the MGI_mini dataframe for the records where the field 'Type=Gene'

MGI_gene_mini = MGI_mini[MGI_mini.Type == 'Gene']

###########################################################################################################

# Read in the NCBI_gene_info data. This file has headers. pd. table will use them by default.
# entrez gene id is independent of tax id, so should be able to get the mouse genes directly, but I have put in a tax_ID check anyway.

log_file.write('Name of infile 2:  ' + os.path.basename(NCBI_gene_info_fname) + '\n')

print 'Processing', os.path.basename(NCBI_gene_info_fname)

# open gene_info.gz into a dataframe

# read the data into a data frame using pd.table
NCBI_gene = pd.read_table(NCBI_gene_info_fname, sep = '\t', dtype = str, compression = 'gzip')
# print NCBI_gene[:10]

# remove the '#' from the #tax_id header
NCBI_gene.rename(columns={'#tax_id': 'tax_id'}, inplace=True)

# extract the  NCBI_gene columns we want
NCBI_gene_mini = NCBI_gene.iloc[:, [0, 1, 2, 4, 9]]

# filter for mouse protein coding genes in this data frame.
mouse_NCBI_pc_genes = NCBI_gene_mini[(NCBI_gene_mini.type_of_gene == 'protein-coding') & (NCBI_gene_mini.tax_id == '10090')]

mouse_NCBI_pc_genes_list = mouse_NCBI_pc_genes['GeneID'].tolist()

log_file.write('number of mouse protein coding genes in NCBI gene_info file = ' + str(len(mouse_NCBI_pc_genes_list)) + '\n')

# merge the two data frames by the NCBI_gene ID
# count how many mouse marker IDs have NCBI IDS  and compare with no of mouse genes in NCBI gene

# merge the mouse marker data frame with the NCBI gene info data frames
MGI_marker_to_entrez_ID = pd.merge(MGI_gene_mini, NCBI_gene_mini, on='GeneID')

# check that all the TAX_IDs are mouse.
# convert tax_ID column to a list

tax_ID_list = MGI_marker_to_entrez_ID['tax_id'].tolist()

for ID in tax_ID_list:
    if re.search(mm_tax_ID, ID):
        continue
    else:
        print "Error! Not all the tax_IDs are mouse!"
        sys.exit()

log_file.write('Tax IDs checked: all genes are mouse' + '\n')

# count the number of rows in the merged data frame

count = len(MGI_marker_to_entrez_ID)

# print "total number of MGI_markers mapped to entrez_IDs : ", count
log_file.write('total number of MGI_markers mapped to entrez gene IDs :  ' + str(count) + '\n')

# make dataframe with the following fields:
MGI_marker_to_entrez_ID_mini = MGI_marker_to_entrez_ID.iloc[:, [0, 1, 2, 4, 5, 6]]

#############################################################################################################

# filter for protein coding genes in the mapped data frame and count them.
pc_genes = MGI_marker_to_entrez_ID_mini[MGI_marker_to_entrez_ID_mini.type_of_gene == 'protein-coding']

# get the unique list of genes in the mapped data frame and count them.
pc_genes_list = pc_genes['GeneID'].tolist()
unique_pc_genes = set(pc_genes_list)

# print unique_pc_genes
log_file.write('number of unique protein coding genes mapped to mouse markers = ' + str(len(unique_pc_genes)) + '\n')

# print all mapped MGI markers_to_entrez_ID mappings data to outfile1
MGI_marker_to_entrez_ID_mini.to_csv(out_file1, sep='\t', index=False, header=True)

# print all mapped MGI markers_to protein coding gene entrez_ID mappings data to outfile2
pc_genes.to_csv(out_file2, sep='\t', index=False, header=True)

print 'process complete at :', time.strftime("%Y-%m-%d %H:%M")

log_file.write('process complete at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

log_file.close()


###############################################################################
# NOTES
# pd.table reads the data into a datframe and automatically uses headers as column labels.
# these are the NCBI headers:
#
# Index([u'# tax_id', u'GeneID', u'Symbol', u'LocusTag', u'Synonyms', u'dbXrefs',
#        u'chromosome', u'map_location', u'description', u'type_of_gene',
#        u'Symbol_from_nomenclature_authority',
#        u'Full_name_from_nomenclature_authority', u'Nomenclature_status',
#        u'Other_designations', u'Modification_date'],
#       dtype='object')

# MGI file MGI_EntrezGene.rpt.txt columns:

# 1.MGI Marker Accession ID
# 2.Marker Symbol
# 3.Status
# 4.Marker Name
# 5.cM Position
# 6.Chromosome
# 7.Type(Gene,DNA, Segment,Other Genome Feature,Complex/Cluster/Region,microRNA
# 8.Secondary Accession IDs(|-delimited)
# 9.Entrez Gene ID
# 10.Synonyms(|-delimited)
# 11.Feature Types(|-delimited)
# 12.Genome Coordinate Start
# 13.Genome Coordinate End
# 14.Strand
# 15.BioTypes(|-delimited)
