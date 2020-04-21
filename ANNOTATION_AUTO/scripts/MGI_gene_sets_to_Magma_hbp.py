#!/usr/bin/env python


import sys
import os
import re
import time
###############################################

# create file paths

# get the curent working directory
cwd = os.getcwd()

# set the root directory
root_dir = os.path.split(cwd)[0]

in_dir_path = root_dir + '/processed/GENE_SETS/MAMMALIAN_PHENOTYPE/MOUSE'

out_dir_path = root_dir + '/processed/GENE_SETS/MAMMALIAN_PHENOTYPE/MOUSE/MOUSE_MAGMA_FORMAT'

outfile_type = 'pheno_to_human_10_MAGMA.txt'

# define log file path
log_fname = root_dir + '/logs/MGI_gene_sets_to_Magma.log'

# open log file for writing
log_file = open(log_fname, 'w')

##########################################
# Main process

# generate the output file name from the input file name

for infilename in os.listdir(in_dir_path):
    count = 0
    if infilename.endswith(".txt"):

        print 'processing infile =' , infilename
        log_file.write('processing infile =  ' + infilename + '\n')

        in_fname = os.path.join(in_dir_path, infilename)

        #  make outfile names from the infile name.
        f_info = re.split('_|\.', infilename, 3)
        #print f_info
        prefix = '_'.join([f_info[0], f_info[1], f_info[2]])
        #print prefix
        outfilename = "_".join([prefix, outfile_type])
        print 'outfile =', outfilename
        log_file.write('outfile = ' + outfilename + '\n')

        out_fname = os.path.join(out_dir_path, outfilename)

        # open outfile
        outfile = open(out_fname, 'w')

##########################################

        for line in open(in_fname,'r'):
            # make each line into a list
            record = line.strip().split('\t')
            #print record
            if len(record) < 4:
                sys.exit('some records have missing data')
            else:
                name_text = '_'.join(record[1].split())
                #print name_text
                name = '_'.join([record[0],name_text])
                genes = record[2].split('|')
                if (10 <= len(genes) <= 2000):
                    gene_output = ' '.join(genes)
                    count += 1
                    outfile.write(name + '\t' + gene_output + '\n' )


        print 'number of gene sets with at least 10 genes in' , os.path.basename(in_fname),'=', count
        log_file.write('number of gene sets at least 10 genes in' + os.path.basename(in_fname) + '=' + str(count)+ '\n')

log_file.write('process complete at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')
print 'end of processing'

##########################################
