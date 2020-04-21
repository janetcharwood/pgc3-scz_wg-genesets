#!/usr/bin/env python

import sys
import os
from collections import defaultdict
import time
import re # for searching for regular expressions

# create file paths

# get the curent working directory
cwd = os.getcwd()

# set the root directory
root_dir = os.path.split(cwd)[0]

annot_dir_path = root_dir + '/processed/GO_ANNOTATIONS'

path_dir_path = root_dir + '/processed/GO_PATHS'

attr_dir_path = root_dir + '/processed/GO_ATTRIBUTES'

out_dir_path = root_dir + '/processed/GO_GENE_CHECK'

# define infile types
infile_type1 = 'attr.txt'
infile_type2 = 'ALL_PC.txt'
#infile_type3 = 'STRICT_PC.txt'

# define outfile type
outfile_type = 'ALL_id_genes_parents.txt'

# define log file
log_fname = root_dir + '/logs/GO_gene_check.log'

# open log file for writing
log_file = open(log_fname,'w')

#########################################

# functions

# print a dictionary out to the screen
def print_dictionary(d):
    for k, v in d.items():
                print k ,"=>", v

#######################################################################################################
#main process

print 'process began at :', time.strftime("%Y-%m-%d %H:%M")

log_file.write('process began at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

# process paths file and make id_parents dictionary

#open the paths.txt files
for pathfilename in os.listdir(path_dir_path):
    if pathfilename.endswith("paths.txt"):

        # make empty id parent dictionary
        id_parents = defaultdict(set)

        # make empty GO_gene dictionary
        GO_gene = defaultdict(set)

        #GO attributes dictionary
        GO_attr =defaultdict(list)

        # make a  master ID list
        master_ID = set()

        GO_ID_attr_set = set()

        # read in the path files
        file_info = re.split('_|\.',pathfilename)
        file_id = re.split('_|\.',pathfilename)[0]
        #print 'file_id =:',file_id

        # derive directory paths
        path_fname =os.path.join(path_dir_path,pathfilename)
        #print 'path infilename =', path_fname

        #get the corresponding attribute file
        attr_infilename = "_".join([file_id, infile_type1])

        # derive the attr file directory path
        attr_fname =os.path.join(attr_dir_path,attr_infilename)
        #print 'attributes file directory path = ', attr_fname

        #get the corresponding annot files
        annot_infilename = "_".join([file_id, infile_type2])
        #print annot_infilename

        # derive the annot file directory path
        annot_fname =os.path.join(annot_dir_path,annot_infilename)
        print 'annot file directory path = ', annot_fname

        # derive the outfile name and path.
        outfilename = "_".join([file_id, outfile_type])
        print 'outfilename = ',outfilename
        outfile_fname = os.path.join(out_dir_path,outfilename )
        print 'outfilepath = ',outfile_fname

        #open the output file
        GO_out = open(outfile_fname,'w')

        print 'processing paths file =',pathfilename
        log_file.write('path infile = ' + pathfilename + '\n')
        print 'attributes infilename = ',attr_infilename
        log_file.write('attributes infile = '+ attr_infilename + '\n')
        print 'annot file name = ', annot_infilename
        log_file.write('annot infile = '+ annot_infilename + '\n')

        ##############################################################################

        #read in the attributes file, extract the GO_ID and name and make the GO_attr  dictionary

        # make an empty set to put all the ids in from the attributes file.. what about alternative ids?
        id_set = set()

        for line in open(attr_fname,'r'):
            #skip the header line
            if not line.startswith('#'):
               #make each line into a list
               attr = line.strip().split('\t')
               GO_ID_attr = attr[0]
               GO_ID_name = attr[1]
               # get a list of GO ids
               id_set.add(GO_ID_attr)
               # make a dictionary from GO ID and name
               GO_attr[GO_ID_attr]=GO_ID_name

        #print_dictionary(GO_attr)

        print "GO_attr dictionary constructed from ", os.path.basename(attr_fname), 'at',time.strftime("%Y-%m-%d %H:%M")

        # #check number of unique IDs in the atrribute file
        print 'number of unique IDs in the atrribute file =  ' , len(id_set)
        log_file.write ('number of unique IDs in the atrribute file =  ' + str(len(id_set))+ '\n')

        ##############################################################################

        # process the paths file : make ID parents dictionary
        id_paths = []
        for line in open(path_fname,'r'):
            #make each line into a list
            data = line.strip().split('\t')
            id_paths.append(data)

        # make the ID parents dictionary
        for node in id_set: id_parents[node] = set([node])

        #print_dictionary(id_parents)

        for path in id_paths: id_parents[path[-1]].update(path)

        #print_dictionary(id_parents)

        #################################################################################################################

        print "id_parents dictionary constructed from ", os.path.basename(path_fname), 'at',time.strftime("%Y-%m-%d %H:%M")
        log_file.write('id_parents dictionary constructed from ' + os.path.basename(path_fname) + '\n\n')
        print 'size of id parents dictionary  =  ', len(id_parents)
        log_file.write ('size of id parents dictionary  =  ' + str(len(id_parents)) + '\n')

        #################################################################################################################
        annot_ID = set()
        annot_genes= set()
        annot_GO_gene = {}
        annot_gene_GO = defaultdict(set)

        # read in the GO_annot.txt file
        for line in open(annot_fname,'r'):
            record = line.strip().split('\t')
            #print record
            ID = record[0].split('_')[0]
            #print ID
            annot_ID.add(ID)

            #print record[0]
            gene_set = set(record[2].split('|'))
            #print gene_set

            # extract a set of all the genes in the file
            annot_genes.update(gene_set)

            #make a GO_gene dictionary
            annot_GO_gene[record[0]]=gene_set

        #print_dictionary(annot_GO_gene)
        print 'size of annot_GO_gene dictionary  =  ', len(annot_GO_gene)
        log_file.write ('size of annot_GO_gene dictionary  =  ' + str(len(annot_GO_gene)) + '\n')


        #make an output file of GO_ID,genes and all the parents?

        for GO_ID,gene_set in annot_GO_gene.items():
            parents = id_parents[GO_ID]
            # print GO_ID
            # print gene_set
            # print parents

            # make an output file: GO_annot_id_genes_parents.txt
            #GO id, genes parents.

            gene_str  = '|'.join(gene_set)
            parent_str = '|'.join(parents)
            #print '', GO_ID,  gene_str, parent_str
            #write to output file.
            GO_out.write(GO_ID + '\t' + gene_str + '\t' + parent_str  + '\n')

#check parents of a couple of terms
# print id_parents['GO:0070256']
# print '\n'
# print id_parents['GO:0070257']


print 'end of processing at : ', time.strftime("%Y-%m-%d %H:%M")

log_file.write('end of processing at :   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

#  END OF PROCESSING.
log_file.close()

########################################################################################################
