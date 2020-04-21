#!/usr/bin/python


import sys
import os
import re # for searching for regular expresssions
from xml.etree.ElementTree import ElementTree
import time

# set directory paths

# get the curent working directory
cwd = os.getcwd()

# get the root directory
root_dir = os.path.split(cwd)[0]

tree_dir_path = root_dir + '/processed/GO_TREES'
attr_dir_path = root_dir + '/processed/GO_ATTRIBUTES'

out_dir_path1 = root_dir + '/processed/GO_PATHS'
out_dir_path2 = root_dir + '/processed/GO_ATTRIBUTES_UPDATED'

# define infile type
infile_type1 = 'attr.txt'

# define outfile type
outfile_type1 = 'paths.txt'
outfile_type2 = 'attr_level.txt'

log_fname = root_dir + '/logs/GO_obo_tree_to_paths.log'

# open log in append mode
log_file = open(log_fname, 'a')

# functions

def read_attr(f_name):
# read in attributes
# f_name = file name
    ####################
    def add(add_dict, line_dict, key):
        if key in add_dict:
            err_str = 'Multiple entries for ' + key
            raise IOError(err_str)
        add_dict[key] = line_dict
    ####################
    attr_dict = {}
    in_file = open(f_name, 'r')

    #read header
    field_str = in_file.readline().lstrip('#').strip('\n').split('\t')

    line = in_file.readline().strip('\n')

    while line:

        line_dict = {}
        for tag, value in zip(field_str, line.split('\t')):
            line_dict[tag] = value

        #add blank 'level' entry
        line_dict['level'] = set()

        add(attr_dict, line_dict, line_dict['id'])

        #read next line
        line = in_file.readline().strip('\n')

    return attr_dict



##############################################################################

# main process

print 'process began at :', time.strftime("%Y-%m-%d %H:%M")

log_file.write('process began at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

# read in the files from the 'in' directory

for infilename in os.listdir(tree_dir_path):
    if infilename.endswith(".txt"):

        # read in tree, create paths
        print 'processing tree file =', infilename
        log_file.write('processing tree file=' + infilename + '\n')

        file_id = re.split('_|\.', infilename)[0]
        #print 'file_id =:',file_id

        tree_fname = os.path.join(tree_dir_path, infilename)
        #print 'tree infilename =', tree_fname

        attr_infilename = "_".join([file_id, infile_type1])
        print 'attributes infilename = ', attr_infilename
        log_file.write('attributes infilename = '+ attr_infilename + '\n')

        attr_fname = os.path.join(attr_dir_path, attr_infilename)
        #print 'attributes infilename = ', attr_fname

        GO_tree = ElementTree(file=tree_fname)

        # create dictionary containing all nodes (key) and their immediate parents (value)
        # each node has a single parent, but multiple nodes may have the same phenotype id
        node_dict = {}

        for parent in GO_tree.getiterator():
            for child in parent:
                node_dict[child] = parent

        ###############################################
        # create second dictionary containing node paths
        node_path = {}

        for node in node_dict.keys():

            #initialise node parent set
            node_path[node] = [node]

            #initialise parent node
            parent_node = node_dict[node]

            while parent_node in node_dict:

                #add parent to path
                node_path[node].append(parent_node)

                #move up one level
                tmp_node = node_dict[parent_node]
                parent_node = tmp_node

            #add root node to path
            node_path[node].append(parent_node)

       ###############################################
        # create list of paths
        id_paths = []

        for p in node_path.values():

            tmp_path = []
            for node in p:
                tmp_path.append(node.get('id'))

            tmp_path.reverse()
            id_paths.append(tmp_path)

        log_file.write('Processing paths for  ' + infilename + '\n')

        log_file.write('Number of paths = ' + str(len(id_paths)) + '\n')

        log_file.write('Maximum path length = ' + str(max([len(p) for p in id_paths])) + '\n\n')

        ###############################################
        # update attributes

        #read attribute file
        attr = read_attr(attr_fname)

        #add level annotation
        for p in id_paths:
            for level, tmp_id in enumerate(p):
                attr[tmp_id]['level'].add(str(level))

        ###############################################

        #  make outfile names from the infile name.

        outfilename1 = "_".join([file_id, outfile_type1])
        outfilename2 = "_".join([file_id, outfile_type2])

        print 'outfile1 =', outfilename1
        print 'outfile2 =', outfilename2

        path_fname = os.path.join(out_dir_path1, outfilename1)
        attr_level_fname = os.path.join(out_dir_path2, outfilename2)

        # print path_fname
        # print attr_level_fname

       ###############################################

        #save data
        #paths
        path_out = open(path_fname, 'w')

        for p in id_paths:
            path_out.write('\t'.join(p) + '\n')

        path_out.close()

        log_file.write('Paths saved to file: ' + os.path.basename(path_fname) + '\n\n')

        ###########
        #attributes
        field_str = ['id', 'name:', 'def', 'alt_id', 'level']
        attr_out = open(attr_level_fname, 'w')

        #header
        attr_out.write('#' + '\t'.join(field_str) + '\n')

        #data
        for line in attr.values():
            #print line
            #convert level annotation to string
            tmp = '|'.join(line['level'])
            line['level'] = tmp

            #save line
            attr_out.write('\t'.join([line[k] for k in field_str]) + '\n')

        attr_out.close()

        log_file.write('Updated attributes saved to file: ' + os.path.basename(attr_level_fname) + '\n')
        log_file.write('\n\n')


       ################

print 'process complete at :', time.strftime("%Y-%m-%d %H:%M")

log_file.write('\n')
log_file.write('process complete at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

#close log file
log_file.close()
