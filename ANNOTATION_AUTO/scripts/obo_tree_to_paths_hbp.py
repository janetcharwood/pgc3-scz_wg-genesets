#!/usr/bin/env python

from xml.etree.ElementTree import ElementTree
import time
import os

# get the curent working directory
cwd = os.getcwd()
# print 'current working directory is', cwd

# get the root directory
root_dir = os.path.split(cwd)[0]
# print 'root directory is', root_dir

# initialise file variables
tree_fname = root_dir + '/processed/MP_TREES/MP_tree.txt'
attr_fname = root_dir + '/processed/MP_TREES/MP_attr.txt'

path_fname = root_dir + '/processed/MP_TREES/MP_paths.txt'
attr_level_fname = root_dir + '/processed/MP_TREES/MP_attr_level.txt'

log_fname = root_dir + '/logs/mouse_obo_tree_to_paths_log.txt'

####################
# read in attributes
# f_name = file name
def read_attr(f_name):


    ####################
    def add(add_dict,line_dict,key):

        if key in add_dict:

            err_str = 'Multiple entries for ' + key
            raise IOError(err_str)

        add_dict[key] = line_dict
    ####################
    attr_dict = {}
    in_file = open(f_name,'r')

    #read header
    field_str = in_file.readline().lstrip('#').strip('\n').split('\t')

    line = in_file.readline().strip('\n')

    while line:

        line_dict = {}
        for tag,value in zip(field_str,line.split('\t')): line_dict[tag] = value

        #add blank 'level' entry
        line_dict['level'] = set()

        add(attr_dict,line_dict,line_dict['id'])

        #read next line
        line = in_file.readline().strip('\n')

    return attr_dict


#############################
# open log for writing
log_file = open(log_fname,'w')
#############################

MP_tree = ElementTree(file = tree_fname)

#############################
# read in tree, create paths
###############################################
# create dictionary containing all nodes (key) and their immediate parents (value)
# each node has a single parent, but multiple nodes may have the same phenotype id
node_dict = {}

for parent in MP_tree.getiterator():
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

log_file.write('Number of paths = ' + str(len(id_paths)) + '\n')

log_file.write('Maximum path length = ' + str(max([len(p) for p in id_paths])) + '\n\n')
###############################################
# update attributes

#read attribute file
attr = read_attr(attr_fname)

#add level annotation
for p in id_paths:

    for level,tmp_id in enumerate(p): attr[tmp_id]['level'].add(str(level))


###############################################
# save data

######
#paths
path_out = open(path_fname,'w')

for p in id_paths: path_out.write('\t'.join(p) + '\n')

path_out.close()

log_file.write('Paths saved to file: ' + path_fname + '\n\n')

###########
#attributes
field_str = ['id','name','def','alt_id','level']
attr_out = open(attr_level_fname,'w')

#header
attr_out.write('#' + '\t'.join(field_str) + '\n')

#data
for line in attr.values():

    #convert level annotation to string
    tmp = '|'.join(line['level'])
    line['level'] = tmp

    #save line
    attr_out.write('\t'.join([line[k] for k in field_str]) + '\n')

attr_out.close()

log_file.write('Updated attributes saved to file: ' + attr_level_fname + '\n')
################
# close log
log_file.close()

print 'process complete at :', time.strftime("%Y-%m-%d %H:%M")
