#!/usr/bin/env python

from xml.etree.ElementTree import ElementTree, Element
from re import split, findall
import time
import os

# get the curent working directory
cwd = os.getcwd()
# print 'current working directory is', cwd

# get the root directory
root_dir = os.path.split(cwd)[0]
# print 'root directory is', root_dir

# initialise variables

ontology_fname = root_dir + '/downloads/MPheno_OBO.ontology'
tree_fname = root_dir + '/processed/MP_TREES/MP_tree.txt'
attr_fname = root_dir + '/processed/MP_TREES/MP_attr.txt'

log_fname = root_dir + '/logs/mouse_obo_to_tree_log.txt'

elem_type = 'term' #name of elements to be extracted
attr_list = ['id','name','alt_id','def','is_a'] #list of element attributes

############################################
# functions

def stanza_type(line):
    type_set = set(['[Term]','[Typedef]','[Instance]'])

    for t in type_set:
        if t in line:
            return t[1:-1].lower()
    return None


def make_element(elem_type,attr_dict):

    new_elem = Element(elem_type)

    #set id
    id_set = attr_dict.get('id')

    el_id = id_set.pop()
    new_elem.set('id',el_id)
    attr_dict['id'] = el_id

    #put any additional ids into 'alt_id'
    attr_dict['alt_id'].update(id_set)
    attr_dict['alt_id'].discard(el_id)

    #process 'name' entry
    tmp = attr_dict['name'].pop()
    attr_dict['name'] = tmp

    #process 'def' entry
    tmp = attr_dict['def']
    if len(tmp) > 0:
        attr_dict['def'] = tmp.pop().split("[",1)[0].strip().strip('\"')
    else:
        attr_dict['def'] = 'None'

    #process 'is_a' entry
    tmp = set()
    for s in attr_dict['is_a']: tmp.update(findall(r'MP:\d+',s))

    attr_dict['is_a'] = tmp

    return new_elem


def element_reader(filename,elem_type,attr_list):

    attr_dict = {'stanza_type' : 'header'}
    attr_dict['is_obsolete'] = set()
    for a in attr_list: attr_dict[a] = set()

    for line in open(filename,'r'):

        s = stanza_type(line)

        #skip comment and white-space lines
        if line[0] == '!' or line.isspace():
            pass
        #if not reached end of stanza, add tag-value pairs
        elif s is None:
            tmp_tag, tmp_value = split(r'\s*:\s+',line.strip(),1)

            #only add tags in attr_list (plus is_obsolete)
            if tmp_tag in attr_list or tmp_tag == 'is_obsolete': attr_dict[tmp_tag].add(tmp_value)


        #otherwise return element and reinitialise dictionary
        else:

            #skip unwanted/obsolete elements
            if (attr_dict.pop('stanza_type') == elem_type) and 'true' not in attr_dict.pop('is_obsolete'):

                yield (make_element(elem_type,attr_dict),attr_dict)

            attr_dict = {}
            attr_dict['stanza_type'] = s
            attr_dict['is_obsolete'] = set()
            for a in attr_list: attr_dict[a] = set()

################################################################
# main process

# open log for writing
log_file = open(log_fname,'w')

##############################
# read elements into dictionary
elem_dict = {}
elem_attr={}
no_parent_list =[]
multi_parent_list =[]

for el,attr_dict in element_reader(ontology_fname,elem_type,attr_list):

    el_id = el.get('id')

    if el_id in elem_dict:

        err_str = 'Element already present: ' + el_id
        raise IOError(err_str)

    else:
        elem_dict[el_id] = el
        if len(attr_dict.get('is_a')) == 0: no_parent_list.append(el_id)
        if len(attr_dict.get('is_a')) > 1: multi_parent_list.append(el_id)
        elem_attr[el_id] = attr_dict

log_file.write('Number of elements = ' + str(len(elem_dict)) + '\n')
log_file.write('Number of nodes with no parents = ' + str(len(no_parent_list)) + '\n')
log_file.write('Number of nodes with multiple parents = ' + str(len(multi_parent_list)) + '\n\n')

#####################
# check for redundancy
redundant_ids = []

for main_id in elem_dict.keys():

    for alt_id in elem_attr[main_id]['alt_id']:

        if alt_id in elem_dict: redundant_ids.append([main_id,alt_id])

log_file.write('Number of redundant ids = ' + str(len(redundant_ids)) + '\n\n')

###################################
# set parent-child links in elements

for el_id in elem_dict.keys():
    for parent_id in elem_attr[el_id].pop('is_a'):
        elem_dict[parent_id].append(elem_dict[el_id])

###################################
# separate roots from orphan_nodes
root_list = []
orphan_list = []

for tmp_id in no_parent_list:

    if len(elem_dict[tmp_id]) > 0: root_list.append(tmp_id)
    else: orphan_list.append(tmp_id)


log_file.write('Number of root elements = ' + str(len(root_list)) + '\n')
log_file.write('Number of orphan elements = ' + str(len(orphan_list)) + '\n')

#######################################################
#output
log_file.write('\n\n##############################################\n')
##########
#save tree

MP_tree = ElementTree(elem_dict[root_list[0]])
MP_tree.write(tree_fname)

log_file.write('Tree saved to file: ' + tree_fname + '\n')
#############################################
#save attributes (only those in attr_sublist)
attr_sublist = ['name','def','alt_id']
null_str = 'None'
attr_out = open(attr_fname,'w')

#write header
attr_out.write('#id\t' + '\t'.join(attr_sublist) + '\n')

#write attributes for each id
for attr_dict in elem_attr.values():

    attr_out.write(attr_dict['id'])

    for s in attr_sublist:

        tmp_attr = attr_dict[s]
        if isinstance(tmp_attr,str): attr_out.write('\t' + tmp_attr)
        elif len(tmp_attr) == 0: attr_out.write('\t' + null_str)
        else: attr_out.write('\t' + '|'.join(tmp_attr))

    attr_out.write('\n')

attr_out.close()

log_file.write('\nSaving: ' + ', '.join(attr_sublist) + '\n')
log_file.write('Filename: ' + attr_fname)

################
# close log
log_file.close()

print 'process complete at :', time.strftime("%Y-%m-%d %H:%M")
