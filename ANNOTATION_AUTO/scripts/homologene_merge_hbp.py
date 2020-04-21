#!/usr/bin/python

import sys
import os
import time
from pandas import DataFrame, Series
import pandas as pd
import collections

########################################################################
# get the curent working directory
cwd = os.getcwd()
# print 'current working directory is', cwd

# get the root directory
root_dir = os.path.split(cwd)[0]
# print 'root directory is', root_dir

indir_path = root_dir + '/processed/HOMOLOGENE'
# define infiles

infname_1 = root_dir + '/processed/HOMOLOGENE/homologene_human_all_minus_pseudo.txt'
infname_2 = root_dir + '/processed/HOMOLOGENE/homologene_mouse_all_minus_pseudo.txt'
infname_3 = root_dir + '/processed/HOMOLOGENE/homologene_human_protein_coding.txt'
infname_4 = root_dir + '/processed/HOMOLOGENE/homologene_mouse_protein_coding.txt'

# define log file
log_fname = root_dir + '/logs/homologene_merge.log'

# define outfiles
out_fname1 = root_dir + '/processed/HOMOLOGENE/hm_one_to_one_homol_ALL.txt'
out_fname2 = root_dir + '/processed/HOMOLOGENE/hm_one_to_many_homol_ALL.txt'
out_fname3 = root_dir + '/processed/HOMOLOGENE/hm_many_to_one_homol_ALL.txt'
out_fname4 = root_dir + '/processed/HOMOLOGENE/hm_many_to_many_homol_ALL.txt'

out_fname5 = root_dir + '/processed/HOMOLOGENE/hm_one_to_one_homol_PC.txt'
out_fname6 = root_dir + '/processed/HOMOLOGENE/hm_one_to_many_homol_PC.txt'
out_fname7 = root_dir + '/processed/HOMOLOGENE/hm_many_to_one_homol_PC.txt'
out_fname8 = root_dir + '/processed/HOMOLOGENE/hm_many_to_many_homol_PC.txt'


#open output files
outfile1 = open (out_fname1, 'w')
outfile2 = open (out_fname2, 'w')
outfile3 = open (out_fname3, 'w')
outfile4 = open (out_fname4, 'w')

outfile5 = open (out_fname5, 'w')
outfile6 = open (out_fname6, 'w')
outfile7 = open (out_fname7, 'w')
outfile8 = open (out_fname8, 'w')

# open log file for writing
log_file = open(log_fname, 'w')

# make lists

human_headers = ['HID','hGENE_ID','hGENE_SYMBOL']
mouse_headers = ['HID','mGENE_ID','mGENE_SYMBOL']

########################################################################
#functions

# process the homologene gene set files for human and mouse into a dataframe
def process_file(file, headers):
    try:
        in_file = open(file, 'r')
    except IOError, message:
        print >> sys.stderr, "infile could not be opened:"
        sys.exit ( 1 )
    print "opened: ", os.path.basename(file)
    records = pd.read_table(in_file , sep= '\t')
    records.columns = headers
    print '\n'
    return records

def get_relationships(human_records, mouse_records, out1, out2, out3, out4):

    #make empty lists
    one_to_one =[]
    one_to_many =[]
    many_to_one =[]
    many_to_many =[]

    # merge the two data frames via their homologene IDs
    human_mouse_merge = pd.merge(human_records,mouse_records, on='HID')
    #print human_mouse_merge

    # extract the unique homologene IDS from the merged list = list_ab
    #list_a=  human
    #list_b=  mouse
    #list_ab = uniqueID list for mouse and human combined

    # convert the human HID column to a list
    list_a = human_records['HID'].tolist()
    #print list_a

    #convert the mouse HID column to a list
    list_b = mouse_records['HID'].tolist()
    #print list_b

    list_ab = list(set (list_a + list_b))
    #print list_ab

    # extract the relationships
    for y in list_ab:
        #get the one to one relationships a to b.
         if list_a.count(y)==1 and list_b.count(y)==1:
           one_to_one.append(y)
           human_mouse_merge.loc[human_mouse_merge['HID'].isin([y])].to_csv(out1, sep='\t', header=False, index = False)

        #get the one to many relationships a to b.
         elif list_a.count(y)==1 and list_b.count(y)>1:
           one_to_many.append(y)
           human_mouse_merge.loc[human_mouse_merge['HID'].isin([y])].to_csv(out2, sep='\t', header=False, index = False)

        # get the many to one relationships  a to b.
         elif list_a.count(y)>1 and list_b.count(y)==1:
           many_to_one.append(y)
           human_mouse_merge.loc[human_mouse_merge['HID'].isin([y])].to_csv(out3, sep='\t', header=False, index = False)

        # get the many to many relationships  a to b.
         elif list_a.count(y)>1 and list_b.count(y)>1:
           many_to_many.append(y)
           human_mouse_merge.loc[human_mouse_merge['HID'].isin([y])].to_csv(out4, sep='\t', header=False ,index = False)


    # print the information to the log file.

    log_file.write("no of IDs in human id list  =  " + str(len(list_a)) + '\n')
    log_file.write("no of IDs in mouse id list =  " + str(len(list_b)) + '\n')
    log_file.write("no of IDs common to both lists =  " + str(len(list_ab)) + '\n')
    log_file.write("no of IDs in one to one list =  " + str(len(one_to_one)) + '\n')
    log_file.write("no of IDs in one to many list =  " + str(len(one_to_many)) + '\n')
    log_file.write("no of IDs in many to one list =  " + str(len(many_to_one)) + '\n')
    log_file.write("no of IDs in many to many list =  " + str(len(many_to_many)) + '\n')

########################################################################
# Main process
print '\n'
print 'process began at :', time.strftime("%Y-%m-%d %H:%M")

log_file.write('process began at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

# process the homologene gene sets for mouse and human
print '\n'
print 'processing human and mouse homologene gene sets'

human_all_records = process_file(infname_1,human_headers)
#print human_all_records
print 'processed human_all_records'
log_file.write('processed human_all_records' + '\n')

mouse_all_records = process_file(infname_2,mouse_headers)
#print mouse_all_records
print 'processed mouse_all_records'
log_file.write('processed mouse_all_records' + '\n')

human_pc_records = process_file(infname_3,human_headers)
#print human_pc_records
print 'processed human_pc_records'
log_file.write('processed human_pc_records' + '\n')

mouse_pc_records = process_file(infname_4,mouse_headers)
#print mouse_pc_records
print 'processed mouse_pc_records'
log_file.write('processed mouse_pc_records' + '\n')

#extract the one-to-one, one_to_many, many_to_one, many_to_many relationships for mouse and human ALL genes

print 'extracting the one-to-one, one_to_many, many_to_one, many_to_many relationships for ALL genes '

log_file.write('extracted data for human to mouse all genes minus pseudo genes relationships' + '\n')

get_relationships(human_all_records, mouse_all_records, outfile1, outfile2, outfile3, outfile4)

print 'extracting the one-to-one, one_to_many, many_to_one, many_to_many relationships for PROTEIN CODING genes '

log_file.write('extracted data for human to mouse PROTEIN CODING gene relationships' + '\n')

get_relationships(human_pc_records, mouse_pc_records, outfile5, outfile6, outfile7, outfile8)

# #########################################################################

log_file.write('process complete at:   ' + time.strftime("%Y-%m-%d %H:%M") + '\n')

print 'end of processing at : ', time.strftime("%Y-%m-%d %H:%M")

log_file.close()

########################################################################
