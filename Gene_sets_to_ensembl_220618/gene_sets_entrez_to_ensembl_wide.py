
#!/usr/bin/python

# Written by Janet Harwood June 21st 2018.

import os
import sys
import re

###############################################################################

# get the curent working directory
cwd = os.getcwd()

# get the root directory
#root_dir = os.path.split(cwd)[0]

ensembl_all_fname = cwd + '/ENSG.genes.txt'

gene_set_path = cwd + '/gene_sets'

ncbi_all_fname = cwd + '/ncbi_human_entrez_to_ensembl_all_18_6_2018.txt'

# make outfile paths

out_dir_path = cwd + '/gene_sets_ensembl'

outfile_type = 'ensembl_magma.txt'

###############################################################################

def print_dictionary(d):
    for k, v in d.items():
        print k, "=>", v

##########################################
# Read in the PGC3FXMHC GRCh37 annotation file and make two dictionaries: protein coding gene dictionary and  'all'genes dictionary.

# headers: in PGC3FXMHC GRCh37 annotation file = ensembl_gene_id	external_gene_name	chromosome_name	start_position	end_position	strand	gene_biotype	hgnc_symbol	mim	entrezID	uniprot

ensembl_all_dict = {}
ensembl_pc_dict = {}

for line in open(ensembl_all_fname,'r'):
    if not line.startswith('ensembl_gene_id'):
        data = line.strip().split('\t')
        #print data
        ensembl_all_dict[data[9]]=[data[0],data[6],data[7],'ensembl_GRCh37']
        if data[6] == 'protein_coding':
            ensembl_pc_dict[data[9]]=[data[0]]

# print_dictionary(ensembl_all_dict)
# print_dictionary(ensembl_pc_dict)

##########################################
# Read in the ncbi gene annotations.

ncbi_all_dict = {}

for line in open(ncbi_all_fname,'r'):
        data = line.strip().split('\t')
        ncbi_all_dict[data[0]]=[data[2],data[3],data[1],'ncbi']

# print_dictionary(ncbi_all_dict)

##########################################

# Read in gene set file and make a default dictionary

for infilename in os.listdir(gene_set_path):
    if infilename.endswith('.txt'):
        gene_set_dict = {}
        gene_set_dict_edit ={}
        unmapped =set()
        not_in_ensembl =[]
        not_in_ncbi =[]
        subset = []
        count = 0
        print 'processing infile =', infilename
        in_fname = os.path.join(gene_set_path, infilename)
        #print '', in_fname
        file_id = re.split('\.', infilename)[0]
        #print 'file_id =',file_id
        #make outfile names from the infile name.
        outfilename1 = "_".join([file_id,outfile_type])
        print 'outfile1 =', outfilename1

        outfilename2 = "_".join([file_id,'entrez_ensembl_unmapped.txt'])
        print 'outfile2 =', outfilename2

        out_fname1 = os.path.join(out_dir_path, outfilename1)
        out_fname2 = os.path.join(out_dir_path, outfilename2)

        # open outfile
        outfile1 = open(out_fname1, 'w')
        outfile2 = open(out_fname2, 'w')

        for line in open(in_fname, 'r'):
            data = line.strip().split('\t')
            #print data
            genes = data[1].split(' ')
            #print genes
            gene_set_dict[data[0]] = genes

        # print_dictionary(gene_set_dict)
        print 'number of gene sets in infile = ',len(gene_set_dict)

#     # remove ids from the gene sets that don't map to Kyoko's annotation file and make an edited dictionary
#     # add the unmapped ids to a set

        for k,v in gene_set_dict.items():
            ensembl_genes = []
            for id in v:
                ensembl = ensembl_pc_dict.get(id)
                if ensembl is None:
                    unmapped.add(id)
                    i = v.index(id)
                    del v[i]
            gene_set_dict_edit[k] = v

        print 'number of gene sets after editing unmapped genes = ',len(gene_set_dict_edit)
        # print out the edited gene sets formatted with ensembl ids.

        for k,v in gene_set_dict_edit.items():
            ensembl_genes = []
            for id in v:
                ensembl = ensembl_pc_dict.get(id)
                # print ensembl
                if ensembl is not None:
                    ensembl_str = " ".join(ensembl)
                    ensembl_genes.append(ensembl_str)

            # output ensembl gene sets to file
            outstring = " ".join(ensembl_genes)
            outfile1.write(k + '\t' + outstring + '\n' )

      # check the unmapped ids in the all genes dictionary from Kyoko's file and in ncbi gene annotations.

        for id in unmapped:
            unmapped_info_ens = ensembl_all_dict.get(id)
            unmapped_info_ncbi = ncbi_all_dict.get(id)
            if unmapped_info_ens is not None:
                unmapped_info_ens_str = " ".join(unmapped_info_ens )
                outfile2.write(id + '\t' + unmapped_info_ens_str + '\n')
                #print '', id, unmapped_info_ens
            else:
                not_in_ensembl.append(id)

            if unmapped_info_ncbi is not None:
                unmapped_info_ncbi_str = " ".join(unmapped_info_ncbi)
                outfile2.write(id + '\t' + unmapped_info_ncbi_str + '\n')
                # print '', id, unmapped_info_ncbi
            else:
                not_in_ncbi.append(id)


        print 'mapped ensembl ids for ',infilename, 'written out to file'
        print 'ids that do not map to ensembl ids and were removed from the gene sets for ',infilename, '=' ,unmapped
        print 'annotations for unmapped id written to file', outfilename2
        print '\n'

print 'processing complete'









###########################################################################################################
