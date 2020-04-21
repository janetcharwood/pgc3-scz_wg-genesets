
#!/usr/bin/python

import os
import sys
import socket
import ftplib # for importing data from an ftp site.
import urllib # for importing data from a web site.

# get the curent working directory
cwd = os.getcwd()
# print 'current working directory is', cwd

# get the root directory
root = os.path.split(cwd)[0]
# print 'root directory is', root

dir_path = root + '/downloads'

# change directory to the path above
os.chdir(dir_path)

print 'Starting downloads from MGI'

testfile = urllib.URLopener()
testfile.retrieve("http://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt", "MGI_PhenoGenoMP.rpt")

print 'MGI_PhenoGenoMP.rpt downloaded'

testfile.retrieve("http://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt", "MGI_EntrezGene.rpt")

print 'MGI_EntrezGene.rpt downloaded'

testfile.retrieve("http://www.informatics.jax.org/downloads/reports/MPheno_OBO.ontology", "MPheno_OBO.ontology")

print 'MPheno_OBO.ontology downloaded'
print '\n'

####################################################################################################
#download data.

print 'starting downloads'

HOST_NCBI = 'ftp.ncbi.nlm.nih.gov' # path to website where NCBI data is hosted
HOST_GO = 'ftp.geneontology.org' # path to website where GO data is hosted


DIRN_NCBI_1 = 'gene/DATA/'   # path to directory where data is hosted
DIRN_NCBI_2 = 'pub/HomoloGene/current'   # path to directory where homologene data is hosted
DIRN_GO = 'go/ontology'   # path to directory where GO data is hosted

# files to download.

gene_info = 'gene_info.gz'
gene2go = 'gene2go.gz'
homologene_info = 'homologene.data'
go_obo = 'go-basic.obo'

# pass the host and the file and the directory to the function below.

def download(HOST,DIRN,FILE):
    try:
        f= ftplib. FTP(HOST)
    except(socket.error, socket.gaierror), e:
        print 'ERROR:cannot reach "%s"' % HOST
        return
    print '*** Connected to host "%s"'  % HOST

    try:
        f.login()
    except ftplib.error_perm:
        print 'ERROR:cannot login anonymously'
        f.quit()
        return
    print '*** Logged in as "anonymous"'

    try:
        f.cwd(DIRN)
    except ftplib.error_perm:
        print 'ERROR:cannot CD to "%s"' %DIRN
        f.quit()
        return
    print '*** Changed to "%s" folder' %DIRN


    try:
        f.retrbinary('RETR %s' %FILE,
         open(FILE,'wb').write)
    except ftplib.error_perm:
        print 'ERROR:cannot read file "%s"' %FILE
        os.unlink(FILE)
    else: print '*** Downloaded "%s" to CWD ' %FILE
    f.quit()
    return


# download NCBI files
download(HOST_NCBI,DIRN_NCBI_1,gene_info)
print 'NCBI gene_info.gz download complete'

download(HOST_NCBI, DIRN_NCBI_1, gene2go)
print 'NCBI gene2go download complete'

download(HOST_NCBI, DIRN_NCBI_2, homologene_info)
print 'NCBI homologene_info download complete'

download(HOST_GO, DIRN_GO, go_obo)
print 'GO go_obo download complete'

print 'all files downloaded, process complete'



###########################################################################################################
