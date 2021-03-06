__author__ = 'Jonathan Rubin'

import os

def run(scriptdir, newpath, tempdir):
    outfile = open(scriptdir + '/runqual.sh', 'w')
    outfile.write("id=" + newpath + "\n")
    outfile.write("for pathandfilename in `ls $id*.fastq`; do\n")
    outfile.write("entry=`basename $pathandfilename`\n")
    outfile.write("qsub -v filename=$entry,indir=$id -N ${fn1}qual " + scriptdir + "/qual.sh\n")
    outfile.write("done")
    outfile.close()
    
    os.system("bash " + scriptdir + "/runqual.sh > " + tempdir + "/Job_ID.txt")
    
    return tempdir + "/Job_ID.txt"