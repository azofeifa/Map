__author__ = 'Jonathan Rubin'

import os

def run(scriptdir, newpath, tempdir):
    outfile = open(scriptdir + '/runsamtobaidir.sh', 'w')
    outfile.write("id=" + newpath + "\n")
    outfile.write("od=" + newpath + "sortedbam/\n")
    outfile.write("mkdir -p $od\n")
    outfile.write("for pathandfilename in `ls $id*.sam`; do\n")
    outfile.write("entry=`basename $pathandfilename .sam`\n")
    outfile.write("echo $entry\n")
    outfile.write("qsub -v outdir=$od,indir=$id,basename=$entry -N ${entry}samtobai " + scriptdir + "/samtobai.sh\n")
    outfile.write("done")
    outfile.close()
    
    os.system("bash " + scriptdir + "/runsamtobaidir.sh > " + tempdir + "/Job_ID.txt")
    
    return newpath + 'sortedbam/'