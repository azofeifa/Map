__author__ = 'Jonathan Rubin'

import os

def run(trimdir,newpath):
    output = newpath + 'trimmed/'
    for file1 in os.listdir(newpath):
        if 'fastq' in file1.split('.')[-1]:
            os.system(trimdir + "trim_galore - o " + output + "/" + file1 + ".trimmed.fastq " + newpath + file1)

    return output