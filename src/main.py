__author__ = 'Joey Azofeifa'
#Scripts adapted from Josephina Hendrix and Mary Allen

#Required modules:
#module load cutadapt_1.2.1
#module load sra_2.3.2-5
#module load samtools_0.1.19
#module load bedtools2_2.22.0
#Also need to git clone TrimGalore if using trim module into Map/ directory:
#git clone https://github.com/FelixKrueger/TrimGalore.git

import sys
import os
import trim_galore
import trimmomatic
import write_scripts
import sra_to_fastq
import check_job
import quality_check
import flip_reads
import fastq_to_sam
import sam_to_bam
import bam_to_bedgraph
import readcount_correction
import igv_create
import millions_mapped


#User-defined input
#======================================================================
#Give full path to SRA files or FASTQ files (remember to include a '/' at the end)
fullpath = sys.argv[1]

#Give email address where job e-mails will be sent
email="joseph.azofeifa@colorado.edu"

#Specify genome
genome = 'hg19'
#genome = 'dm3'
# genome = 'mm10'
# genome = 'ERCC'

#Specify bowtie options
#Used for ChIP-Seq
# bowtieoptions = "-k 1 -n 2 -l 36 --best"
#Used for GRO-Seq
bowtieoptions = "--very-sensitive -p 32"
# bowtieoptions = "-p 4 -S"
#Used for ATAC-Seq
#bowtieoptions = "-X2000"
#Used for RNAPII-ChIP
# bowtieoptions = "-n 1 -m 1-best-strata"

#Trim adaptors?
trimgalore = True
trimmomaticbool = False
#If no options desired use "" else needs a space at the end.
trimgaloreoptions = ""
# trimgaloreoptions = "--clip_R1 15 "

#Check read quality?
quality=False

#Flip reads? This is used for some GRO-Seq protocols
flip = False

#Check for Spike-In controls? Only True if you added spike-in controls from Jonathan Rubin to your GRO-Seq samples
spike= True
#======================================================================

#Return parent directory
def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

#Full path to genome construct and bowtie indexes. Genome files contain two tab separated columns: Chromosome, Length of chromosome.
#Bowtie indexes can be created with bowtie and a fasta file of your genome using the command:
#bowtie2-build genomefasta.fa basename
#genomefasta.fa = fasta file of entire genome
#basename = base filename given to bowtie index files
if genome == 'hg19':
    genomedir='/projects/dowellLab/groseq/forJoey/human.hg19.genome'
    bowtieindex='/projects/Down/Dowellseq/genomes/bowtiebwaindexs/hg19_Bowtie2_indexp32'
elif genome == 'dm3':
    genomedir='/projects/dowellLab/groseq/forJoey/dro/dm3.fa.genome'
    bowtieindex='/projects/Down/Dowellseq/genomes/bowtiebwaindexs/dm3.fa.Bowtie2'
elif genome == 'mm10':
    genomedir='/projects/dowellLab/groseq/forJoey/mm10.genome'
    bowtieindex='/projects/Down/Dowellseq/genomes/bowtiebwaindexs/mm10_Bowtie2_index'
elif genome == 'ERCC':
    bowtieindex='/projects/Down/Dowellseq/genomes/bowtiebwaindexs/ERCC92_Bowtie2_index'
    genomedir='/projects/Down/Dowellseq/genomes/bowtiebwaindexs/ERCC.genome'
else:
    sys.exit("Genome not found, exiting...")

#Path to spike-in control genomes
SpikeIngenomes=['/projects/Down/Dowellseq/genomes/LBS-1.genome','/projects/Down/Dowellseq/genomes/A-region.genome','/projects/Down/Dowellseq/genomes/LBS-3.genome','/projects/Down/Dowellseq/genomes/C-unit.genome','/projects/Down/Dowellseq/genomes/TRNAS23.genome']

#Path to spike-in control bowtie indexes
SpikeInbowtieindexes=['/projects/Down/Dowellseq/genomes/bowtiebwaindexs/LBS-1','/projects/Down/Dowellseq/genomes/bowtiebwaindexs/A-region','/projects/Down/Dowellseq/genomes/bowtiebwaindexs/LBS-3','/projects/Down/Dowellseq/genomes/bowtiebwaindexs/C-unit','/projects/Down/Dowellseq/genomes/bowtiebwaindexs/TRNAS23']


#Home directory
homedir = os.path.dirname(os.path.realpath(__file__))

#Scripts directory
scriptdir = parent_dir(homedir) + '/scripts'

#Temporary files directory
tempdir = parent_dir(homedir) + '/temp'

#Directory to temporary job file
job = tempdir + "/Job_ID.txt"

#Trim Galore directory
trimdir = parent_dir(homedir) + '/TrimGalore/'

#Trimmomatic directory
trimmomaticdir = parent_dir(homedir) + '/Trimmomatic-0.36/trimmomatic-0.36.jar'


def run():
    print "Filepath: ", fullpath
    newpath = fullpath
    
    #Converts SRA to Fastq format
    print "Converting SRA to Fastq..."
    boolean = sra_to_fastq.run(scriptdir, newpath, tempdir)
    if boolean:
        check_job.run(job,tempdir)
        print "done"
    else:
        print "No SRA files in filepath"

    #Trim adaptor sequences from fastq files
    if trimgalore:
        newpath = trim_galore.run(trimdir,trimgaloreoptions,newpath)
    elif trimmomaticbool:
        newpath = trimmomatic.run(trimmomaticdir,newpath)

    
    #Checks read quality
    if quality:
        print "Checking quality..."
        quality_check.run(scriptdir, newpath, tempdir)
        check_job.run(job,tempdir)
    
    #Flips reads (use for some GRO-Seq protocols)
    if flip:
        print "done\nFlipping Reads..."
        newpath = flip_reads.run(scriptdir, newpath, tempdir)
        check_job.run(job,tempdir)
    
    #Checks reads mapped to spike-in control genomes
    if spike:
        print "done\nChecking Spike-in Controls..."
        if not os.path.exists(newpath + 'bowtie2/'):
            os.mkdir(newpath + 'bowtie2/')
        if not os.path.exists(newpath + 'bowtie2/spikeins/'):
            os.mkdir(newpath + 'bowtie2/spikeins/')
        for i in range(len(SpikeIngenomes)):
            print "Checking reads mapped to: " + SpikeInbowtieindexes[i].split('/')[-1]
            g = SpikeIngenomes[i]
            b = SpikeInbowtieindexes[i]
            write_scripts.run(scriptdir,g,b,bowtieoptions,email)
            fastq_to_sam.run(scriptdir, newpath, tempdir, g, boolean=False)
            check_job.run(job,tempdir)
        for file1 in [i for i in os.listdir(newpath + 'bowtie2/spikeins/') if '.stderr' in i]:
            print file1.split('.')[0] + file1.split('.')[2]
            with open(newpath + 'bowtie2/spikeins/' + file1) as F:
                for line in F:
                    print line

    #Writes script files based on genome and bowtie index
    print "done\nWriting script files..."
    write_scripts.run(scriptdir,genomedir,bowtieindex,bowtieoptions,email)

    #Converts Fastq to SAM format
    print "done\nConverting Fastq to SAM..."
    newpath = fastq_to_sam.run(scriptdir, newpath, tempdir, genomedir)
    check_job.run(job,tempdir)
    
    #Converts SAM to BAM format
    print "done\nConverting SAM to BAM..."
    newpath = sam_to_bam.run(scriptdir, newpath, tempdir)
    check_job.run(job,tempdir)
    
    #Converts BAM to Bedgraph format
    print "done\nConverting BAM to Bedgraph..."
    bam_to_bedgraph.run(scriptdir, newpath, tempdir)
    check_job.run(job,tempdir)
    
    #Normalizes bedgraphs to millions mapped
    print "done\nCorrecting for readcounts..."
    readcount_correction.run(scriptdir, newpath)
    
    #Creates IGV scripts 
    print "done\nCreating IGV files..."
    igv_create.run(scriptdir, newpath, genome)
    
    
    print "done\nGetting millions mapped reads..."
    millions_mapped.run(scriptdir, fullpath)
    
    print "done"
    
    
    
    
