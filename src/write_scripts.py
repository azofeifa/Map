__author__ = 'Jonathan Rubin'

import os

def run(scriptdir, genomedir, bowtieindex, bowtieoptions, email):
    outfile = open(scriptdir + '/bam_to_5primebed.sh', 'w')
    outfile.write("### Run in desired queue\n")
    outfile.write("#PBS -q long8gb\n")
    outfile.write("### Use the bourne shell\n")
    outfile.write("#PBS -S /bin/sh\n")
    outfile.write("### Specify the number of nodes and processors for your job\n")
    outfile.write("#PBS -l nodes=1:ppn=1\n")
    outfile.write("#PBS -o /projects/dowellLab/groseq/pubgro/e_and_o/\n")
    outfile.write("#PBS -e /projects/dowellLab/groseq/pubgro/e_and_o/\n")
    outfile.write("### Set your email address\n")
    outfile.write("#PBS -m ae\n")
    outfile.write("#PBS -M " + email + "\n")
    outfile.write("### Switch to the working directory; by default TORQUE launches processes\n")
    outfile.write("### from your home directory.  This is a good idea because your -o and -e files\n")
    outfile.write("### will go here\n")
    outfile.write("cd $PBS_O_WORKDIR\n")
    outfile.write("echo Working directory is $PBS_O_WORKDIR\n")
    outfile.write("### Retrieve/use all modules loaded ###\n")
    outfile.write("#PBS -V\n")
    outfile.write("genome=" + genomedir + "\n")
    outfile.write("echo $infile\n")
    outfile.write("echo $genome\n")
    outfile.write("echo $outfile1\n")
    outfile.write("mkdir -p  $outdir/genomecoveragebed\n")
    outfile.write("mkdir -p $outdir/forFstitch\n")
    outfile.write("mkdir -p $outdir/genomecoveragebed/fortdf\n")
    outfile.write("/opt/bedtools/2.22.0/genomeCoverageBed -5 -bg -strand + -ibam $infile -g $genome > $outdir/forFstitch/$outfile1\n")
    outfile.write("/opt/bedtools/2.22.0/genomeCoverageBed -5 -bg -strand - -ibam $infile -g $genome > $outdir/forFstitch/$outfile2\n")
    outfile.write("/opt/bedtools/2.22.0/genomeCoverageBed -bg -strand + -ibam $infile -g $genome > $outdir/genomecoveragebed/$outfile3\n")
    outfile.write("/opt/bedtools/2.22.0/genomeCoverageBed -bg -strand - -ibam $infile -g $genome | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }'> $outdir/genomecoveragebed/$outfile4\n")
    outfile.write("cat $outdir/genomecoveragebed/$outfile4 $outdir/genomecoveragebed/$outfile3 > $outdir/genomecoveragebed/fortdf/$outfile5.bed\n")
    outfile.write("/opt/bedtools/2.22.0/sortBed -i $outdir/genomecoveragebed/fortdf/$outfile5.bed >$outdir/genomecoveragebed/fortdf/$outfile5\n")
    outfile.close()
    
    outfile2 = open(scriptdir + '/bowtieafastq.sh','w')
    outfile2.write("#PBS -q long2gb\n")
    outfile2.write("### Use the bourne shell\n")
    outfile2.write("#PBS -S /bin/sh\n")
    outfile2.write("### Specify the number of nodes and processors for your job\n")
    outfile2.write("#PBS -l nodes=1:ppn=32\n")
    outfile2.write("### Set your email address\n")
    outfile2.write("#PBS -m ae\n")
    outfile2.write("#PBS -M " + email + "\n")
    outfile2.write("### Switch to the working directory; by default TORQUE launches processes\n")
    outfile2.write("### from your home directory.  This is a good idea because your -o and -e files\n")
    outfile2.write("### will go here\n")
    outfile2.write("cd $PBS_O_WORKDIR\n")
    outfile2.write("echo Working directory is $PBS_O_WORKDIR\n")
    outfile2.write("### Retrieve/use all modules loaded ###\n")
    outfile2.write("#PBS -V\n")
    outfile2.write("echo $fastq1pathandfile\n")
    outfile2.write("echo ${outdir}${outfile}.sam\n")
    # outfile2.write("/opt/bowtie/bowtie2-2.0.2/bowtie2 -p32 " + bowtieoptions + " -un ${outdir}${outfile}.unmapped.fastq \\\n")
    outfile2.write("/opt/bowtie/bowtie2-2.0.2/bowtie2 -p32 " + bowtieoptions + " \\\n")
    outfile2.write(bowtieindex + " \\\n")
    outfile2.write("-U $fastq1pathandfile \\\n")
    outfile2.write("> ${outdir}${outfile}.sam \\\n")
    outfile2.write("2> ${outdir}${outfile}.stderr\n")
    outfile2.close()
