from __future__ import division
import glob
import os
import sys


def main2(directory_of_paper):
	save_path = directory_of_paper
	directory_of_sortedbams = os.path.join(directory_of_paper, "flipped/bowtie2/sortedbam")
	keyword = directory_of_paper.split("/")[-1]
	completeName = os.path.join(save_path, "millions_mapped.txt")
	file2 = open(completeName, "w")
	file2.truncate()
	#print directory_of_sortedbams
	for sorted_bam_file_and_path in glob.glob(os.path.join(directory_of_sortedbams, '*sorted.bam.flagstat')):
			bamfileroot = sorted_bam_file_and_path.split("/")[-1]
			n = 1
			while n < len(bamfileroot):
				character = bamfileroot[n]
				if str(character) == ".":
					break
				else:
					n = n+1
			cropped_bamfileroot = bamfileroot[0:n]
			f = open(sorted_bam_file_and_path)
			lines = f.readlines()
			total_reads = lines[0]
			total_reads = int(total_reads.split(" ")[0])
			mapped_reads = lines[2]
			mapped_reads = int(mapped_reads.split(" ")[0])
			percent_mapped = mapped_reads/total_reads
                        #print bamfileroot, total_reads, mapped_reads, percent_mapped   
			file2 = open(completeName, "a")
			value = keyword, cropped_bamfileroot, total_reads, mapped_reads, percent_mapped
			#cropped_bamfileroot, total_reads, mapped_reads, percent_mapped
			s = str(value)
			file2.write(s)
			file2.write('\n')
                        #file2.write(bamfileroot + ' ' + total_reads)
			#file2.write( total_reads)
			#file2.write( mapped_reads)
			#file2.write( percent_mapped+ "\n")
			file2.close()

if __name__=="__main__":
	main2(sys.argv[1])
	#f2.close()
