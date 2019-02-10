'''
Created on Oct 19, 2018

@author: zhentrinh
'''

'''
This program read in a BLASTp file, which maps DNA against DNA, for example 
gene sequences against a reference genome. 
Since there are multiple columns in the file, and we are only interested 
in getting a few columns, so this program will export the desired columns 
as a separate file.
'''

# Read in file 
blastfile = open("scratch/RNASeq/blastp.outfmt6")
all_lines = blastfile.readlines()

# Open file for writing
output = open("Outputs/parsed_blast.txt", "w")

# Grab the transcript ID, isoform, SwissProt ID, and percent of identical matches from each line
for line in all_lines:
    element = line.replace('|', '\t').split()
    transcriptID = element[0]
    isoform = element[1]
    SwissProtID = element[5].split(".")[0] 
    pident = element[7] 
    # Print the 4 fields on tab-separated format to a file
    output.write(transcriptID  + '\t' + isoform +  '\t' + SwissProtID +  '\t' + pident + '\n')

blastfile.close()
output.close()
