import sys, re

'''
This program will accept a tabular BLAST output file (XXX.outfmt6) and 
a Trinity differential expression FPKM matrix (diffExp.XXX.matrix) file as inputs, 
and perform a transcript-to-protein lookup. 

inputs: XXX.outfmt6, diffExp.XXX.matrix
output: a file that contains transcript information and its corresponding protein 

To run this program, type:
"python3 parse_blast_match_DEprofile.py XXX.outfmt6 diffExp.XXX.matrix",
where the last two arguments are the file names.
'''

class BLAST_record():
    def __init__(self, blastInput):
        all_fields = blastInput.rstrip('\n').split('\t')
        
        self.transcript_ID = all_fields[0].split('|')[0]
        self.swissProt_ID  = all_fields[1].split('|')[3].split('.')[0]
        self.identity      = all_fields[2]
      

class DE_MATRIX():
    def __init__(self, DE_FILE_PATH):
        self.line = [x.split('\t') for x in open(DE_FILE_PATH).readlines()]


# function that will accept a BLAST object and return whether its identity attribute is > 95.
def analyzeBLAST(pident):
    if float(pident) > 95:
        return(True)
    else:
        return(False)
    

# function that accept a tuple and return it as a tab-separated string.
def tuple_to_string(tupleInput):
    return('\t'.join(tupleInput))


# Modify the program to accept command-line arguments and have it accept 2 file inputs: the blast file and the matrix file. 
if __name__== "__main__":
    
    # Ensure user inouts are correct first 
    try:
        blast_de_files = sys.argv
        blastfile, defile = sys.argv[1], sys.argv[2]
        assert blastfile.lower().endswith('.outfmt6')       
        assert re.match('(.*)diffExp.(.*).matrix', defile)
    
    except IndexError: 
        print('Incorrect number of inputs: 3 expected, but only ' + str(len(sys.argv)) + ' provided!\n' +
              'Note: This program accepts a tabular BLAST output file (XXX.outfmt6) and a \n \
     Trinity differential expression FPKM matrix (diffExp.XXX.matrix) file as inputs.')
    
    except AssertionError:
        print('Incorrect file type(s)!\n' +
              'Note: This program accepts a tabular BLAST output file (XXX.outfmt6) and a \n \
     Trinity differential expression FPKM matrix (diffExp.XXX.matrix) file as inputs.') 
    
    else:    
        # Loads the BLAST objects into a dictionary.
        with open(blastfile, "r") as blastf:
            blastfile_all_lines = blastf.readlines()
                       
            filteredDict = {blastob.transcript_ID:blastob.swissProt_ID 
                            for blastob in (BLAST_record(line) for line in blastfile_all_lines)
                            if analyzeBLAST(blastob.identity)}             
         
        # Performs a transcript-to-protein lookup.  
        DE = DE_MATRIX(defile)
        with open("Outputs/matched.txt", "w") as out, open("Outputs/notmatched.txt", "w") as e:
            for x in DE.line:
              
                if any(val == '' for val in x):
                    e.write("Skipping line: " + tuple_to_string(x).rstrip('\n') + " with error: line is missing fields\n")
                      
                else:
                    try:
                        out.write(filteredDict[x[0]] + "\t" + tuple_to_string(x[1:]))
                    except KeyError:
                        e.write("Skipping line: " + tuple_to_string(x).rstrip('\n') + "\t with error: no match in BLAST file\n")
