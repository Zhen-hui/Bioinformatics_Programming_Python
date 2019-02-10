'''
Created on Oct 30, 2018

@author: cathytrinh
'''

'''
This program uses regex to parse a gene ontology file "go-basic.obo"
'''

import re
import os


# define a class that will contain the necessary attribute values of a single GO term record.
class GO_attributes():
    
    def __init__(self, term):
        ID_pattern        = re.compile(r"^id:\s+(GO:[0-9]+)", re.M)
        name_pattern      = re.compile(r"^name:\s+(.*)\s+", re.M)
        namespace_pattern = re.compile(r"^namespace:\s+(.*?)\s+", re.M)
        is_a_pattern      = re.compile(r"^is_a:\s+(GO:[0-9]+.*?\n)", re.M)
         
        self.ID        = re.findall(ID_pattern, term)  
        self.name      = re.findall(name_pattern, term)
        self.namespace = re.findall(namespace_pattern, term)
        self.is_as      = re.findall(is_a_pattern, term)


    def writeResult(self, outfile = "Outputs/parsed_go_terms.txt"):    
        with open (outfile, "a") as p:
            p.write(self.ID[0] + "\t")
            p.write(self.namespace[0] + "\n")
            p.write("\t" + self.name[0] + "\n")
            for is_a in self.is_as:
                p.write("\t" + is_a)
            p.write("\n")
              

# define a function for splitting the file into records.
def splitRecords(file_path):
    with open(file_path) as f:
        allData = f.read()        
        go_terms = re.split("\[Term\]", allData)       

        # add each valid GO term object to a dictionary.
        GO_dict = {}
        for go_term in go_terms:   
            info = GO_attributes(go_term)
            if info.ID:
                GO_dict.update({info.ID[0]: info})
                
        # iterate through the dictionary and invoke the object's printing method
        for k in sorted (GO_dict.keys()):
            GO_dict[k].writeResult()        


# Export results
if os.path.exists("Outputs/parsed_go_terms.txt"):
    os.remove("Outputs/parsed_go_terms.txt")
    splitRecords("scratch/go-basic.obo")
else:          
    splitRecords("scratch/go-basic.obo")