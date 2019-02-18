"""
Created on Oct 30, 2018

@author: cathytrinh
"""
import re
import os

'''
This program uses regex to parse a gene ontology file "go-basic.obo"
'''


# define a class that will contain the necessary attribute values of a single GO term record.
class GOAttributes:
    
    def __init__(self, term):
        id_pattern = re.compile(r"^id:\s+(GO:[0-9]+)", re.M)
        name_pattern = re.compile(r"^name:\s+(.*)\s+", re.M)
        namespace_pattern = re.compile(r"^namespace:\s+(.*?)\s+", re.M)
        is_a_pattern = re.compile(r"^is_a:\s+(GO:[0-9]+.*?\n)", re.M)
         
        self.ID = re.findall(id_pattern, term)
        self.name = re.findall(name_pattern, term)
        self.namespace = re.findall(namespace_pattern, term)
        self.is_as = re.findall(is_a_pattern, term)

    def write_result(self, outfile="Outputs/parsed_go_terms.txt"):
        with open(outfile, "a") as p:
            p.write(self.ID[0] + "\t")
            p.write(self.namespace[0] + "\n")
            p.write("\t" + self.name[0] + "\n")
            for is_a in self.is_as:
                p.write("\t" + is_a)
            p.write("\n")
              

def split_records(file_path):
    """
    # define a function for splitting the file into records.
    :param file_path: file path for reading
    :return:
    """
    with open(file_path) as f:
        all_data = f.read()
        go_terms = re.split("\[Term\]", all_data)

        # add each valid GO term object to a dictionary.
        go_dict = {}
        for go_term in go_terms:   
            info = GOAttributes(go_term)
            if info.ID:
                go_dict.update({info.ID[0]: info})
                
        # iterate through the dictionary and invoke the object's printing method
        for k in sorted(go_dict.keys()):
            go_dict[k].write_result()


# Export results
if os.path.exists("Outputs/parsed_go_terms.txt"):
    os.remove("Outputs/parsed_go_terms.txt")
    split_records("scratch/go-basic.obo")
else:          
    split_records("scratch/go-basic.obo")