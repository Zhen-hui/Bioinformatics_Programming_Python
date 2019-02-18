'''
Created on Nov 3, 2018

@author: zhentrinh
'''
import re
from typing import Tuple

'''
This program will parse a gene ontology file to map a GO term to its parent terms
'''


def split_file(file_path):
    """
    Splitting the file into records
    :param file_path: the path of your gene ontology file
    :return: data after splitting
    """
    with open(file_path) as f:
        data = f.read()        
        data_1 = re.split("\[Term\]", data)        
        return data_1


def parse_go(term):
    """
    # Parse and capture the GO ID and its parent terms from a single GO term
    :param term: each gene record
    :return: ID and is_a terms from each gene record
    """
    id_pattern = re.compile(r"^id:\s+(GO:[0-9]+)\n", re.M)
    is_a_pattern = re.compile(r"^is_a:\s+(GO:[0-9]+)", re.M) 
    ID = re.findall(id_pattern, term)
    is_a = re.findall(is_a_pattern, term)
    return ID, is_a


def to_dict(all_records):
    """
    # Store GO id and is_as returned from parse_go into a dictionary
    :param all_records: all records returned from parse_go function
    :return: a dictionary with GO term ID as key, and is_a as values
    """
    go_dict = {}
    for record in all_records:
        (go_id, is_a) = parse_go(record)
        key_in_dict = "".join(go_id)
        go_dict[key_in_dict] = is_a
    return go_dict
        

def find_parent(child_id, parent_list):
    """
    # Recursively looks for all of the parent terms of a given GO term
    :param child_id: a GO term ID
    :param parent_list: a list of all parents of a given GO term
    :return: a list of all parents of a given GO term
    """
    global all_dicts
    parents = all_dicts[child_id]
    if parents == []:
        return
    else:
        for parent in parents:
            if parent not in parent_list:
                parent_list.append(parent)
            find_parent(parent, parent_list)


def map_protein_to_go(annotation_file):
    """
    # Map a protein ID to all its associated GO terms
    :param annotation_file: a file path to write results to
    :return: a file containing all the mapped protein and GO term
    """
    with open(annotation_file) as f:
        annotation_set = set()
        for line in f:
            all_fields = line.split("\t")
            protein_id = all_fields[1]
            go_id = all_fields[4]
            annotation_set.add((protein_id, go_id))
        return sorted(annotation_set)
         

def print_report(outfile):
    """
    # Write result to a file
    :param outfile: a file path to write results to
    :return: an output file
    """
    with open(outfile, "w") as result:
        current_protein = None
        annotation: Tuple[str, str]
        for annotation in all_annotation: 
            (protein, GO) = annotation
            if protein != current_protein:
                current_protein = protein
                result.write("%s" % protein)
            if GO in all_dicts:                      
                parent_list = []
                find_parent(GO, parent_list)
                result.write("\t%s\n" % GO)
                for parent in parent_list:           
                    result.write("\t\t%s\n" % parent)
      

if __name__ == "__main__":
    test_records = split_file("scratch/go-basic.obo")
    all_dicts = to_dict(test_records)
    all_annotation = map_protein_to_go("scratch/gene_association_subset.txt")
    print_report("Outputs/gene_to_go.txt")
