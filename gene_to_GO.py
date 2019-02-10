'''
Created on Nov 3, 2018

@author: zhentrinh
'''
import re

'''
This program will parse a gene ontology file to map a GO term to its parent terms
'''


# Define a function for splitting the file into records
def splitfile(file_path):
    with open(file_path) as f:
        data = f.read()        
        data_1 = re.split("\[Term\]", data)        
        return data_1
    

# Define a function to parse and capture the GO ID and its parent terms from a single GO term
def parseGO(term):
    ID_pattern = re.compile(r"^id:\s+(GO:[0-9]+)\n", re.M)
    is_a_pattern = re.compile(r"^is_a:\s+(GO:[0-9]+)", re.M) 
    ID = re.findall(ID_pattern, term)
    is_a = re.findall(is_a_pattern, term)
    return ID, is_a


# Define a function to store GO id and is_as returned from parseGO into a dictionary
def to_dict(allrecords):
    GO_dict = {}
    for record in allrecords:
        (go_id, is_a) = parseGO(record)
        key_inDict = "".join(go_id)
        GO_dict[key_inDict] = is_a
    return(GO_dict)
        

# Define a function to recursively looks for all of the parent terms of a given GO term
def findParent(child_ID, parent_list):   
    global all_dicts
    parents = all_dicts[child_ID]   
    if parents == []:
        return
    else:
        for parent in parents:
            if parent not in parent_list:
                parent_list.append(parent)
            findParent(parent, parent_list) 


# Define a function to map a protein ID to all its associated GO terms
def mapProteinToGO(annotationFile):
    with open(annotationFile) as f:
        annotation_set = set()
        for line in f:
            allfields = line.split("\t") 
            protein_ID = allfields[1]
            GO_ID = allfields[4]
            annotation_set.add((protein_ID, GO_ID))
        return(sorted(annotation_set)) 
         

# Define a function to write result to a file          
def print_report(outfile):
    with open(outfile, "w") as result:
        current_protein = None
        for annotation in all_annotation: 
            (protein, GO) = annotation
            if protein != current_protein:
                current_protein = protein
                result.write("%s" % (protein))
            if GO in all_dicts:                      
                parent_list = []
                findParent(GO, parent_list) 
#                 if protein != "" and GO != "" and parent_list != []:  
                result.write("\t%s\n" % (GO))
                for parent in parent_list:           
                    result.write("\t\t%s\n" % (parent))
      

if __name__ == "__main__":
    test_records = splitfile("scratch/go-basic.obo")    
    all_dicts = to_dict(test_records)
    all_annotation = mapProteinToGO("scratch/gene_association_subset.txt")
    print_report("Outputs/gene_to_go.txt")
    