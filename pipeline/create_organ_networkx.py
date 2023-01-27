#the overall role of this file is to read the mesh ascii file, downloadable from tehir website
#and output the largest network possible - in includes all "main letters" A-N, V and Z
#at a later point in the snakemake procedure, the network will be truncated
#the network that we create is really a number of separated networks, where the number is equal to the number of 
#letter+number firt combinations, that is, A01, A02,... Z13,

#for a single "string" such as "Leukocyte"
#the mesh network can have paths diverge (that reconverge) for a single A01
#it can also have multiple paths that never touch between different A01, A15, etc
#we create our network such that paths in a single A01 never reconverge
#this will create conveniences when we fetch instances of organs later down the line


import os
import networkx as nx
from collections import defaultdict
from pprint import pprint
import matplotlib.pyplot as plt
import sys

#author Uli Koehler provides function readMeSh
#__author__ = "Uli Koehler"
#__copyright__ = "Copyright 2015, Uli Kö221hler"
#__license__ = "CC0 1.0 Universal"
#__version__ = "1.0"
def readMeSH(fin):
    """
    Given a file-like object, generates MeSH objects, i.e.
    dictionaries with a list of values for each qualifier.
    Example: {"MH": ["Acetylcysteine"]}
    """
    currentEntry = None
    for line in fin:
        line = line.strip()
        if not line:
            continue
        # Handle new record. MeSH explicitly marks this
        if line == "*NEWRECORD":
            # Yiel old entry, initialize new one
            if currentEntry:
                yield currentEntry
            currentEntry = defaultdict(list)
            continue
        # Line example: "MH = Acetylcysteine"
        key, _, value = line.partition(" = ")
        # Append to value list
        currentEntry[key].append(value)
    # If there is a non-empty entry left, yield it
    if currentEntry:
        yield currentEntry

def add_nodepath_and_label_to_endnode_to_networkx(temp_nx,temp_mesh_entry):
    '''
    We receive a networkx label and a single mesh entry
    we split the MN into multiple labels
    we split each label into a list of perpetually growing strings (A01, A01.032, A01.032,047)
    we add the "word label" at the end from teh MH
    '''

    #MN and MH are 'attributes' in the ascii text file
    #MN is all paths
    #MH is the end node
    nodepath_string_path_list=temp_mesh_entry['MN']
    
    #confirm that we are adding the right label always because there is only one
    if (len(temp_mesh_entry['MH']))>1:
        print(temp_mesh_entry['MH'])
        hold=input('found an entry with multiple labels')
    end_node_label=temp_mesh_entry['MH'][0]


    for temp_string_path in nodepath_string_path_list:
        node_path_elements=temp_string_path.split('.')
        node_paths=list()

        for i in range(0,len(node_path_elements)):
            node_paths.append('.'.join(node_path_elements[0:i+1]))

        #if 'A11' in node_paths:
            nx.add_path(temp_nx,node_paths)
        temp_nx.nodes[node_paths[-1]]['mesh_label']=end_node_label


if __name__ == "__main__":

    min_fold_change=sys.argv[1]
    #cores_available=sys.argv[2]
    mesh_file_address='../resources/mesh_ascii/d2021.bin'
    output_organ_networkx_address='../results/'+str(min_fold_change)+'/step_2a_create_organ_and_disease_networkx/mesh_organ_networkx.bin'
    output_disease_networkx_address='../results/'+str(min_fold_change)+'/step_2a_create_organ_and_disease_networkx/mesh_disease_networkx.bin'
    os.system('mkdir -p ../results/'+str(min_fold_change)+'/step_2a_create_organ_and_disease_networkx/')
    os.system('touch ../results/'+str(min_fold_change)+'/step_2a_create_organ_and_disease_networkx/dummy.txt')


    organ_networkx=nx.DiGraph()
    
    #Example of how to use readMeSH()
    #import argparse
    #parser = argparse.ArgumentParser()
    #parser.add_argument("file")
    #args = parser.parse_args()
    with open(mesh_file_address, "r") as infile:
        # readMeSH() yields MeSH objects, i.e. dictionaries
        for entry in readMeSH(infile):
             add_nodepath_and_label_to_endnode_to_networkx(organ_networkx,entry)

    #we add a parent node to all the 'C' trees because those are the disease trees, and we 
    #want to conveniently take those and make them their own networkx
    #this list was manually obtained by going to the mesh website and looking at all of the headnodes that they had
    for temp in ['C01','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23','C24','C25','C26']:
        nx.add_path(organ_networkx,['C',temp])
    organ_networkx.nodes['C']['mesh_label']='Disease'

    #split the total network into two networks, one for diseases and one for organs
    disease_networkx=organ_networkx.subgraph(nx.algorithms.dag.descendants(organ_networkx,'C')).copy()
    organ_networkx.remove_nodes_from(nx.algorithms.dag.descendants(organ_networkx,'C'))
    organ_networkx.remove_nodes_from('C')

    #prove that all 'C' are gone from organ hierarchy
    for temp in organ_networkx.nodes:
        if 'C' in temp:
            print('theres still a disease ehre')
    print(disease_networkx.nodes)
    #prove that the disease hierarchy is filled with nodes (albeit not showing edges)
    #hold=input('check diesases')

    nx.readwrite.gpickle.write_gpickle(organ_networkx,output_organ_networkx_address)
    nx.readwrite.gpickle.write_gpickle(disease_networkx,output_disease_networkx_address)