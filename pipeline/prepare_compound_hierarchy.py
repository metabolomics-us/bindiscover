import sys
import obonet
import matplotlib.pyplot as plt
import networkx as nx
import pandas
import numpy as np
import os
from pprint import pprint

def obtain_deepest_classyfire_class_per_bin(temp_panda):
    '''
    we want to obtain the deepest class possible for each bin to append each bin
    to the classyfire network

    so we scroll through each bin's class listing, starting with the deepest
    and if one is not null, then we set that value to the deepest class column
    '''
    
    temp_panda['deepest_class']='pre_analysis'

    class_list=[
        'direct_parent_5',
        'direct_parent_4',
        'direct_parent_3',
        'direct_parent_2',
        'direct_parent_1',
        'Subclass',
        'Class',
        'Superclass',
        'Kingdom'
    ]

    for index, series in temp_panda.iterrows():
        for temp_column in class_list:
            ##only nan is not equal to itself
            if series[temp_column] == series[temp_column]:
                temp_panda.at[index,'deepest_class']=series[temp_column]
                break

def make_class_to_node_name_dict(temp_nx):
    #get all node names
    #put in list
    #for each node name
    #declare the value to be that node's name and the key to be the name
    #because we get the classes (names) from classyfire
    node_name_list=list(temp_nx.nodes)
    class_to_node_dict=dict()
    for temp_node_name in node_name_list:
        class_to_node_dict[temp_nx.nodes[temp_node_name]['name']]=temp_node_name

    return class_to_node_dict

def add_one_node_to_classyfire_network(temp_nx,temp_bin,temp_class_to_node_dict):
    '''

    '''
    #get the class to use as the key for the class:node_name dict
    current_bin_name=temp_bin['deepest_class']

    #get the name of the node that this bin will connect to
    #plb 2-6-2022
    #basically, if the compound is unknown, then we connect directly to the root node
    ##update 7-4-2022
    ##too much data to do this analysis for unknowns
    try:
        node_to_connect_to=temp_class_to_node_dict[current_bin_name]
    except KeyError:
        return

    #add this bin as a node in the network
    #the id number is the name of the node
    temp_nx.add_node(
        int(temp_bin['id']),
        inchikey=temp_bin['inchikey'],
        type_of_node='from_binvestigate',
        common_name=temp_bin['name']
    )

    #add a connection between the added node and the class that it is most specifically
    #identified as
    temp_nx.add_edge(temp_bin['id'],node_to_connect_to,'is_a')

def add_all_bins_to_network(temp_nx,temp_panda,temp_class_to_node_dict):
    for index,series in temp_panda.iterrows():
        print(index)
        add_one_node_to_classyfire_network(temp_nx,series,temp_class_to_node_dict)

def visualize_added_classes(temp_nx,temp_original_classyfire_nodecount):
    color_list_original=['#1f78b4' for i in range(0,temp_original_classyfire_nodecount)]
    color_list_new=['#32cd32' for i in range(0,len(temp_nx.nodes)-temp_original_classyfire_nodecount)]
    total_color_list=color_list_original+color_list_new
    nx.draw(temp_nx,with_labels=True,node_color=total_color_list,node_size=50)
    plt.show()    
    
def remove_branches_without_fold_matrices(temp_nx):
    '''
    we want to remove unrelated branches so that we can do the recursive analysis

    there are two conditions for the keeping of a node, either

    the node itself has a numerical name (meaning that it came from a bin)
    or the ancestors of the bin have a number (meaning that it will be a node merge point later)
    '''

    nodes_to_remove=list()

    for i, temp_node in enumerate(temp_nx.nodes):

        if type(temp_node) != str:
            continue

        temp_ancestor_list=nx.algorithms.dag.ancestors(temp_nx,temp_node)
        has_numerical_value_in_ancestors=False

        for temp_ancestor in temp_ancestor_list:
            if type(temp_ancestor) != str:
                has_numerical_value_in_ancestors=True
                break
        
        if has_numerical_value_in_ancestors==False:
            nodes_to_remove.append(temp_node)

    temp_nx.remove_nodes_from(nodes_to_remove)

def do_everything(temp_nx,temp_node_keep_address,temp_hierarchy_type):
    
    if temp_hierarchy_type=='compound':
        
        temp_nx=nx.DiGraph.reverse(temp_nx)
        #draw_nx_for_analysis(temp_nx,'compound')
        compounds_to_keep_panda=pandas.read_csv(temp_node_keep_address)
        keep_list=compounds_to_keep_panda['nodes_to_keep'].to_list()    
        remove_unwanted_nodes(temp_nx,keep_list,temp_hierarchy_type)
        #draw_nx_for_analysis(temp_nx,'compound')

    return temp_nx

def remove_unwanted_nodes(temp_nx,temp_nodes_to_keep,temp_hierarchy_type):
    '''
    for species and compounds, we automatically keep a node if it is a leaf
    for organs and diseases, we require that a node be explicitly listed to be kept
    '''
    
    nodes_to_drop=list()
    if temp_hierarchy_type=='compound' or temp_hierarchy_type=='species':
        for temp_node in temp_nx.nodes:
            if (temp_node not in temp_nodes_to_keep) and (len(list(temp_nx.successors(temp_node)))>0):
                nodes_to_drop.append(temp_node)
    elif temp_hierarchy_type=='organ' or temp_hierarchy_type=='disease':
        for temp_node in temp_nx.nodes:
            if (temp_node not in temp_nodes_to_keep):
                nodes_to_drop.append(temp_node)


    for temp_node in nodes_to_drop:

        list_of_predecessors=list(temp_nx.predecessors(temp_node))
        list_of_successors=list(temp_nx.successors(temp_node))
        for temp_predecessor in list_of_predecessors:
            print(temp_predecessor)
            for temp_successor in list_of_successors:
                print(temp_successor)
                temp_nx.add_edge(temp_predecessor,temp_successor)
        temp_nx.remove_node(temp_node)

    return temp_nx

def draw_nx_for_analysis(temp_nx,temp_hierarchy_type,set_of_binvestigate_nodes=None):
    '''

    '''
    #for compound matrix, it is obvious that bins => source of data
    #for species matrix, we put that information into an attribute called "type of node" where the values
    #   are "from binvestigate" or "combination"
   

    total_color_list=list()
    
    if temp_hierarchy_type=='compound':
        label_dict=get_labels_for_drawing(temp_nx,'compound')
        for i,temp_node in enumerate(temp_nx.nodes):
            if temp_nx.nodes[temp_node]['type_of_node']=='combination':
                total_color_list.append('#ff0000')
            elif temp_nx.nodes[temp_node]['type_of_node']=='from_binvestigate':
                total_color_list.append('#32cd32')
    
    elif temp_hierarchy_type=='species':
        label_dict=get_labels_for_drawing(temp_nx,'species')
        for i,temp_node in enumerate(temp_nx.nodes):
            if temp_node in set_of_binvestigate_nodes:
                total_color_list.append('#32cd32')
            else:
                total_color_list.append('#ff0000')
    
    elif temp_hierarchy_type=='organ':
        label_dict=get_labels_for_drawing(temp_nx,'organ')
        for i,temp_node in enumerate(temp_nx.nodes):
            if temp_nx.nodes[temp_node]['mesh_label'] in set_of_binvestigate_nodes:
                total_color_list.append('#32cd32')
            else:
                total_color_list.append('#ff0000')

    elif temp_hierarchy_type=='disease':
        label_dict=get_labels_for_drawing(temp_nx,'disease')
        for i,temp_node in enumerate(temp_nx.nodes):
            if temp_nx.nodes[temp_node]['mesh_label'] in set_of_binvestigate_nodes:
                total_color_list.append('#32cd32')
            else:
                total_color_list.append('#ff0000')

    pos = nx.nx_agraph.pygraphviz_layout(temp_nx, prog='dot')
    nx.draw(temp_nx, pos,labels=label_dict,node_color=total_color_list)
    plt.show()

def get_labels_for_drawing(temp_nx,temp_hierarchy_type):
    label_dict=dict()
    if temp_hierarchy_type=='compound':
        for temp_node in temp_nx.nodes:
            try:
                label_dict[temp_node]=(temp_node+'\n'+temp_nx.nodes[temp_node]['name'])
            except TypeError:
                label_dict[temp_node]=(str(temp_node))

    elif temp_hierarchy_type=='species':
        for temp_node in temp_nx.nodes:
            pprint(temp_nx.nodes[temp_node])
            try:
                temp_sci_name=temp_nx.nodes[temp_node]['scientific_name']
            except KeyError:
                temp_sci_name='no scientific name found'

            try:
                if type(temp_nx.nodes[temp_node]['common_name']) == list:
                    temp_common_name=temp_nx.nodes[temp_node]['common_name'][0]
                else:
                    temp_common_name=temp_nx.nodes[temp_node]['common_name']
            except KeyError:
                temp_common_name='no common name found' 

            label_dict[temp_node]=(temp_node+'\n'+temp_sci_name+'\n'+temp_common_name)

    elif temp_hierarchy_type=='organ' or temp_hierarchy_type=='disease':
        for temp_node in temp_nx.nodes:
            print(temp_node)
            label_dict[temp_node]=temp_node+'\n'+temp_nx.nodes[temp_node]['mesh_label']

    return label_dict


if __name__ == "__main__":
    
    min_fold_change=sys.argv[1]
    os.system('mkdir -p ../results/'+str(min_fold_change)+'/step_7_prepare_compound_hierarchy/')
    os.system('touch ../results/'+str(min_fold_change)+'/step_7_prepare_compound_hierarchy/dummy.txt')
    obo_file_address='../resources/classyfire_files/ChemOnt_2_1.obo'
    parsed_obo=obonet.read_obo(obo_file_address)

    for temp_node in parsed_obo.nodes:
        parsed_obo.nodes[temp_node]['type_of_node']='combination'

    pipeline_input_panda_directory='../results/'+str(min_fold_change)+'/step_6_b_generate_signifigance_test_matrices/'
    pipeline_output_directory='../results/'+str(min_fold_change)+'/step_7_prepare_compound_hierarchy/'

    file_list=os.listdir(pipeline_input_panda_directory)
    file_list.remove('dummy.txt')
    
    for file_counter,temp_file in enumerate(file_list):
        if file_counter==0:
            master_panda=pandas.read_pickle(pipeline_input_panda_directory+temp_file)    
            master_panda=master_panda.loc[master_panda['inchikey'] != '@@@@@@@',:]

        else:
            temp_panda=pandas.read_pickle(pipeline_input_panda_directory+temp_file) 
            temp_panda=temp_panda.loc[temp_panda['inchikey'] != '@@@@@@@',:]
            master_panda=pandas.concat(
                [master_panda,temp_panda],axis='index',ignore_index=True
            )

    #get dict 
    #for this dict, the keys are the classes and the values are the chemontid
    #in this way, we use the cfb tool to get the classes for all inchikeys
    #and get the nodes that they connect to with this dict
    class_to_node_dict=make_class_to_node_name_dict(parsed_obo)
    print(class_to_node_dict)
    #hold=input('dict')

    #make a column with the deepest classyfire class possible - in this way we can 
    #find the class to use on the above dict 
    obtain_deepest_classyfire_class_per_bin(master_panda)

    #draw the network before adding the single node
    print(len(parsed_obo.nodes))
    original_classyfire_node_count=len(parsed_obo.nodes)
    #draw classyfire with no adjustments or anything
    #nx.draw(parsed_obo,node_color=['#1f78b4' for i in range(0,len(parsed_obo.nodes))],node_size=50)
    #plt.show()

    add_all_bins_to_network(parsed_obo,master_panda,class_to_node_dict)
    #for each bin assign as a child node of the deepest node that is possible

    remove_branches_without_fold_matrices(parsed_obo)

    #update 220926 plb
    #originally from "reduce hierarchy complexity post dash"
    #compound_nx_address='../results/'+str(min_fold_change)+'/step_8_perform_compound_hierarchical_analysis/classyfire_analysis_results.bin'
    compound_node_keep_address='../resources/species_organ_maps/networkx_shrink_compound.txt'
    parsed_obo=do_everything(parsed_obo,compound_node_keep_address,'compound')
    


    print(master_panda)
    print('*'*50)
    for temp_node in parsed_obo.nodes:
        print(parsed_obo.nodes[temp_node])
    
    nx.readwrite.gpickle.write_gpickle(parsed_obo,pipeline_output_directory+'classyfire_ont_with_bins_added.bin')