#remember that this does more than preapre the species networkx
#it also swaps the species names with the ncbi taxids to make it so that all the different sources
#work (ete3, metagenompy, our curation)
#this swap occurs in the binvestigate panda

from cProfile import label
import os
import pandas
import metagenompy
import matplotlib.pyplot as plt
import networkx as nx
from pprint import pprint
import pydot
from networkx.drawing.nx_pydot import graphviz_layout
import sys

def get_all_strings_in_list_across_panda_column(temp_panda,temp_column_name):
    '''
    this is a dupe of the funciton foudn in transform written species to ncbi species
    couldnt import because of conda versions, ete3, etc
    '''
    total_elements=set()

    for index, series in temp_panda.iterrows():
        temp_set=set(series[temp_column_name])
        total_elements.update(temp_set)

    return total_elements

def swap_panda_column_values_per_mapping(binvestigate_panda,species_taxid_panda):
    '''
    '''
    swapping_dict={pair[0]:pair[1] for pair in zip(species_taxid_panda['species'],species_taxid_panda['tax_id'].astype(str))}

    for index, series in binvestigate_panda.iterrows():
        temp_list=list(map(swapping_dict.get,series['species']))
        binvestigate_panda.at[index,'species']=temp_list

def contract_irrelevant_nodes(temp_nx,temp_binvestigate_taxid_list):
    #add node to list if (its name is NOT in the species taxid panda) AND (it has <2 children)
    #when those two conditions are met, then the node's species results can be inferred as being exactly that of the child

    #the contract function takes pairs of nodes
    #we find the predeccessor dynamically because the graph changes with evey contraction 
    #binvestigate_taxid_list=temp_panda['tax_id'].astype(str).to_list()
    
    nodes_to_contract=list()

    #i believe that the headnode is kept only because it is actually self referential. i believe that
    #we should probably have another condition (len(predecessors) != 0)
    for temp_node in temp_nx.nodes:
        if (temp_node not in temp_binvestigate_taxid_list) and (len(list(temp_nx.successors(temp_node)))<2):
            nodes_to_contract.append(temp_node)

    #first node is kept. so the predecessor of each to remove
    for temp_node in nodes_to_contract:
        temp_edge=temp_nx.add_edge(list(temp_nx.predecessors(temp_node))[0],list(temp_nx.successors(temp_node))[0])
        temp_nx.remove_node(temp_node)

    #remove all self loops. artifact that has no meaning
    for temp_node in temp_nx.nodes:
        try:
            temp_nx.remove_edge(temp_node,temp_node)
        except nx.exception.NetworkXError:
            continue


def visualize_nodes_on_a_list(temp_nx,temp_list,temp_attribute_name):
    '''
    temp_attribute_name=None when you dont want to use an attribute name, rather the 
    node names themselves
    '''
    if temp_attribute_name != None:
        color_list=list()
        for temp_node in temp_nx.nodes:
            if temp_nx.nodes[temp_node][temp_attribute_name] in temp_list:
                color_list.append('#ffa500')
            else:
                color_list.append('#00ff00')
        pos = graphviz_layout(temp_nx, prog="dot")
        labels = nx.get_node_attributes(temp_nx, 'scientific_name') 
        nx.draw(temp_nx, pos,with_labels=labels,node_color=color_list)
        plt.show()    

    elif temp_attribute_name == None:
        color_list=list()
        for temp_node in temp_nx.nodes:
            if temp_node in temp_list:
                color_list.append('#ffa500')
            else:
                color_list.append('#00ff00')
        pos = nx.nx_agraph.pygraphviz_layout(temp_nx, prog='dot')

        labels = nx.get_node_attributes(temp_nx, 'scientific_name') 
        print(labels)
        nx.draw(temp_nx, pos,labels=labels)#,node_color=color_list)
        plt.show()  



if __name__ == "__main__":

    #min_fold_change=10
    min_fold_change=sys.argv[1]
    #the input panda here is the taxid<->species mapping.
    #we do not need to get a taxonomy mapping because metagenompy prepares that automatically
    #from a local database
    input_mapping_panda_address='../results/'+str(min_fold_change)+'/step_8_a_create_species_taxid_mapping/species_tax_id_mapping.bin'
    
    pipeline_input_panda_directory='../results/'+str(min_fold_change)+'/step_6_generate_fold_matrices/'
    file_list=os.listdir(pipeline_input_panda_directory)
    file_list.remove('dummy.txt')
    input_binvestigate_panda_address=pipeline_input_panda_directory+file_list[0]
    
    output_networkx_address='../results/'+str(min_fold_change)+'/step_8_b_prepare_species_networkx/species_networkx.bin'
    output_binvestigate_panda_address='../results/'+str(min_fold_change)+'/step_8_b_prepare_species_networkx/binvestigate_species_as_taxid.bin'
    os.system('mkdir -p ../results/'+str(min_fold_change)+'/step_8_b_prepare_species_networkx/')
    os.system('touch ../results/'+str(min_fold_change)+'/step_8_b_prepare_species_networkx/dummy.txt')

    #read in things or make the species network from databse
    #plb edit 2-6-2022
    #even tho we added a signifigance testing step, we only get the species form step six, so we can leave "step_6" above alone
    binvestigate_panda=pandas.read_pickle(input_binvestigate_panda_address)
    total_ncbi_networkx = metagenompy.generate_taxonomy_network(auto_download=True)
    species_taxid_panda=pandas.read_pickle(input_mapping_panda_address)

    #get a list of taxids for each species foudn in this binvestigate panda
    #all species in the google sheets transform
    #shrinkt the total possibel list to those that appear in this particular binvestigate panda
    species_in_this_binvestigate_panda_set=get_all_strings_in_list_across_panda_column(binvestigate_panda,'species')
    taxid_for_this_binvestigate_panda_list=list()
    
    output_for_species_translation={
        'english':[],
        'ncbi_id':[]
    }
    for temp_species in species_in_this_binvestigate_panda_set:
        print(temp_species)
        print(species_taxid_panda.loc[species_taxid_panda['species']==temp_species,'tax_id'].values[0])
        taxid_for_this_binvestigate_panda_list.append(str(species_taxid_panda.loc[species_taxid_panda['species']==temp_species,'tax_id'].values[0]))
        output_for_species_translation['english'].append(temp_species)
        output_for_species_translation['ncbi_id'].append(species_taxid_panda.loc[species_taxid_panda['species']==temp_species,'tax_id'].values[0])

    pandas.DataFrame.from_dict(output_for_species_translation).to_pickle('../results/'+str(min_fold_change)+'/step_8_b_prepare_species_networkx/for_index_panda_for_dash_species_translation.bin')
    #shrink the entire ncbi networkx according to the sublist of nodes 
    #create a copy that we work with
    binvestigate_species_networkx_view=metagenompy.highlight_nodes(total_ncbi_networkx,taxid_for_this_binvestigate_panda_list)
    binvestigate_species_networkx=binvestigate_species_networkx_view.copy()

    #back when we thought that the classic King Philliip Went To Greece etc set would work
    #now we must include instraclasses and chordates and blah blah
    #metagenompy.condense_taxonomy(binvestigate_species_networkx)

    #render the networkx before further work is done
    # fig, ax = plt.subplots(figsize=(10, 10))
    # metagenompy.plot_network(binvestigate_species_networkx, ax=ax, labels_kws=dict(font_size=10))
    # fig.tight_layout()
    # plt.show()
    
    print('------------------------------------------------')
    for temp_taxid in taxid_for_this_binvestigate_panda_list:
        if temp_taxid not in binvestigate_species_networkx.nodes:
            print(temp_taxid)
    print('------------------------------------------------')
    
    #at this point i realized that the species labels from ete3 were not playing nicely with the species labels from metagenompy
    #however, the taxid were
    #i also realized that even when using the taxid, stupid things like "selachii" aka sharks being technically an "intraclass"
    #prevented them from being mapped onto the simplified taxonomy
    #therefore, we are going to switch the species names with taxid using the "mapping file produced in the previous step"
    #also, we are going to reduce the networkx to just nodes that have direct entries or are branches. this will greatly greatly 
    #speed the computation. later, if desired, we can fill in nodes directly copying the results of child nodes
    #switch the species names with taxid values in binvestigate panda
    swap_panda_column_values_per_mapping(binvestigate_panda,species_taxid_panda)
    
    #read in cleaned panda to see how many species dont map to somewhere
    all_species_set=get_all_strings_in_list_across_panda_column(binvestigate_panda,'species')
    for temp_species in all_species_set:
        if temp_species not in binvestigate_species_networkx.nodes:
            print(type(temp_species))
            print(temp_species)

    for temp_node in binvestigate_species_networkx.nodes:
        binvestigate_species_networkx.nodes[temp_node]['ncbi_number']=temp_node


    #render the network after node contraction
    node_labels = {
        n: data.get('scientific_name', '') for n, data in binvestigate_species_networkx.nodes(data=True)
    }

    #render the network with the "binvestigate nodes" being a special color
    species_id_set=get_all_strings_in_list_across_panda_column(binvestigate_panda,'species')
    # visualize_nodes_on_a_list(binvestigate_species_networkx,species_id_set,None)

    #print the network and the binvestigate panda with the species as taxid
    nx.readwrite.gpickle.write_gpickle(binvestigate_species_networkx,output_networkx_address)
    binvestigate_panda.to_pickle(output_binvestigate_panda_address)
