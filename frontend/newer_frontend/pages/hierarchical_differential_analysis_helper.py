import networkx as nx

'''
together these functions handle the network analysis that dynamically selects species, organs, and diseases
'''

def extract_networkx_selections_species():
    species_networkx=nx.read_gpickle('../newer_datasets/species_networkx.bin')
    species_node_dict=dict()
    #when we make the options dict, the strategy depends on what common nodes are available
    #note that some have more than one name, so we just choose the first listed
    #so we check if its a string. if its not, then it is a list
    for temp_node in species_networkx.nodes:
        if temp_node=='1':
            species_node_dict[temp_node]='All Species'
        elif 'common_name' in species_networkx.nodes[temp_node].keys():
            if isinstance(species_networkx.nodes[temp_node]['common_name'],str):
                species_node_dict[temp_node]=species_networkx.nodes[temp_node]['scientific_name']+' AKA '+species_networkx.nodes[temp_node]['common_name']
            else:
                species_node_dict[temp_node]=species_networkx.nodes[temp_node]['scientific_name']+' AKA '+species_networkx.nodes[temp_node]['common_name'][0]
        elif 'genbank_common_name' in species_networkx.nodes[temp_node].keys():
            if isinstance(species_networkx.nodes[temp_node]['genbank_common_name'],str):
                species_node_dict[temp_node]=species_networkx.nodes[temp_node]['scientific_name']+' AKA '+species_networkx.nodes[temp_node]['genbank_common_name']
            else:
                species_node_dict[temp_node]=species_networkx.nodes[temp_node]['scientific_name']+' AKA '+species_networkx.nodes[temp_node]['genbank_common_name'][0]
        else:
            species_node_dict[temp_node]=species_networkx.nodes[temp_node]['scientific_name']
        #humans has an unusual common name 
        species_node_dict['9606']='Homo sapiens AKA human'

    return species_networkx,species_node_dict
    

def extract_networkx_selections_organ():
    organ_networkx=nx.read_gpickle('../newer_datasets/organ_networkx.bin')
    organ_node_dict=dict()
    for temp_node in organ_networkx.nodes:
        if temp_node=='organ':
            organ_node_dict[temp_node]='All Organs'
        else:
            organ_node_dict[temp_node]=organ_networkx.nodes[temp_node]['mesh_label']+' - '+temp_node

    return organ_networkx,organ_node_dict

def extract_networkx_selections_disease():
    disease_networkx=nx.read_gpickle('../newer_datasets/disease_networkx.bin')
    disease_node_dict=dict()
    for temp_node in disease_networkx.nodes:
        if temp_node=='disease':
            disease_node_dict[temp_node]='All diseases'
        else:
            disease_node_dict[temp_node]=disease_networkx.nodes[temp_node]['mesh_label']+' - '+temp_node

    return disease_networkx,disease_node_dict