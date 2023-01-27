import sys
import os
import networkx as nx
import pandas as pd


def make_index_panda_for_dash_app():
    starting_point=pd.read_pickle('../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/sod_combo.bin')
    starting_point['species_english']=starting_point['species']
    starting_point['organ_english']=starting_point['organ']
    starting_point['disease_english']=starting_point['disease']
    species_mapping_panda=pd.read_pickle('../results/'+str(min_fold_change)+'/step_8_b_prepare_species_networkx/for_index_panda_for_dash_species_translation.bin')

    disease_networkx=nx.read_gpickle('../results/'+str(min_fold_change)+'/step_8_c_prepare_organ_and_disease_networkx/disease_networkx.bin')
    species_mapping_dict=dict(zip(species_mapping_panda.english,species_mapping_panda.ncbi_id))
    starting_point['species']=starting_point['species'].map(species_mapping_dict.get)

    organ_networkx=nx.read_gpickle('../results/'+str(min_fold_change)+'/step_8_c_prepare_organ_and_disease_networkx/organ_networkx.bin')
    keys=set()
    for temp_node in organ_networkx.nodes:
        keys.add(organ_networkx.nodes[temp_node]['mesh_label'])
    organ_mapping_dict={temp_key:list() for temp_key in keys}
    for temp_node in organ_networkx.nodes:
        organ_mapping_dict[organ_networkx.nodes[temp_node]['mesh_label']].append(temp_node)
    starting_point['organ']=starting_point['organ'].map(organ_mapping_dict.get)
    starting_point=starting_point.explode('organ')
    disease_networkx=nx.read_gpickle('../results/'+str(min_fold_change)+'/step_8_c_prepare_organ_and_disease_networkx/disease_networkx.bin')
    keys=set()
    for temp_node in disease_networkx.nodes:
        keys.add(disease_networkx.nodes[temp_node]['mesh_label'])
    disease_mapping_dict={temp_key:list() for temp_key in keys}
    for temp_node in disease_networkx.nodes:
        disease_mapping_dict[disease_networkx.nodes[temp_node]['mesh_label']].append(temp_node)
    starting_point['disease']=starting_point['disease'].map(disease_mapping_dict.get)
    starting_point=starting_point.explode('disease')
    starting_point.set_index(keys=['organ','species','disease'],drop=False,inplace=True)
    starting_point.to_pickle('../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/index_panda.bin')


def get_unique_sod_combinations():
    unique_sod_combinations_address = "../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/unique_sod_combinations.bin"
    unique_sod_combinations_panda = pd.read_pickle(unique_sod_combinations_address)
    unique_sod_combinations_dict = {
        temp:temp for temp in unique_sod_combinations_panda.keys().to_list()
    }
    return unique_sod_combinations_dict

def extract_networkx_selections_species():
    species_networkx=nx.read_gpickle('../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/species_networkx.bin')
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
        species_node_dict['9606']='Homo sapiens AKA human'

    return species_networkx,species_node_dict

def get_index_panda():
    temp=pd.read_pickle('../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/index_panda.bin')
    return temp

if __name__ == "__main__":

    min_fold_change=sys.argv[1]
    os.system('mkdir -p ../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/')
    os.system('touch ../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/dummy.txt')

    #add non ratio stuff
    non_ratio_dropdown_address='../results/'+str(min_fold_change)+'/step_5_b_make_non_ratio_table/unique_sod_combinations.bin'
    non_ratio_dropdown_address_output='../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/unique_sod_combinations.bin'    
    os.system(f'cp {non_ratio_dropdown_address} {non_ratio_dropdown_address_output}')

    #networkxs
    compound_networkx_address='../results/'+str(min_fold_change)+'/step_7_prepare_compound_hierarchy/classyfire_ont_with_bins_added.bin'
    compound_networkx_address_output='../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/compounds_networkx.bin'
    compound_networkx=nx.readwrite.read_gpickle(compound_networkx_address)
    nx.readwrite.write_gpickle(compound_networkx,compound_networkx_address_output)

    species_networkx_address='../results/'+str(min_fold_change)+'/step_8_b_prepare_species_networkx/species_networkx.bin'
    species_networkx_address_output='../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/species_networkx.bin'    
    os.system(f'cp {species_networkx_address} {species_networkx_address_output}')

    organ_networkx_address='../results/'+str(min_fold_change)+'/step_8_c_prepare_organ_and_disease_networkx/organ_networkx.bin'
    organ_networkx_address_output='../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/organ_networkx.bin'    
    os.system(f'cp {organ_networkx_address} {organ_networkx_address_output}')

    disease_networkx_address='../results/'+str(min_fold_change)+'/step_8_c_prepare_organ_and_disease_networkx/disease_networkx.bin'
    disease_networkx_address_output='../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/disease_networkx.bin'    
    os.system(f'cp {disease_networkx_address} {disease_networkx_address_output}')

    #an index of a fold matrix so that we can choose subsets of S,O,D on the graphical tool
    random_fold_matrix_for_index_directory_file_list=os.listdir('../results/'+str(min_fold_change)+'/step_8_perform_compound_hierarchical_analysis/all_matrices/fold_change_matrix_average/')
    random_fold_matrix_for_index_address='../results/'+str(min_fold_change)+'/step_8_perform_compound_hierarchical_analysis/all_matrices/fold_change_matrix_average/'+random_fold_matrix_for_index_directory_file_list[0]
    temp=pd.read_pickle(random_fold_matrix_for_index_address)
    output_panda=temp.index.to_frame()
    print(output_panda)
    output_panda.to_pickle('../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/sod_combo.bin')

    triplet_mapping_address='../results/'+str(min_fold_change)+'/step_9_generate_extras_for_db_and_api/triplet_translation_panda.bin'
    triplet_mapping_output='../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/triplet_translation_panda.bin' 
    os.system(f'cp {triplet_mapping_address} {triplet_mapping_output}')

    compound_mapping_address='../results/'+str(min_fold_change)+'/step_9_generate_extras_for_db_and_api/compound_translation_panda.bin'
    compound_mapping_output='../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/compound_translation_panda.bin' 
    os.system(f'cp {compound_mapping_address} {compound_mapping_output}')


    index_panda=make_index_panda_for_dash_app()


    species_translation_address='../results/'+str(min_fold_change)+'/step_8_b_prepare_species_networkx/for_index_panda_for_dash_species_translation.bin'
    species_translation_output='../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/for_index_panda_for_dash_species_translation.bin'
    os.system(f'cp {species_translation_address} {species_translation_output}')

    #after transferring, we glue on this step in order to make the clickable labels for differential + upset
    unique_sod_combinations_address = "../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/unique_sod_combinations.bin"
    _,species_code_to_english_dict_map=extract_networkx_selections_species()
    unique_sod_combinations_panda = pd.DataFrame(
        pd.read_pickle(unique_sod_combinations_address)
    )
    unique_sod_combinations_panda.reset_index(inplace=True)
    output_panda=unique_sod_combinations_panda['index'].str.split(' - ',expand=True)
    species_english_to_code_map=dict(zip(
        index_panda['species_english'],index_panda['species']
    ))
    output_panda[0]=output_panda[0].map(species_english_to_code_map.get)
    output_panda[0]=output_panda[0].astype(str).map(species_code_to_english_dict_map.get)

    output_panda['triplet_with_common_name']=output_panda[[0,1,2]].agg(' - '.join,axis=1)
    output_panda['triplet']=unique_sod_combinations_panda['index']
    output_panda.drop([0,1,2],axis='columns',inplace=True)
    output_panda.to_pickle('../results/'+str(min_fold_change)+'/step_11_organize_files_for_dash_app/unique_sod_combinations_common_and_database_names.bin')