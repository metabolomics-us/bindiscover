import sys
import os
import pandas as pd
import networkx as nx

def choose_all_bins(directory_address):
    full_list=os.listdir(directory_address)
    return full_list

def create_translation_panda_for_compounds(networkx_address,directory_address):
    
    compound_networkx=nx.read_gpickle(networkx_address)
    compound_networkx_nodes_as_str=[str(element) for element in compound_networkx.nodes]
    print(compound_networkx_nodes_as_str)
    print(compound_networkx.nodes)
    full_list=choose_all_bins(directory_address)
    compound_translation_panda=pd.DataFrame.from_dict(
        {
            'compound_identifier':[element[:-4] for i,element in enumerate(full_list)],
            'integer_representation':[i for i,element in enumerate(full_list)]
        }
    )
    [print(compound_translation_panda)]

    bin_type=list()
    for temp_identifier in compound_translation_panda.compound_identifier.to_list():
        if temp_identifier not in compound_networkx_nodes_as_str:
            bin_type.append('unknown')
        else:
            try:
                if compound_networkx.nodes[int(temp_identifier)]['type_of_node']=='from_binvestigate':
                    bin_type.append('known')
                else:
                    bin_type.append('class')
            except ValueError:
                if compound_networkx.nodes[(temp_identifier)]['type_of_node']=='from_binvestigate':
                    bin_type.append('known')
                else:
                    bin_type.append('class')         

    english_name=list()
    for temp_identifier in compound_translation_panda.compound_identifier.to_list():    
        if temp_identifier not in compound_networkx_nodes_as_str:
            english_name.append('Unknown: Bin '+str(temp_identifier))
        else:
            if 'CHEMONTID' in temp_identifier:
                english_name.append(compound_networkx.nodes[temp_identifier]['name'])
            else:
                english_name.append(compound_networkx.nodes[int(temp_identifier)]['common_name'])

    print(english_name)

    identifier=list()
    for temp_identifier in compound_translation_panda.compound_identifier.to_list():    
        if temp_identifier not in compound_networkx_nodes_as_str:
            identifier.append('unknown')
        else:
            if 'CHEMONTID' in temp_identifier:
                identifier.append(temp_identifier)
            else:
                identifier.append(compound_networkx.nodes[int(temp_identifier)]['inchikey'])



    compound_translation_panda['bin_type']=bin_type
    compound_translation_panda['english_name']=english_name
    compound_translation_panda['identifier']=identifier
    print(compound_translation_panda)
    return compound_translation_panda

def create_translation_dict_for_triplets(temp_bin_address):
    temp=pd.read_pickle(temp_bin_address)
    triplet_translation_panda=pd.DataFrame.from_dict(
        {
            'triplet_identifier_tuple':[temp.index[i] for i in range(len(temp.index))],
            'triplet_identifier_string':[temp.index[i][1]+' - '+temp.index[i][0]+' - '+temp.index[i][2] for i in range(len(temp.index))],
            'integer_representation':[i for i in range(len(temp.index))]
        }
    )    

    triplet_translation_panda.set_index('triplet_identifier_string',inplace=True,drop=False)

    #add the count column
    pipeline_input_panda_directory='../results/'+str(min_fold_change)+'/step_5_panda_cleaned/'
    file_list=os.listdir(pipeline_input_panda_directory)
    file_list.remove('dummy.txt')
    binvestigate_panda_address=pipeline_input_panda_directory+file_list[0]
    binvestigate_panda=pd.read_pickle(binvestigate_panda_address)
    organ_list=binvestigate_panda.at[0,'organ']
    species_list=binvestigate_panda.at[0,'species']
    disease_list=binvestigate_panda.at[0,'special_property_list']
    count_list=binvestigate_panda.at[0,'count']
    count_list_dict=dict()
    for i in range(len(count_list)):
        count_list_dict[species_list[i]+' - '+organ_list[i]+' - '+disease_list[i]]=count_list[i]

    print(count_list_dict)
    triplet_translation_panda['count']=triplet_translation_panda['triplet_identifier_string']
    print(triplet_translation_panda)
    triplet_translation_panda['count']=triplet_translation_panda['count'].map(count_list_dict.get)
    print(triplet_translation_panda)
    print('*'*50)

    return triplet_translation_panda

if __name__ == "__main__":

    min_fold_change=sys.argv[1]
    os.system('mkdir -p ../results/'+str(min_fold_change)+'/step_9_generate_extras_for_db_and_api/')
    os.system('touch ../results/'+str(min_fold_change)+'/step_9_generate_extras_for_db_and_api/dummy.txt')

    #make compound translation panda
    compound_networkx_address='../results/'+str(min_fold_change)+'/step_7_prepare_compound_hierarchy/classyfire_ont_with_bins_added.bin'
    compound_panda_directory_address='../results/'+str(min_fold_change)+'/step_8_perform_compound_hierarchical_analysis/all_matrices/fold_change_matrix_average'
    compound_translation_panda=create_translation_panda_for_compounds(compound_networkx_address,compound_panda_directory_address)
    compound_translation_panda.to_pickle('../results/'+str(min_fold_change)+'/step_9_generate_extras_for_db_and_api/compound_translation_panda.bin')
    print(compound_translation_panda)

    #make triplet translation panda
    any_file=os.listdir('../results/'+str(min_fold_change)+'/step_8_perform_compound_hierarchical_analysis/all_matrices/fold_change_matrix_average/')[0]
    triplet_translation_panda=create_translation_dict_for_triplets(
        '../results/'+str(min_fold_change)+'/step_8_perform_compound_hierarchical_analysis/all_matrices/fold_change_matrix_average/'+any_file
    )
    triplet_translation_panda.to_pickle('../results/'+str(min_fold_change)+'/step_9_generate_extras_for_db_and_api/triplet_translation_panda.bin')
    print(triplet_translation_panda)