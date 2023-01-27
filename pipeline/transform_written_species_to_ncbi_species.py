import numpy as np
import pandas
from ete3 import NCBITaxa
import os
import multiprocessing
import sys
from sqlalchemy import create_engine
from sqlalchemy import Table, String
from sqlalchemy.dialects import postgresql
import re

def get_all_strings_in_list_across_panda_column(temp_panda,temp_column_name):
    '''
    in the binvestigate panda, each row is a bin
    each bin is associated with organs and species that are saved as lists
    
    this function iterates over each rows lists and attempts to add all elements to a cumulative set
    '''
    total_elements=set()

    for index, series in temp_panda.iterrows():

        temp_set=set(series[temp_column_name])

        total_elements.update(temp_set)

    return total_elements

def identify_elements_not_in_ncbi(temp_set):
    '''
    takes a set of strings and returns a set that do not map to exactly 1 taxonomical identification
    '''

    #creat an ncbi taxonomy object. theres a sql database hidden in the home directory that this references
    ncbi=NCBITaxa()
    has_0_id_set=set()
    has_more_than_1_id_set=set()

    #we scroll over each instead of sending the entire set because those that do not map are not returned in the dict
    for temp_element in temp_set:

        returned_dict=ncbi.get_name_translator([temp_element])

        #if the dict is emtpy
        if not any(returned_dict):
            has_0_id_set.add(temp_element)

        #if the string returns more than 1 element
        elif len(returned_dict[temp_element]) > 1:
            has_more_than_1_id_set.add(temp_element)

        #otherwise do nothing, we only care about problems

    return has_0_id_set, has_more_than_1_id_set

def print_some_set(temp_set):
    '''
    prints elements of a set, 1 element per line
    '''
    for temp_element in temp_set:
        print(temp_element)

def transform_species_column(temp_bin_panda):
    '''
    This function takes the binvestigate panda and the address to a mapping .tsv

    It applies the rules from the mapping tsv to the binvestigate panda

    There are two types of rules: drop rules and transform rules
    Drop rules exist because we dont want all "species" in our final analysis. that is because
    the species is something nonsensical like "nist standards" or "delete me"
    In drop rules, we find the index of all species equal to the original, then drop the species at those indices
    and the corresponding organs at those indices

    in transfomrs, we take the original text and map it to the most specific thing that we can that is in the ncbi database
    '''
    temp_mapping_address=species_mapping_address
    mapping_panda=pandas.read_csv(temp_mapping_address,sep='@')

    #apply each transformation to each row
    for mapping_index, mapping_series in mapping_panda.iterrows():

        #declare mapping rule
        #drop rules work differently. must get all indices and then drop from mapping as well as organs
        if mapping_series['most_specific'] == 'drop' or mapping_series['most_specific'] == 'Drop':
            
            for bin_index, bin_series in temp_bin_panda.iterrows():
                indices_to_drop=[i for i in range(0,len(bin_series['species'])) if bin_series['species'][i] == mapping_series['list_of_species_that_had_zero_ncbi_id']]
                species_list_with_indices_removed=list(np.delete(np.array(bin_series['species'],dtype=object),indices_to_drop))
                organ_list_with_indices_removed=list(np.delete(np.array(bin_series['organ'],dtype=object),indices_to_drop))
                total_intensity_list_with_indices_removed=list(np.delete(np.array(bin_series['total_intensity'],dtype=object),indices_to_drop))
                median_intensity_list_with_indices_removed=list(np.delete(np.array(bin_series['median_intensity'],dtype=object),indices_to_drop))
                count_list_with_indices_removed=list(np.delete(np.array(bin_series['count'],dtype=object),indices_to_drop))
                annotation_distribution_list_with_indices_removed=list(np.delete(np.array(bin_series['annotation_distribution'],dtype=object),indices_to_drop))
                percent_present_list_with_indices_removed=list(np.delete(np.array(bin_series['percent_present'],dtype=object),indices_to_drop))

                temp_bin_panda.at[bin_index,'species']=species_list_with_indices_removed
                temp_bin_panda.at[bin_index,'organ']=organ_list_with_indices_removed
                temp_bin_panda.at[bin_index,'total_intensity']=total_intensity_list_with_indices_removed
                temp_bin_panda.at[bin_index,'median_intensity']=median_intensity_list_with_indices_removed
                temp_bin_panda.at[bin_index,'count']=count_list_with_indices_removed
                temp_bin_panda.at[bin_index,'annotation_distribution']=annotation_distribution_list_with_indices_removed
                temp_bin_panda.at[bin_index,'percent_present']=percent_present_list_with_indices_removed


        #if we have a transformation
        elif mapping_series['most_specific'] != 'drop' or mapping_series['most_specific'] != 'Drop':

            for bin_index, bin_series in temp_bin_panda.iterrows():

                for i in range(0,len(bin_series['species'])):
                    if (bin_series['species'][i].lower()) == (mapping_series['list_of_species_that_had_zero_ncbi_id'].lower()):
                        bin_series['species'][i] = mapping_series['most_specific'].lower()
                    
                temp_bin_panda.at[bin_index,'species']=bin_series['species']
                
    #quick fix to make everythin the same case
    for index,series in temp_bin_panda.iterrows():
        temp_bin_panda.at[index,'species']=[element.lower() for element in temp_bin_panda.at[index,'species']]

    return temp_bin_panda

if __name__ == "__main__":
    min_fold_change=sys.argv[1]
    cores_available=int(sys.argv[2])
    start_from_aws='False'
    species_mapping_address='../resources/species_organ_maps/species_map.txt'

    pipeline_input_panda_directory='../results/'+str(min_fold_change)+'/step_0_c_complete_pipeline_input/'
    pipeline_output_directory='../results/'+str(min_fold_change)+'/step_1_species_transformed/'
    os.system('mkdir -p ../results/'+str(min_fold_change)+'/step_1_species_transformed/')
    os.system('touch ../results/'+str(min_fold_change)+'/step_1_species_transformed/dummy.txt')
    
    file_list=os.listdir(pipeline_input_panda_directory)
    file_list.remove('dummy.txt')
    pandas_list=list()        
    
    for temp_file_counter,temp_file in enumerate(file_list):
        print(f'we are startin to compute file {temp_file_counter}')
        print(f'we are starting to compute file {temp_file}')

        binvestigate_panda=pandas.read_pickle(pipeline_input_panda_directory+temp_file)

        #identify species names that do not map to the NCBI database
        #this is done one time not as part of the "informatics pipeline" but rather
        #so that we can build the "species_mapping.txt"
        all_species_set=get_all_strings_in_list_across_panda_column(binvestigate_panda,'species')
        #we want species to map to exactly 1 place in the taxonomy. sometimes
        #species map to 0 and sometimes they map to >1
        #the "id" in this case is a number associated with the ncbi taxonomy
        has_0_id_set,has_more_than_1_id_set=identify_elements_not_in_ncbi(all_species_set)
        
        #later, after doing a curatin of this list, we put the transform file in the 
        #'species_mapping_address' location
        #translate species names according to some transform list and drop things that either have no translation or are identified with a translation 'drop'
        num_processes= cores_available-1
        chunk_size = len(binvestigate_panda.index)//num_processes
        panda_chunks=list()
        for i in range(0,num_processes):
            panda_chunks.append(binvestigate_panda.iloc[i*chunk_size+i:(i+1)*chunk_size+i+1])
        pool = multiprocessing.Pool(processes=cores_available)
        transformed_chunks=pool.map(transform_species_column,panda_chunks)
        pool.close()
        pool.join()

        #recombine_chunks
        for i in range(len(transformed_chunks)):
            binvestigate_panda.iloc[transformed_chunks[i].index]=transformed_chunks[i]
        binvestigate_panda=pandas.concat(transformed_chunks)
        

        #confirm that has_0_id_set,has_more_than_1_id_set are empty
        all_species_set=get_all_strings_in_list_across_panda_column(binvestigate_panda,'species')
        has_0_id_set,has_more_than_1_id_set=identify_elements_not_in_ncbi(all_species_set)
        print_some_set(has_0_id_set)
        print_some_set(has_more_than_1_id_set)

        #output the result
        temporary_file_integer=re.findall(r'\d+', temp_file)[0]

        print('we are about to write to file')
        binvestigate_panda.to_pickle(pipeline_output_directory+'binvestigate_species_transformed_'+str(temporary_file_integer)+'.bin')