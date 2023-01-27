import numpy as np
import pandas
import os
#import transform_written_organs
import multiprocessing
from pprint import pprint
import sys
from itertools import repeat
import scipy.stats
import re


def show_all_organ_species_disease_triplets(temp_panda):
    set_of_organ_species_disease_tuples=set()
    for index, series in temp_panda.iterrows():
        this_bins_triplets=zip(series['organ'],series['species'],series['special_property_list'])
        for this_bins_triplets in this_bins_triplets:
            set_of_organ_species_disease_tuples.add(this_bins_triplets)
    return set_of_organ_species_disease_tuples


def calculate_one_signifigance_matrix_trip(temp_bin,temp_MultiIndex,signifigance_type,organ_species_disease_tuple_list):
    '''
    asdf
    '''

    #this is the fold change matrix that we start with
    temp_DataFrame=pandas.DataFrame(data=np.nan,index=temp_MultiIndex,columns=temp_MultiIndex)
    tuple_list=zip(temp_bin['organ'],temp_bin['species'],temp_bin['special_property_list'])
    distribution_dict=dict(zip(tuple_list,temp_bin['annotation_distribution']))

    #we iterate through the rows in the fold change matrix
    #we couldnt do neat .apply or other vectorized approaches because the data required
    #were in lists for each bin (series)
    #plb edit 2-7-2022
    #we put the if statement on the outside so we do it once not 1 billion times
    if signifigance_type=='mannwhitney':
        for index,series in temp_DataFrame.iterrows():
            from_distribution=distribution_dict[series.name]
            for temp_column in temp_DataFrame.columns:
                #if we are on a diagonal
                if index == temp_column:
                    temp_DataFrame.at[series.name,temp_column]=np.nan
                    continue
                else:
                    #placeholder while scipy is down
                    _,p=scipy.stats.mannwhitneyu(from_distribution,distribution_dict[temp_column],use_continuity=False,alternative='two-sided')
                    temp_DataFrame.at[series.name,temp_column]=p
    elif signifigance_type=='welch':
        for index,series in temp_DataFrame.iterrows():
            from_distribution=distribution_dict[series.name]
            for temp_column in temp_DataFrame.columns:
                #if we are on a diagonal
                if index == temp_column:
                    temp_DataFrame.at[series.name,temp_column]=np.nan
                    continue
                else:
                    #placeholder while scipy is down
                    _,p=scipy.stats.ttest_ind(np.log10(np.array(from_distribution)),np.log10(np.array(distribution_dict[temp_column])),equal_var=False,alternative='two-sided')
                    temp_DataFrame.at[series.name,temp_column]=p
    
    #plb 2-7-2022
    #some of these things are so signifigant that they are numerically incalculable
    #like, the p value is zero
    #therefore, we impute the min p value for that dataframe
    temp=temp_DataFrame.stack().stack().stack()
    impute_value=temp.loc[temp != 0].min()
    temp_DataFrame=temp_DataFrame.applymap(lambda x: impute_value if x==0 else x)
    
    return temp_DataFrame


def calculate_all_signifigance_matrices_trip(temp_panda,signifigance_type,organ_species_disease_tuple_list):
    '''
    asdf
    '''

    temp_organ_species_disease_tuple_list=organ_species_disease_tuple_list

    temp_panda['signifigance_'+signifigance_type]='pre_analysis'

    #we use multiindex so that we can do cute things switching the ordering later
    my_MultiIndex=pandas.MultiIndex.from_tuples(tuples=temp_organ_species_disease_tuple_list,sortorder=None,names=['organ','species','disease'])

    for index,series in temp_panda.iterrows():
        #we print merely to see how long this is taking
        print(index)
        print(series['name'])
        temp_panda.at[index,'signifigance_'+signifigance_type]=calculate_one_signifigance_matrix_trip(series,my_MultiIndex,signifigance_type,organ_species_disease_tuple_list)    

    return temp_panda


if __name__ == "__main__":
    multiprocessing.set_start_method("spawn")
    #2-06-2022 plb
    #deleted the "non trip version" of the functions as they only seemed to handle species/organ

    #min_fold_change=snakemake.params.min_fold_change
    min_fold_change=sys.argv[1]
    cores_available=int(sys.argv[2])

    input_panda_file=sys.argv[3]
    temp_file_numbers=temporary_file_integer=re.findall(r'\d+', input_panda_file)
    input_panda_address='../results/'+str(min_fold_change)+'/step_6_generate_fold_matrices/'+input_panda_file
    output_panda_address='../results/'+str(min_fold_change)+'/step_6_b_generate_signifigance_test_matrices/'+\
        'binvestigate_with_signifigance_matrices_'+str(temp_file_numbers[0])+'_'+str(temp_file_numbers[1])+'.bin'

    input_panda=pandas.read_pickle(input_panda_address)
    organ_species_disease_tuple_list=list(show_all_organ_species_disease_triplets(input_panda))
    #all fold change matrices have the same row/column labels (see single for further logic)
    organ_species_disease_tuple_list.sort(key=lambda temp_tup: (temp_tup[0],temp_tup[1],temp_tup[2]))
    
    #update 7-4-22 plb
    #basically, we are only going to do the volcano plot stuff for knowns.
    #so we grab a subset of the entire panda, those with inchikeys, do the comparison for them
    #and then merge
    #update 220926 we are going to try for all of the bins
    #so now, only_identified is a misnomer. felt easier than rewriting everything
    input_panda_only_identified=input_panda.copy()
    
    ####
    temp_signifigance_type='mannwhitney'
    
    input_panda_only_identified=calculate_all_signifigance_matrices_trip(input_panda_only_identified,temp_signifigance_type,organ_species_disease_tuple_list)
    ####

    ####
    temp_signifigance_type='welch'
    
    input_panda_only_identified=calculate_all_signifigance_matrices_trip(input_panda_only_identified,temp_signifigance_type,organ_species_disease_tuple_list)

    input_panda_only_identified.to_pickle(output_panda_address)