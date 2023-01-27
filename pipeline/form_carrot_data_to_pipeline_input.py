#this program takes the output from carrot and transforms it into something pretty close to what we got from the orginal 
#binvestigate pull.
#the basic idea is to obtain the total list of compounds and the total list of (species, organ)
#then we check the soc annotation-only files for each species organ compound
#if the file is found we impute whats missing and add to the input panda
#if the file isnt found we impute the entire list using the noise level
#somewhat unfortunately, retrospectively, the information is stored in "parallel lists" for each column for each bin
#so, the length of species, organ, count, median intensity, annotation distribution, etc should be equal 
import sys
import pandas as pd
import numpy as np
import os
import random

def get_compound_properties(compound_properties_panda):

    a=set(zip(
            compound_properties_panda.loc[compound_properties_panda.target_type=='CONFIRMED'].target_id.astype(int),
            compound_properties_panda.loc[compound_properties_panda.target_type=='CONFIRMED'].name
    ))

    unique_bin_and_name=list(a)

    return unique_bin_and_name


def get_all_species_organs():
    pass

def make_empty_input_panda(pipeline_input_panda_columns,bin_and_name_list):
    temp_dict_of_list={
        temp_column:[np.nan for i in range(len(bin_and_name_list))] for temp_column in pipeline_input_panda_columns
    }

    temp_dict_of_list['id']=[temp[0] for temp in bin_and_name_list]
    temp_dict_of_list['name']=[temp[1] for temp in bin_and_name_list]

    temp_dict_of_list['species']=[list() for temp in bin_and_name_list]
    temp_dict_of_list['organ']=[list() for temp in bin_and_name_list]
    temp_dict_of_list['median_intensity']=[list() for temp in bin_and_name_list]
    temp_dict_of_list['annotation_distribution']=[list() for temp in bin_and_name_list]
    temp_dict_of_list['count']=[list() for temp in bin_and_name_list]
    temp_dict_of_list['total_intensity']=[list() for temp in bin_and_name_list]
    temp_dict_of_list['percent_present']=[list() for temp in bin_and_name_list]

    pipeline_input_panda=pd.DataFrame.from_dict(temp_dict_of_list)

    return pipeline_input_panda

def insert_combination(species,organ,compound):
    pass

def impute_annotations_at_least_one_found(temp_annotations_list,temp_so_count,temp_average_fame_intensity):
    percent_present=len(temp_annotations_list)/temp_so_count
    imputed_value_base=percent_present*min(temp_annotations_list)
    imputed_value_sigma=0.05*imputed_value_base
    return percent_present,temp_annotations_list+[
        (imputed_value_base+random.gauss(0,imputed_value_sigma)) for i in range(temp_so_count-len(temp_annotations_list))
    ]

def impute_annotations_zero_found(temp_average_fame_intensity,temp_so_count):
    #what we impute results from normalizing the imputed value
    #(noise_intensity/1e6)=(imputed/avg fame intensity)
    #by avg fame intensity
    #indeed, what we get is invariant to the fame, which is what we want
    #we want all 0% annotations to have basically the same intensities
    return [((noise_intensity+random.gauss(0,gauss_sigma))/mid_fame_intensity) for i in range(temp_so_count)]

def insert_combination_wrapper(pipeline_input_panda,bin_and_name_list,species_organ_pair_list,data_base_address,species_organ_properties_panda,so_to_skip_set):
    #
    bins_to_remove=list()
    print('inside combination wrapper')
    for temp_counter,temp_bin_and_name in enumerate(bin_and_name_list):
        print(f'we are on directory {temp_counter}')
        print('some row in the input panda')
        print(temp_bin_and_name)
        #we add this to avoid checking groups.
        if len(os.listdir(data_base_address+str(temp_bin_and_name[0])))==0:
            print(data_base_address+str(temp_bin_and_name[0]))
            print('the length was zero')
            bins_to_remove.append(temp_bin_and_name[0])
            continue

        for temp_species_organ_pair in species_organ_pair_list:

            this_species_organ_count=species_organ_properties_panda.loc[
                (species_organ_properties_panda.species==temp_species_organ_pair[0]) &
                (species_organ_properties_panda.organ==temp_species_organ_pair[1])
            ]['count'].item()
            try:
                temp_panda=pd.read_pickle(data_base_address+str(temp_bin_and_name[0])+'/'+temp_species_organ_pair[0]+'¬'+temp_species_organ_pair[1]+'¬'+str(temp_bin_and_name[0])+'.bin')
                temp_panda=temp_panda.to_frame()

                #there is a special case where the annotation is made but the value is still zero
                #why this occurs is a binbase issue
                #if this is the case, then we resort to the noise imputation strategy
                if temp_panda.normalized_intensity.sum()==0:
                    this_species_organ_average_fame_intensity=species_organ_properties_panda.loc[
                        (species_organ_properties_panda.species==temp_species_organ_pair[0]) &
                        (species_organ_properties_panda.organ==temp_species_organ_pair[1])
                    ]['average_fame_intensity'].item()
                    temp_imputed_annotations_list=impute_annotations_zero_found(this_species_organ_average_fame_intensity,this_species_organ_count)
                    temp_percent_present=0
                else:
                    #in general, the lists can have the property where some annotation has a magnitude of zero
                    #the "all zero"        strategy is different, so we keep it by itself above'
                    this_species_organ_average_fame_intensity=species_organ_properties_panda.loc[
                        (species_organ_properties_panda.species==temp_species_organ_pair[0]) &
                        (species_organ_properties_panda.organ==temp_species_organ_pair[1])
                    ]['average_fame_intensity'].item()

                    temp_panda.drop(labels=temp_panda.index[temp_panda.normalized_intensity==0], axis='index',inplace=True)
                    temp_panda.reset_index(inplace=True,drop=True)
                    temp_percent_present,temp_imputed_annotations_list=impute_annotations_at_least_one_found(temp_panda['normalized_intensity'].to_list(),this_species_organ_count,this_species_organ_average_fame_intensity)
                
            except FileNotFoundError:

                this_species_organ_average_fame_intensity=species_organ_properties_panda.loc[
                    (species_organ_properties_panda.species==temp_species_organ_pair[0]) &
                    (species_organ_properties_panda.organ==temp_species_organ_pair[1])
                ]['average_fame_intensity'].item()

                temp_imputed_annotations_list=impute_annotations_zero_found(this_species_organ_average_fame_intensity,this_species_organ_count)
                temp_percent_present=0

            pipeline_input_panda.at[temp_bin_and_name[0],'species'].append(temp_species_organ_pair[0])
            pipeline_input_panda.at[temp_bin_and_name[0],'organ'].append(temp_species_organ_pair[1])
            pipeline_input_panda.at[temp_bin_and_name[0],'count'].append(this_species_organ_count)
            pipeline_input_panda.at[temp_bin_and_name[0],'median_intensity'].append(np.median(temp_imputed_annotations_list))
            pipeline_input_panda.at[temp_bin_and_name[0],'annotation_distribution'].append(temp_imputed_annotations_list)
            pipeline_input_panda.at[temp_bin_and_name[0],'total_intensity'].append(sum(temp_imputed_annotations_list))
            pipeline_input_panda.at[temp_bin_and_name[0],'percent_present'].append(temp_percent_present)


    return pipeline_input_panda,bins_to_remove      


if __name__=="__main__":

    min_fold_change=sys.argv[1]
    #output_pickle_address='../results/'+str(min_fold_change)+'/step_0_b_shape_aws_pull_to_pipeline_input/binvestigate_species_transformed.bin'
    os.system('mkdir -p ../results/'+str(min_fold_change)+'/step_0_b_shape_aws_pull_to_pipeline_input/')
    os.system('touch ../results/'+str(min_fold_change)+'/step_0_b_shape_aws_pull_to_pipeline_input/dummy.txt')
    ##when looking at the fames it was discovered that a number of metadata had "bad fame appearance frequency"
    #that is, for whatever reason, the fames only appeared like 50 percent of the time
    #these are the things in the "elbow plot"
    #this file tells us how to skip systematic bias
    so_to_skip_csv_address='../resources/pull_from_carrot/intermediates/so_to_skip.csv'
    noise_intensity=1
    gauss_sigma=0.05
    #if we dont do this, then we will create artifacts based on the variations in fame intensity
    #histogram the ave fame intensities to see how we got 1M
    #so the imputed value for none found is (noise_intensity/1e6)=(imputed/avg fame intensity)
    #the values of 1 and e6 were chosen by looking at intensities less than 1e6 for bin 4 (valine?)
    #200/1e6, or 1e4, was too big of a signal
    #when we plotted the median intensities of urea, biphenyl, and hexose, 1e6 looked too big.
    #instead now gong with e13
    #then we saw that the total intensity divided by the count didnt have this bimodalityy, and that 1e6
    #was appropriate. so we are going to get rid of median
    #no, before, we were erroneously dividng the already normalized imputed value by the fame intensity, but this was unnecessary
    #everything seems good now
    mid_fame_intensity=1e6
    ##data_base_address='../results/'+str(min_fold_change)+'/step_0_a_pull_distributions_from_aws/soc_data/'
    ##we are not making the pull from carrot part of the chain anymore
    data_base_address='../resources/pull_from_carrot/intermediates/soc_data/'
    pipeline_input_panda_columns=['id','name','species','organ','count','total_intensity','median_intensity','group','inchikey','annotation_distribution','percent_present']
    compound_properties_panda_address='../resources/pull_from_carrot/intermediates/all_species_organs_compounds_panda_2.tsv'
    species_organ_properties_panda_address='../resources/pull_from_carrot/original_from_carrot_input/species_organ_sample_count_and_average_fame.tsv'

    output_address='../results/'+str(min_fold_change)+'/step_0_b_shape_aws_pull_to_pipeline_input/pipeline_input'
    
    so_to_skip_panda=pd.read_csv(so_to_skip_csv_address,sep='¬')
    so_to_skip_set=set(zip(so_to_skip_panda['species'],so_to_skip_panda['organ']))

    compound_properties_panda=pd.read_csv(compound_properties_panda_address,sep='\t')
    species_organ_properties_panda=pd.read_csv(species_organ_properties_panda_address,sep='\t')    


    bin_and_name_list_full=get_compound_properties(compound_properties_panda)
    bin_and_name_list_full=sorted(bin_and_name_list_full)
    print(bin_and_name_list_full)
    print('---------------------')

    #split into 500 bin chunks
    bin_and_name_list_list=[bin_and_name_list_full[x:x+500] for x in range(0, len(bin_and_name_list_full), 500)]

    for bin_and_name_list_counter,bin_and_name_list in enumerate(bin_and_name_list_list):
        species_organ_pair_list=list(zip(species_organ_properties_panda.species,species_organ_properties_panda.organ))

        pipeline_input_panda=make_empty_input_panda(pipeline_input_panda_columns,bin_and_name_list)
        #basically, we discvoered (see notebook step 3) that some groups dont have a bin assocaited with the group name
        #therefore, we made some directories 

        pipeline_input_panda.set_index(keys='id',drop=False,inplace=True)
        
        pipeline_input_panda,bins_to_remove=insert_combination_wrapper(pipeline_input_panda,bin_and_name_list,species_organ_pair_list,data_base_address,species_organ_properties_panda,so_to_skip_set)
        pipeline_input_panda.reset_index(inplace=True,drop=True)

        print(pipeline_input_panda)

        pipeline_input_panda=pipeline_input_panda.loc[~pipeline_input_panda.id.isin(bins_to_remove),:]

        print(pipeline_input_panda)

        pipeline_input_panda.to_pickle(output_address+'_'+str(bin_and_name_list_counter)+'.bin')
