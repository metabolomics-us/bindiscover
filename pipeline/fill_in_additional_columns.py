#certain properties arent accounted for in the carrot database
#things like "groups" or "inchikeys" or "whether a bin "is a real compound""
#this script adds those properties to make them congruent with the input to the pipeline

#the strategy is to read through every line in the additional property panda and add group/inchikey info as available
#to the input panda
#later down the pipeline, if there is no inchikey, then the compound will be connected to root in the networkx?
import pandas as pd
import numpy as np
import sys
import os
import re

def input_addtional_properties(pipeline_input_panda, additional_property_panda):

    #set the index of each to be the bin id momentarily
    #iterate through the additional property panda, putting values into the pipeline input panda
    #reset the index
    #dont need to do this for one actually, but in a rush so leave for now

    pipeline_input_panda.set_index(keys='id',drop=False,inplace=True)
    additional_property_panda.set_index(keys='bin_id',drop=False,inplace=True)

    pipeline_input_panda['group']=pipeline_input_panda['group'].astype(str)
    pipeline_input_panda['inchikey']=pipeline_input_panda['inchikey'].astype(str)

    for index, series in additional_property_panda.iterrows():
        #in case our input file has extra bins that our carrot data does not
        if index in pipeline_input_panda.index:
            pipeline_input_panda.at[index,'group']=series['group_id']
            pipeline_input_panda.at[index,'inchikey']=series['inchi_key']

    pipeline_input_panda.reset_index(inplace=True,drop=True)
    additional_property_panda.reset_index(inplace=True,drop=True)

    pipeline_input_panda['group'].replace('nan',np.nan,inplace=True)
    pipeline_input_panda['inchikey'].replace('nan',np.nan,inplace=True)

    return

def concat_many_pipeline_inputs(pipeline_input_panda_directory,named,addtional_property_csv_address):
    '''
    takes the output of step 0_b and concats it into one big file
    local desktop cant handle ram of all compounds, hence named property
    so if named property is true, read in all the pandas, and maintain only those with inchikeys 
    in the harmonized bin group inchi pandas
    '''

    file_list=os.listdir(pipeline_input_panda_directory)
    file_list.remove('dummy.txt')


    if named=='only_named':
        additional_property_panda=pd.read_csv(addtional_property_csv_address)
        bins_of_interest=additional_property_panda.loc[additional_property_panda.inchi_key.isnull()==False].bin_id.tolist()
        pandas_list=list()        
        for temp_file in file_list:
            temp_panda=pd.read_pickle(pipeline_input_panda_directory+temp_file)
            temp_panda=temp_panda.loc[
                temp_panda['id'].isin(bins_of_interest)
            ]
            pandas_list.append(temp_panda)

    elif named=='all':
        pandas_list=list()        
        for temp_file in file_list:
            temp_panda=pd.read_pickle(pipeline_input_panda_directory+temp_file)
            pandas_list.append(temp_panda)       


    total_pipeline_input_panda=pd.concat(pandas_list,axis='index',ignore_index=True)
    total_pipeline_input_panda.to_pickle('../results/'+str(min_fold_change)+'/step_0_b_shape_aws_pull_to_pipeline_input/overall_pipeline_input.bin')


if __name__ == "__main__":

    min_fold_change=sys.argv[1]
    named_or_all=sys.argv[2]
    #output_pickle_address='../results/'+str(min_fold_change)+'/step_0_c_complete_pipeline_input/binvestigate_species_transformed.bin'
    os.system('mkdir -p ../results/'+str(min_fold_change)+'/step_0_c_complete_pipeline_input/')
    os.system('touch ../results/'+str(min_fold_change)+'/step_0_c_complete_pipeline_input/dummy.txt')

    addtional_property_csv_address='../resources/pull_from_carrot/intermediates/bins_groups_inchi_from_gert_inchikey_group_harmonized.csv'
    additional_property_panda=pd.read_csv(addtional_property_csv_address,sep=',')##,index_col=0)
 
    pipeline_input_panda_directory='../results/'+str(min_fold_change)+'/step_0_b_shape_aws_pull_to_pipeline_input/'
    pipeline_output_directory='../results/'+str(min_fold_change)+'/step_0_c_complete_pipeline_input/'
    
    file_list=os.listdir(pipeline_input_panda_directory)
    file_list.remove('dummy.txt')
    
    for file_counter,temp_file in enumerate(file_list):
        print('we are on file number: '+str(file_counter))
        print(f'we are on file {temp_file}')
        temporary_input_panda=pd.read_pickle(pipeline_input_panda_directory+temp_file)

        temporary_input_panda=temporary_input_panda.loc[
            ~((temporary_input_panda['name'].str[0].str.lower()=='z') & (temporary_input_panda['name'].str[1]==' ')),
            :
        ]
 
        input_addtional_properties(temporary_input_panda,additional_property_panda)
        temporary_file_integer=re.findall(r'\d+', temp_file)[0]
        temporary_input_panda.to_pickle(pipeline_output_directory+'pipeline_input_group_properties_added_'+str(temporary_file_integer)+'.bin',protocol=0)
