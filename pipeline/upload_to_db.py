import sys
import os
from sqlalchemy import create_engine
from sqlalchemy import Table, String
from sqlalchemy.dialects import postgresql
import pandas as pd
import psycopg2
import numpy as np
#https://stackoverflow.com/questions/50626058/psycopg2-cant-adapt-type-numpy-int64
from psycopg2.extensions import register_adapter, AsIs
#https://stackoverflow.com/questions/50626058/psycopg2-cant-adapt-type-numpy-int64
def adapt_numpy_int64(numpy_int64):
    return AsIs(numpy_int64)
import time

def choose_all_bins(directory_address):
    full_list=os.listdir(directory_address)
    return full_list

def prepare_one_bin_for_upload(temp_bin,mapping_dict,compound_mapping_dict):
    '''
    prepares one bin/class entry for the database
    we open all 4 fold/sig for one compund, stack them, typecast to dataframe, concat, change header names, change header order

    '''
    matrix_type_list=[
        '../results/'+str(min_fold_change)+'/step_8_perform_compound_hierarchical_analysis/all_matrices/fold_change_matrix_average/',
        '../results/'+str(min_fold_change)+'/step_8_perform_compound_hierarchical_analysis/all_matrices/fold_change_matrix_median/',
        '../results/'+str(min_fold_change)+'/step_8_perform_compound_hierarchical_analysis/all_matrices/signifigance_matrix_mannwhitney/',
        '../results/'+str(min_fold_change)+'/step_8_perform_compound_hierarchical_analysis/all_matrices/signifigance_matrix_welch/'
    ]
    
    pandas_list=[
        pd.read_pickle(temp_location+temp_bin) for temp_location in matrix_type_list
    ]
    
    for i,panda in enumerate(pandas_list):

        pandas_list[i].columns=pandas_list[i].columns.to_flat_index()
        pandas_list[i].index=pandas_list[i].index.to_flat_index()
        pandas_list[i].rename(columns=mapping_dict,inplace=True,errors='raise')
        pandas_list[i].rename(index=mapping_dict,inplace=True)
        
        pandas_list[i].values[np.tril_indices_from(pandas_list[i], 0)] = np.nan
        
        pandas_list[i]=pd.DataFrame(panda.stack())
        pandas_list[i].index.set_names(['triplet_from','triplet_to'],inplace=True)
        
    pandas_list[0].rename({0:'fold_change_average'},inplace=True,axis='columns')
    pandas_list[1].rename({0:'fold_change_median'},inplace=True,axis='columns')
    pandas_list[2].rename({0:'significance_mwu'},inplace=True,axis='columns')
    pandas_list[3].rename({0:'significance_welch'},inplace=True,axis='columns')
    total_panda=pd.concat(pandas_list,axis='columns')
    
    total_panda.reset_index(inplace=True)
    
    total_panda.insert(loc=0,column='compound_id',value=compound_mapping_dict[temp_bin[:-4]])
    total_panda=total_panda[[
        'compound_id','triplet_from', 
        'triplet_to','fold_change_average', 'fold_change_median',
        'significance_mwu', 'significance_welch'
    ]]
    
    return total_panda


def upload_fold_change_panda(temp_panda_for_upload,bin_iteration,connection):

    if bin_iteration==0:
        temp_panda_for_upload.to_sql(
            'differential_analysis',
            connection,
            index=False,
            dtype={
                'compound_id':postgresql.SMALLINT,
                'triplet_from':postgresql.SMALLINT,
                'triplet_to':postgresql.SMALLINT,
                'fold_change_average':postgresql.REAL,
                'fold_change_median':postgresql.REAL,
                'significance_mwu':postgresql.FLOAT,
                'significance_welch':postgresql.FLOAT,
            },
            if_exists='replace',
            method='multi',
            chunksize=90000
        )
  
    elif bin_iteration!=0:
        temp_panda_for_upload.to_sql(
            'differential_analysis',
            connection,
            index=False,
            dtype={
                'compound_id':postgresql.SMALLINT,
                'triplet_from':postgresql.SMALLINT,
                'triplet_to':postgresql.SMALLINT,
                'fold_change_average':postgresql.REAL,
                'fold_change_median':postgresql.REAL,
                'significance_mwu':postgresql.FLOAT,
                'significance_welch':postgresql.FLOAT,
            },
            if_exists='append',
            method='multi',
            chunksize=70000
        )

def upload_non_ratio_table(temp_panda,connection):
    temp_panda.to_sql(
        'non_ratio_table',
        connection,
        index=False,
        dtype={
            'bin':postgresql.INTEGER,
            'compound':postgresql.TEXT,
            'species':postgresql.TEXT,
            'organ':postgresql.TEXT,
            'disease':postgresql.TEXT,
            'intensity_average':postgresql.FLOAT,
            'intensity_median':postgresql.FLOAT,
            'percent_present':postgresql.FLOAT
        },
        if_exists='replace',
        method='multi',
        chunksize=90000
    )    
    
def upload_compound_translation_table(temp_panda,connection):
    temp_panda.to_sql(
        'compound_translation_table',
        connection,
        index=True,
        dtype={
            'compound_identifier':postgresql.TEXT,
            'integer_representation':postgresql.INTEGER
        },
        if_exists='replace',
        method='multi',
        chunksize=90000
    )    

def upload_triplet_translation_table(temp_panda,connection):
    temp_panda.to_sql(
        'compound_translation_table',
        connection,
        index=True,
        dtype={
            'triplet_identifier':postgresql.TEXT,
            'integer_representation':postgresql.INTEGER
        },
        if_exists='replace',
        method='multi',
        chunksize=90000
    )  

def upload_spectral_bin_table(spectral_bin_panda,connection):
    spectral_bin_panda.to_sql(
        'bin_table',
        connection,
        index=False,
        dtype={
            'retentionIndex':postgresql.FLOAT, 
            'kovats':postgresql.FLOAT, 
            'quantMass':postgresql.FLOAT, 
            'splash':postgresql.TEXT, 
            'purity':postgresql.FLOAT, 
            'uniqueMass':postgresql.FLOAT, 
            'spectrum':postgresql.TEXT, 
            'compound_identifier':postgresql.INTEGER, 
            'english_name':postgresql.TEXT
        },
        if_exists='replace',
        method='multi',
        chunksize=90000,
    )      


if __name__ == "__main__":

    min_fold_change=sys.argv[1]
    use_aws=(sys.argv[2])
    os.system('mkdir -p ../results/'+str(min_fold_change)+'/step_10_upload_to_db/')
    os.system('touch ../results/'+str(min_fold_change)+'/step_10_upload_to_db/dummy.txt')

    if use_aws=='False':
        my_server='localhost'
        my_database='binvestigate_second'
        my_dialect='postgresql'
        my_driver='psycopg2'
        my_username='rictuar'
        my_password='elaine123'
        my_port='5432'

    elif use_aws=='True':
        my_server='fold-result-database.czbab8f7pgfj.us-east-2.rds.amazonaws.com'
        my_database='foldresults'
        my_dialect='postgresql'
        my_driver='psycopg2'
        my_username='postgres'
        my_password='elaine123'
        my_port='5430'

    my_connection=f'{my_dialect}+{my_driver}://{my_username}:{my_password}@{my_server}:{my_port}/{my_database}'
    engine=create_engine(my_connection)#,echo=True)
    connection=engine.connect()

    #upload non-ratio table
    table_5_address='../results/'+str(min_fold_change)+'/step_5_b_make_non_ratio_table/non_ratio_table.bin'
    non_ratio_panda=pd.read_pickle(table_5_address)
    upload_non_ratio_table(non_ratio_panda,connection)
    start_time=time.time()
    #create our index
    connection.execute(
        f'''
        ALTER TABLE non_ratio_table ADD PRIMARY KEY (bin,species,organ,disease);
        '''
    )      
    end_time=time.time()
    print('time to create non ratio index: '+str(end_time-start_time))

    #upload the fold change matrices
    #get list of compounds and classes (listdir on one) (basically just bins)
    full_list=choose_all_bins('../results/'+str(min_fold_change)+'/step_8_perform_compound_hierarchical_analysis/all_matrices/fold_change_matrix_average')
    triplet_mapping_panda=pd.read_pickle('../results/'+str(min_fold_change)+'/step_9_generate_extras_for_db_and_api/triplet_translation_panda.bin')
    triplet_mapping_dict={
        triplet_mapping_panda.at[i,'triplet_identifier_tuple']:triplet_mapping_panda.at[i,'integer_representation'] for i in triplet_mapping_panda.index
    }
    compound_mapping_panda=pd.read_pickle('../results/'+str(min_fold_change)+'/step_9_generate_extras_for_db_and_api/compound_translation_panda.bin')
    compound_mapping_dict={
        compound_mapping_panda.at[i,'compound_identifier']:compound_mapping_panda.at[i,'integer_representation'] for i in compound_mapping_panda.index
    }  

    #220926 plb dont upload these, do the "math" in the api
    # #upload mapping pandas
    # upload_compound_translation_table(compound_mapping_panda,connection)
    # upload_triplet_translation_table(triplet_mapping_panda,connection)

    #for each bin, prepare each then upload each
    for i,temp_bin in enumerate(full_list):
        start_time=time.time()

        temp_panda_for_upload=prepare_one_bin_for_upload(temp_bin,triplet_mapping_dict,compound_mapping_dict)
        upload_fold_change_panda(temp_panda_for_upload,i,connection)

        end_time=time.time()
        print(temp_bin+', iteration '+str(i)+': '+str(end_time-start_time))


    start_time=time.time()
    #create our index
    connection.execute(
        f'''
        ALTER TABLE differential_analysis ADD PRIMARY KEY (compound_id, triplet_from, triplet_to);
        '''
    )      
    end_time=time.time()
    print('time to create complete differential analysis index: '+str(end_time-start_time))

    #upload non-ratio table
    table_9b_address='../results/'+str(min_fold_change)+'/step_9_b_generate_bin_spectral_panda/bin_spectral_panda.bin'
    spectral_bin_panda=pd.read_pickle(table_9b_address)
    upload_spectral_bin_table(spectral_bin_panda,connection)
    # start_time=time.time()
    # #create our index
    connection.execute(
        f'''
        ALTER TABLE bin_table ADD PRIMARY KEY (compound_identifier);
        '''
    )      
