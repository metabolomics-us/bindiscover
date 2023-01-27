#parker wrote a number of views/mat-views on carrot that transform the data into something that we want to use for this analysis
#this script takes data from those views and puts it into local files for statistical analysis
#the vast majority of those files are just lists of normalized intensities with the species organ compound as file name
#however there are two files that give us properties from carrot also. the properties are things like "how many total samples for a species organ"

import os
import sys
from pprint import pprint
import json
import pandas as pd

from sqlalchemy import create_engine
from sqlalchemy import Table, String
from sqlalchemy.dialects import postgresql

my_server='restore-parker.czbqhgrlaqbf.us-west-2.rds.amazonaws.com:5432'
my_database='carrot'
my_dialect='postgresql'
my_driver='psycopg2'
my_username='postgres'
my_password='Fiehnlab2020'
my_connection=f'{my_dialect}+{my_driver}://{my_username}:{my_password}@{my_server}/{my_database}'
my_engine=create_engine(my_connection)#,echo=True)



def obtain_all_species_organs_compounds():
    all_species_organs_compounds_cursor=connection.execute('''
        select
            psdsvasl.species,
            psdsvasl.organ,
            psdsvasl.target_id,
            foo.target_type,
            foo.name
        from
            plb_short_data_soc_vs_average_stddev psdsvasl
        inner join (
            select distinct on 
                (target_id)
                target_id,
                "name",
                target_type
            from plb_short_intensites_normalized_by_sum_sample_fame psinbssf 
        ) foo
        on
        psdsvasl.target_id = foo.target_id
        '''
    )
    all_species_organs_compounds_result=json.dumps([dict(r) for r in all_species_organs_compounds_cursor])
    all_species_organs_compounds_panda=pd.read_json(
        all_species_organs_compounds_result,
        orient='records'
    )

    return all_species_organs_compounds_panda


def obtain_all_species_organ_sample_count_and_average_fames():
    species_organ_sample_count_and_average_fame_cursor=connection.execute('''
        select
            *
        from
            plb_short_data_so_vs_sample_count_and_average_fame
        '''
    )
    species_organ_sample_count_and_average_fame_result=json.dumps([dict(r) for r in species_organ_sample_count_and_average_fame_cursor])
    species_organ_sample_count_and_average_fame_panda=pd.read_json(
        species_organ_sample_count_and_average_fame_result,
        orient='records'
    )
    return species_organ_sample_count_and_average_fame_panda

def aquire_normalized_intensities_for_species_organ_compound(temp_species,temp_organ,temp_compound):
    all_samples_for_soc_cursor=connection.execute(
        f'''
        select
            normalized_intensity
        from
            plb_short_intensites_normalized_by_sum_sample_fame
        where
            species=\'{temp_species}\' AND
            organ=\'{temp_organ}\' AND
            target_id=\'{temp_compound}\'
        '''
    )

    all_samples_for_soc_result=json.dumps([dict(r) for r in all_samples_for_soc_cursor])
    all_samples_for_soc_panda=pd.read_json(
        all_samples_for_soc_result,
        orient='records'
    )

    return all_samples_for_soc_panda

def aquire_normalized_intensities_for_species_organ_compound_wrapper(temp_panda):
    for index,series in temp_panda.iterrows():
        if index%10==0:
            print(index)
        temp_results=aquire_normalized_intensities_for_species_organ_compound(
            series['species'],
            series['organ'],
            series['target_id']
        )

        temp_file_address='../results/'+str(min_fold_change)+'/step_0_a_pull_distributions_from_aws/soc_data/'+series['species']+'¬'+series['organ']+'¬'+str(series['target_id'])+'.bin'

        temp_results.to_pickle(temp_file_address)


if __name__=="__main__":
    

    min_fold_change=sys.argv[1]
    output_pickle_address='../results/'+str(min_fold_change)+'/step_0_a_pull_distributions_from_aws/binvestigate_species_transformed.bin'
    os.system('mkdir -p ../results/'+str(min_fold_change)+'/step_0_a_pull_distributions_from_aws/')
    os.system('touch ../results/'+str(min_fold_change)+'/step_0_a_pull_distributions_from_aws/dummy.txt')
    os.system('mkdir -p ../results/'+str(min_fold_change)+'/step_0_a_pull_distributions_from_aws/so_count_data')
    os.system('mkdir -p ../results/'+str(min_fold_change)+'/step_0_a_pull_distributions_from_aws/soc_data')


    connection=my_engine.connect()

    all_species_organs_compounds_panda=obtain_all_species_organs_compounds()
    all_species_organs_compounds_panda.to_pickle('../results/'+str(min_fold_change)+'/step_0_a_pull_distributions_from_aws/so_count_data/all_species_organs_compounds_panda.bin')
    print(all_species_organs_compounds_panda)

    species_organ_sample_count_and_average_fame_panda=obtain_all_species_organ_sample_count_and_average_fames()
    species_organ_sample_count_and_average_fame_panda.to_pickle('../results/'+str(min_fold_change)+'/step_0_a_pull_distributions_from_aws/so_count_data/species_organ_sample_count_and_average_fame.bin')
    print(species_organ_sample_count_and_average_fame_panda)

    aquire_normalized_intensities_for_species_organ_compound_wrapper(all_species_organs_compounds_panda)

    connection.close()
    my_engine.dispose()