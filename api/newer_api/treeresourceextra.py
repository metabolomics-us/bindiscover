from sqlalchemy import create_engine
from flask import Flask,request 
from flask_restful import Api, Resource, reqparse
import json
import pandas as pd
import numpy as np

my_server='binvestigate-parker.czbqhgrlaqbf.us-west-2.rds.amazonaws.com'
my_database='postgres'
my_dialect='postgresql'
my_driver='psycopg2'
my_username='postgres'
my_password='S7kB93DIT46$'
my_port='5432'
my_connection=f'{my_dialect}+{my_driver}://{my_username}:{my_password}@{my_server}/{my_database}'
my_engine=create_engine(my_connection)#,echo=True)


class TreeResource(Resource):

    def post(self):
        '''
        what we really want to do here is cycle through the query for each triplet requested
        '''
        metadata_triplets=request.json['metadata_triplets']
        datatype=request.json['data_type']
        bintype=request.json['bin_type']

        connection=my_engine.connect()
        
        temp_panda_list=list()

        for temp_triplet in metadata_triplets:
            temp_species=temp_triplet.split(' - ')[0]
            temp_organ=temp_triplet.split(' - ')[1]
            temp_disease=temp_triplet.split(' - ')[2]
        
            if datatype=='percent_present' and bintype=='knowns':
                temp_cursor=connection.execute(
                    f'''
                    select bin,percent_present from non_ratio_table_phylo_knowns where (species='{temp_species}') and (organ='{temp_organ}') and (disease='{temp_disease}')
                    '''
                )

                input_records_panda=pd.read_json(json.dumps([dict(r) for r in temp_cursor]), orient="records")
                input_records_panda['triplet']=temp_triplet
                temp_panda_list.append(input_records_panda)


            elif datatype=='percent_present' and bintype=='all':
                temp_cursor=connection.execute(
                    f'''
                    select bin,percent_present from non_ratio_table where (species='{temp_species}') and (organ='{temp_organ}') and (disease='{temp_disease}')
                    '''
                )

                input_records_panda=pd.read_json(json.dumps([dict(r) for r in temp_cursor]), orient="records")
                input_records_panda['triplet']=temp_triplet
                temp_panda_list.append(input_records_panda)

        output_panda=pd.concat(temp_panda_list,axis='index',ignore_index=True)
        output_panda=output_panda.pivot(
            index='triplet',
            columns='bin',
            values='percent_present'
        )
        output_panda=output_panda.reindex(
            metadata_triplets,

        )
        #what we need to do now is to pivot the triplet column into indexes and make each of the cmpounds into a column
        #note to self: need to make query significantly faster. can do this by adjusting the underlying database. perhaps make a new table
       
        if (temp_cursor.rowcount <= 0):
            connection.close()
            #https://stackoverflow.com/questions/8645250/how-to-close-sqlalchemy-connection-in-mysql
            my_engine.dispose()
            print('row count of final result cursor less than 1')
            return 'fail'
        else:
            #temp_result=json.dumps([dict(r) for r in temp_cursor])
            connection.close()
            #https://stackoverflow.com/questions/8645250/how-to-close-sqlalchemy-connection-in-mysql
            my_engine.dispose()
            #print(temp_result)

            return output_panda.to_json(orient='records')