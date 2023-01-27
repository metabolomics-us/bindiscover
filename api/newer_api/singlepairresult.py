from leafquery import LeafQuery
from sqlalchemy import create_engine
from flask import Flask,request
from flask_restful import Api, Resource, reqparse
import json
import pandas as pd
import time

triplet_translation_panda=pd.read_pickle('../newer_datasets/triplet_translation_panda.bin')
compound_translation_panda=pd.read_pickle('../newer_datasets/compound_translation_for_all_components.bin')

class SinglePairResult:

    def __init__(self,triplet_from,triplet_to,bin_type,connection):

        #extremely annoyingly, the data were put into the database with the triplet order O,S,D
        #yet we receive them with order S,O,D
        #so we need to check alphabetical based on a rearrangement
        triplet_from_osd=triplet_from.split(' - ')[1]+' - '+triplet_from.split(' - ')[0]+' - '+triplet_from.split(' - ')[2]
        triplet_to_osd=triplet_to.split(' - ')[1]+' - '+triplet_to.split(' - ')[0]+' - '+triplet_to.split(' - ')[2]

        #if triplet from comes before triplet to in the alphabet
        #then we are going for something in the "top right corder of the fold matrix"
        #and so we can fetch directly
        if triplet_from_osd<triplet_to_osd:
            need_to_swap_fold_sign=False
        elif triplet_from_osd>triplet_to_osd:
            need_to_swap_fold_sign=True

        triplet_from_int=triplet_translation_panda.at[triplet_from,'integer_representation']
        triplet_to_int=triplet_translation_panda.at[triplet_to,'integer_representation']

        if need_to_swap_fold_sign==False:
            
            temp_LeafQuery=LeafQuery(
                triplet_from_int,triplet_to_int,bin_type
            )

            temp_cursor=connection.execute(temp_LeafQuery.query)

            if (temp_cursor.rowcount <= 0):
                connection.close()
                #https://stackoverflow.com/questions/8645250/how-to-close-sqlalchemy-connection-in-mysql
                #my_engine.dispose()
                print('row count of final result cursor less than 1')
                return 'fail'
            else:
                temp_result=json.dumps([dict(r) for r in temp_cursor])

        elif need_to_swap_fold_sign==True:
            
            temp_LeafQuery=LeafQuery(
                #note the swapped order compared to ==False
                triplet_to_int,triplet_from_int,bin_type
            )

            temp_cursor=connection.execute(temp_LeafQuery.query)

            if (temp_cursor.rowcount <= 0):
                connection.close()
                #https://stackoverflow.com/questions/8645250/how-to-close-sqlalchemy-connection-in-mysql
                #my_engine.dispose()
                print('row count of final result cursor less than 1')
                return 'fail'
            else:
                #read into panda form, change sign of folds, then return to records form
                temp_panda=pd.read_json(
                    json.dumps([dict(r) for r in temp_cursor]),
                    orient='records'
                )

                temp_panda['fold_change_average']=-1*temp_panda['fold_change_average']
                temp_result=temp_panda.to_json(orient='records')

        #if nothing has failed so far, we translate the compouns column then return
        temp_panda=pd.read_json(temp_result,orient='records')

        identifier_dict=dict(zip(compound_translation_panda.integer_representation,compound_translation_panda.identifier))
        english_name_dict=dict(zip(compound_translation_panda.integer_representation,compound_translation_panda.english_name))
        bin_type_dict=dict(zip(compound_translation_panda.integer_representation,compound_translation_panda.bin_type))


        temp_panda['identifier']=temp_panda['compound_id'].map(identifier_dict.get)
        temp_panda['english_name']=temp_panda['compound_id'].map(english_name_dict.get)
        temp_panda['bin_type_dict']=temp_panda['compound_id'].map(bin_type_dict.get)

        self.result_panda=temp_panda
        
        

