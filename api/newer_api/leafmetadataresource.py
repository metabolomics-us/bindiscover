from flask import Flask,request
from flask_restful import Api, Resource, reqparse
import json

import pandas as pd 

triplet_translation_panda=pd.read_pickle('../newer_datasets/triplet_translation_panda.bin')

class LeafMetadataResource(Resource):

    def post(self):
        '''
        should check if each species,organ,disease headnode in headnodes_to_triplets
        should verify types...
        should verify valid ranges (like p value between 0 and 1).....
        '''

        from_triplets=request.json['triplet_from']
        to_triplets=request.json['triplet_to']

        output_dict={
            'from_or_to':[],
            'triplet_id':[],
            'sample_count':[]
        }

        for triplet in from_triplets:
            output_dict['from_or_to'].append('from')
            output_dict['triplet_id'].append(triplet_translation_panda.at[triplet,'triplet_identifier_string'])
            output_dict['sample_count'].append(triplet_translation_panda.at[triplet,'count'])
        for triplet in to_triplets:
            output_dict['from_or_to'].append('to')
            output_dict['triplet_id'].append(triplet_translation_panda.at[triplet,'triplet_identifier_string'])
            output_dict['sample_count'].append(triplet_translation_panda.at[triplet,'count'])        


        return pd.DataFrame.from_dict(output_dict).to_json(orient='records')
