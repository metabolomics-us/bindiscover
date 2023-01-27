from operator import index
from flask import Flask,request 
from flask_restful import Api, Resource, reqparse
import json
from itertools import product
from metadatapanda import MetadataPanda

import networkx as nx
import pandas as pd


class HGDAMetadataResource(Resource):

    def post(self):

        #print('in resource')
        species_from=request.json['from_species']
        organ_from=request.json['from_organ']
        disease_from=request.json['from_disease']
        species_to=request.json['to_species']
        organ_to=request.json['to_organ']
        disease_to=request.json['to_disease']      

        #check if they sent the same things
        if species_from==species_to and organ_from==organ_to and disease_from==disease_to:
            return 'fail - positions vs self'

        my_MetadataPanda=MetadataPanda(species_from,organ_from,disease_from,species_to,organ_to,disease_to)

        return my_MetadataPanda.result_panda.to_json(orient='records')