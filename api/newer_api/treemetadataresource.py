#from operator import index
from flask import Flask,request 
from flask_restful import Api, Resource, reqparse
import json
from itertools import product
#from metadatapanda import MetadataPanda

import networkx as nx
import pandas as pd


index_panda=pd.read_pickle('../newer_datasets/index_panda.bin')
index_panda['species']=index_panda['species'].astype(str)
species_networkx=nx.read_gpickle('../newer_datasets/species_networkx.bin')
organ_networkx=nx.read_gpickle('../newer_datasets/organ_networkx.bin')
disease_networkx=nx.read_gpickle('../newer_datasets/disease_networkx.bin')
triplet_translation_panda=pd.read_pickle('../newer_datasets/triplet_translation_panda.bin')

species_to_disallow_for_tanglegram={'arabidopsis', 'coffea', 'brassica', 'chlamydomonas', 'lactobacillales', 'gossypium', 'synechococcus', 'chlorophyta', 'citrus', 'carabidae', 'clostridium', 'phaseolus', 'bacteria'}


class TreeMetadataResource(Resource):

    def post(self):

        species=request.json['species']
        organ=request.json['organ']
        disease=request.json['disease']
        
        
        if (species is not None) & (organ is not None) & (disease is not None):
            species_mapped_to_set=self.determine_mapped_to_nodes(index_panda,species,species_networkx,'species')
            try:
                species_mapped_to_set.discard('ralstonia eutropha')
            except:
                pass
            try:
                species_mapped_to_set.discard('chromobacterium')
            except:
                pass
            try:
                species_mapped_to_set.discard('helicobacter pylori')
            except:
                pass
            species_mapped_to_set={element for element in species_mapped_to_set if (element not in species_to_disallow_for_tanglegram)}
            organ_mapped_to_set=self.determine_mapped_to_nodes(index_panda,organ,organ_networkx,'organ')
            disease_mapped_to_set=self.determine_mapped_to_nodes(index_panda,disease,disease_networkx,'disease')
            
            final_set_of_triplets=self.determine_valid_triplets(triplet_translation_panda,species_mapped_to_set,organ_mapped_to_set,disease_mapped_to_set)
            self.onto_triplets_panda=self.generate_metadata_results(triplet_translation_panda,final_set_of_triplets)
        else:
            self.onto_triplets_panda=pd.DataFrame.from_dict({'triplet_id':[],'sample_count':[]})
        #################################

        ##deal with individually named things
        output_dict={
            'triplet_id':[],
            'sample_count':[]
        }
        if request.json['species_organ_disease'] is not None:
            species_organ_disease=request.json['species_organ_disease']
            triplet_translation_panda.set_index('triplet_identifier_string',drop=False,inplace=True)

            for triplet in species_organ_disease:
                output_dict['triplet_id'].append(triplet)
                output_dict['sample_count'].append(triplet_translation_panda.at[triplet,'count'])

        self.individually_named_triplets_panda=pd.DataFrame.from_dict(output_dict)

        output_panda=pd.concat([self.onto_triplets_panda,self.individually_named_triplets_panda],axis='index',ignore_index=True)
        output_panda.drop_duplicates(inplace=True)
        return output_panda.to_json(orient='records')



    def determine_mapped_to_nodes(self,index_panda,temp_node,temp_nx,temp_type):
        node_identifier_set=nx.algorithms.dag.descendants(temp_nx,temp_node)
        node_identifier_set.add(temp_node)
        if temp_type=='species':
            to_english_translation_dict=dict((zip(index_panda.species,index_panda.species_english)))
        elif temp_type=='organ':
            to_english_translation_dict=dict((zip(index_panda.organ,index_panda.organ_english)))
        elif temp_type=='disease':
            to_english_translation_dict=dict((zip(index_panda.disease,index_panda.disease_english)))

        node_identifier_set=set(map(lambda x:to_english_translation_dict.get(x),node_identifier_set))
        node_identifier_set.discard(None)
        return node_identifier_set
        
    def determine_valid_triplets(self,triplet_translation_panda,species_set,organ_set,disease_set):
        all_combinations_requested=set(product(organ_set,species_set,disease_set))
        mapped_to_combinations=set(triplet_translation_panda.triplet_identifier_tuple.unique())
        
        valid_queried_combinations=mapped_to_combinations.intersection(all_combinations_requested)
        return valid_queried_combinations
        

    def generate_metadata_results(self,triplet_translation_panda,triplets):

        triplet_translation_panda.set_index('triplet_identifier_tuple',drop=False,inplace=True)

        output_dict={
            'triplet_id':[],
            'sample_count':[]
        }
        for triplet in triplets:
            output_dict['triplet_id'].append(triplet_translation_panda.at[triplet,'triplet_identifier_string'])
            output_dict['sample_count'].append(triplet_translation_panda.at[triplet,'count'])

        return pd.DataFrame.from_dict(output_dict)