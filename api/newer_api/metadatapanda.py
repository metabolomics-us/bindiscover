import pandas as pd
import networkx as nx
from itertools import product
from pprint import pprint

index_panda=pd.read_pickle('../newer_datasets/index_panda.bin')
index_panda['species']=index_panda['species'].astype(str)
species_networkx=nx.read_gpickle('../newer_datasets/species_networkx.bin')
organ_networkx=nx.read_gpickle('../newer_datasets/organ_networkx.bin')
disease_networkx=nx.read_gpickle('../newer_datasets/disease_networkx.bin')
triplet_translation_panda=pd.read_pickle('../newer_datasets/triplet_translation_panda.bin')


class MetadataPanda:

    def __init__(self,species_from,organ_from,disease_from,species_to,organ_to,disease_to):
        #determine parent status
        
        species_parent_status=self.determine_if_from_or_to_are_parents(species_from,species_to,species_networkx)
        organ_parent_status=self.determine_if_from_or_to_are_parents(organ_from,organ_to,organ_networkx)
        disease_parent_status=self.determine_if_from_or_to_are_parents(disease_from,disease_to,disease_networkx)

        #determine full mapped-to nodes for each from,to
        from_species_mapped_to_set=self.determine_mapped_to_nodes(index_panda,species_from,species_networkx,'species')
        to_species_mapped_to_set=self.determine_mapped_to_nodes(index_panda,species_to,species_networkx,'species')
        from_organ_mapped_to_set=self.determine_mapped_to_nodes(index_panda,organ_from,organ_networkx,'organ')
        to_organ_mapped_to_set=self.determine_mapped_to_nodes(index_panda,organ_to,organ_networkx,'organ')
        from_disease_mapped_to_set=self.determine_mapped_to_nodes(index_panda,disease_from,disease_networkx,'disease')
        to_disease_mapped_to_set=self.determine_mapped_to_nodes(index_panda,disease_to,disease_networkx,'disease')

        from_species_mapped_to_set,to_species_mapped_to_set=self.reduce_set_via_parent(from_species_mapped_to_set,to_species_mapped_to_set,species_parent_status)

        from_organ_mapped_to_set,to_organ_mapped_to_set=self.reduce_set_via_parent(from_organ_mapped_to_set,to_organ_mapped_to_set,organ_parent_status)
        from_disease_mapped_to_set,to_disease_mapped_to_set=self.reduce_set_via_parent(from_disease_mapped_to_set,to_disease_mapped_to_set,disease_parent_status)

        final_set_of_from_triplets=self.determine_valid_triplets(triplet_translation_panda,from_species_mapped_to_set,from_organ_mapped_to_set,from_disease_mapped_to_set)
        final_set_of_to_triplets=self.determine_valid_triplets(triplet_translation_panda,to_species_mapped_to_set,to_organ_mapped_to_set,to_disease_mapped_to_set)

        self.result_panda=self.generate_metadata_results(triplet_translation_panda,final_set_of_from_triplets,final_set_of_to_triplets)

        

    def determine_if_from_or_to_are_parents(self,from_node,to_node,temp_nx):
        if from_node==to_node:
            return 'neither'
        elif to_node in nx.algorithms.dag.descendants(temp_nx,from_node):
            return 'from_parent_of_to'
        elif from_node in nx.algorithms.dag.descendants(temp_nx,to_node):
            return 'to_parent_of_from'
        else:
            return 'neither'

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
        if temp_type=='species':
            try:
                node_identifier_set.discard('ralstonia eutropha')
            except:
                pass
            try:
                node_identifier_set.discard('chromobacterium')
            except:
                pass
            try:
                node_identifier_set.discard('helicobacter pylori')
            except:
                pass
        return node_identifier_set
        
    def reduce_set_via_parent(self,from_set,to_set,parent_status):
        if parent_status=='neither':
            return from_set,to_set
        elif parent_status=='from_parent_of_to':
            from_reduced_set=from_set.difference(from_set.intersection(to_set))
            return from_reduced_set,to_set
        elif parent_status=='to_parent_of_from':
            to_reduced_set=to_set.difference(to_set.intersection(from_set))
            return from_set,to_reduced_set

    def determine_valid_triplets(self,triplet_translation_panda,species_set,organ_set,disease_set):
        all_combinations_requested=set(product(organ_set,species_set,disease_set))
        mapped_to_combinations=set(triplet_translation_panda.triplet_identifier_tuple.unique())
        
        valid_queried_combinations=mapped_to_combinations.intersection(all_combinations_requested)
        return valid_queried_combinations
        

    def generate_metadata_results(self,triplet_translation_panda,from_triplets,to_triplets):

        triplet_translation_panda.set_index('triplet_identifier_tuple',drop=False,inplace=True)

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

        return pd.DataFrame.from_dict(output_dict)