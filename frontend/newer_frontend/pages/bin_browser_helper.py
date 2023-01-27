from enum import unique
import pandas as pd

def generate_bin_dropdown_options():
    final_curations=pd.read_pickle('../newer_datasets/compound_translation_for_all_components.bin')
    final_curations.loc[final_curations.bin_type=='known','english_name']='Known: '+final_curations.loc[final_curations.bin_type=='known']['english_name'].astype(str)
    final_curations.drop(['bin_type','identifier','integer_representation'],axis='columns',inplace=True)
    final_curations.rename(columns={'compound_identifier':'value','english_name':'label'},inplace=True)
    return final_curations.to_dict(
        'records'
    )

def generate_compound_classes():
    compound_classes=pd.read_csv('../newer_datasets/classes_curated_map.txt',sep='\t')
    compound_classes.drop_duplicates(subset='InChIKey',inplace=True)
    compound_classes.set_index('InChIKey',drop=True,inplace=True)
    return compound_classes
