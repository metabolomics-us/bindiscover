import sys
import os
import pandas as pd
#import networkx as nx


if __name__ == "__main__":

    #min_fold_change=0
    min_fold_change=sys.argv[1]
    os.system('mkdir -p ../results/'+str(min_fold_change)+'/step_9_b_generate_bin_spectral_panda/')
    os.system('touch ../results/'+str(min_fold_change)+'/step_9_b_generate_bin_spectral_panda/dummy.txt')

    binvestigate_spectral_data_panda_address='../../archives/binvestigate_pull/post_sync/total_binvestigate_panda.bin'
    binvestigate_spectral_panda=pd.read_pickle(binvestigate_spectral_data_panda_address)

    compound_translation_panda_address='../results/'+str(min_fold_change)+'/step_9_generate_extras_for_db_and_api/compound_translation_panda.bin'
    compound_translation_panda=pd.read_pickle(compound_translation_panda_address)

    binvestigate_spectral_panda_output_address='../results/'+str(min_fold_change)+'/step_9_b_generate_bin_spectral_panda/bin_spectral_panda.bin'

    #remove rows not in translation panda (compound identifier
    binvestigate_spectral_panda=binvestigate_spectral_panda.loc[binvestigate_spectral_panda.id.isin(compound_translation_panda.compound_identifier)]
    
    #reshape spectra_mz and spectra_intensity into array (string)
    spectrum_list=list()
    for row,series in binvestigate_spectral_panda.iterrows():
        temp_mzs=series['spectra_mz'][0]
        temp_intensities=[float(x) for x in series['spectra_intensity'][0]]
        temp_max_intensity=max(temp_intensities)
        temp_intensities_normalized=[x/temp_max_intensity for x in temp_intensities]
        normalized_spectrum=[str(temp_mzs[x])+':'+str(temp_intensities_normalized[x]) for x in range(len(temp_mzs))]
        uploaded_string=' '.join(normalized_spectrum)
        #[uploaded_string:=uploaded_string+str(temp_mzs[x])+':'+str(temp_intensities_normalized[x])+' ' for x in range(len(temp_mzs))]
        spectrum_list.append(uploaded_string)
    ##normalize spectrum intensity
    binvestigate_spectral_panda['spectrum']=spectrum_list

    #keep only a handful of columns
    binvestigate_spectral_panda.drop(['spectra_mz','spectra_intensity','name','sample','_id','species','organ','count','intensity','group','inchikey'],axis='columns',inplace=True)
    #get enlgish name from compound translation panda
    #get inchikey from compound translation panda
    compound_translation_panda.drop(['integer_representation','bin_type','identifier'],axis='columns',inplace=True)
    binvestigate_spectral_panda=binvestigate_spectral_panda.merge(
        right=compound_translation_panda,
        how='left',
        left_on='id',
        right_on='compound_identifier'
    )
    
    binvestigate_spectral_panda.drop('id',axis='columns',inplace=True)
    binvestigate_spectral_panda.compound_identifier=binvestigate_spectral_panda.compound_identifier.astype(int)
    
    binvestigate_spectral_panda.to_pickle(binvestigate_spectral_panda_output_address)