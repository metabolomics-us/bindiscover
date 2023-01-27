import pandas as pd
import sys
import os

if __name__ == "__main__":

    
    min_fold_change=sys.argv[1]

    os.system('mkdir -p ../results/'+str(min_fold_change)+'/step_5_b_make_non_ratio_table/')
    os.system('touch ../results/'+str(min_fold_change)+'/step_5_b_make_non_ratio_table/dummy.txt')

    pipeline_input_panda_directory='../results/'+str(min_fold_change)+'/step_5_panda_cleaned/'
    
    output_panda_address='../results/'+str(min_fold_change)+'/step_5_b_make_non_ratio_table/non_ratio_table.bin'
    output_dict={
        'bin':[],
        'compound':[],
        'species':[],
        'organ':[],
        'disease':[],
        'intensity_average':[],
        'intensity_median':[],
        'percent_present':[]
    }

    file_list=os.listdir(pipeline_input_panda_directory)
    file_list.remove('dummy.txt')
    for temp_file in file_list:
        temp=pd.read_pickle(pipeline_input_panda_directory+temp_file)
 

        for index,series in temp.iterrows():
            output_dict['bin']+=[series['id'] for i in range(len(series['species']))]
            output_dict['compound']+=[series['name'] for i in range(len(series['species']))]
            output_dict['species']+=series['species']
            output_dict['organ']+=series['organ']
            output_dict['disease']+=series['special_property_list']
            output_dict['intensity_average']+=series['total_intensity']
            output_dict['intensity_median']+=series['median_intensity']
            output_dict['percent_present']+=series['percent_present']

    temp_2=pd.DataFrame.from_dict(output_dict)
    
    #added because oliver says that very small decimals will freak people out
    temp_2['intensity_average']=1e9*temp_2['intensity_average']
    temp_2['intensity_median']=1e9*temp_2['intensity_median']
    
    #for some reason bin 13110 was appearing twice
    #hotfix    
    temp_2.drop_duplicates(subset=['bin','species','organ','disease'],keep='first',inplace=True,ignore_index=True)

    temp_2.to_pickle(output_panda_address)

    temp_2['combined_sod']=temp_2.species+' - '+temp_2.organ+' - '+temp_2.disease
    temp_2.combined_sod.value_counts().to_pickle(
        '../results/'+str(min_fold_change)+'/step_5_b_make_non_ratio_table/unique_sod_combinations.bin'
    )


