import sys
import pandas
import os
import re

#a note - the classyfire tool is weird - for certain compounds, in will repeat the lowest class a few times
#sometimes with blanks between levels. example - http://classyfire.wishartlab.com/entities/WJNGQIYEQLPJMN-IOSLPCCCSA-N
#where there is no subclass, but a direct parent is reported
#for the purposes of our analysis, this is OK. that is because we only take the most specific thing and map to the classyfire taxonomy
#we are not concerned with weird spaces in the tool's reporting system

#the general logic for this script existing in isolation compared to the bins
#is the fact that the count dropping method is slow
#so if we make changes to the classes we want to be able to update quickly

#the general logic for this is like other transforms - check a mapping.txt and replace

#if the bin did not have an inchikey identity, then the classes are all "pre_curation_file"

def print_bin_information_for_classes_curated(temp_panda):
    '''
    this prints the id and inchikey_curated for creation of the curated_class transform mapping.tsv
    '''
    pandas.options.display.max_rows=10000
    #plb update 2-6-2022
    #we now want the pipeline to handle unknown compounds, so we print things where the inchikey curated is "junk"
    #plb update 7-4-22
    ##just for reference, we want the first half of the pipeline to handle unknowns. unknowns will not be "compared" in a volcano plot way tho
    ##no reason to visit this function
    #print(temp_panda.loc[(temp_panda['inchikey_curated']!='@@@@@@@') & (temp_panda['inchikey_curated']!='pre_curation_file') ][['id','inchikey_curated']])
    print(temp_panda[['id','inchikey_curated']])
    pandas.reset_option('display.max_rows')

def update_curated_classes_from_mapping(temp_panda, temp_class_mapping_address):
    '''
    here we add columns for classes results that come from inchikeys

    we do the same matching style based on 'id' that we do for the inchikeys
    '''
    temp_panda['class_from_curation_not_ML']=True

    class_mapping_panda=pandas.read_csv(temp_class_mapping_address,sep='\t')

    temp_panda=temp_panda.merge(
        right=class_mapping_panda,
        how='left',
        left_on='inchikey',
        right_on='InChIKey'
    )

    temp_panda.drop(
        'Status',
        axis='columns',
        inplace=True
    )

    return temp_panda


if __name__ == "__main__":

    #if snakemake in globals():
    #min_fold_change=snakemake.params.min_fold_change
    min_fold_change=sys.argv[1]
    class_mapping_address='../resources/species_organ_maps/classes_curated_map.txt'
    os.system('mkdir -p ../results/'+str(min_fold_change)+'/step_4_classes_transformed/')
    os.system('touch ../results/'+str(min_fold_change)+'/step_4_classes_transformed/dummy.txt')
    
    pipeline_input_panda_directory='../results/'+str(min_fold_change)+'/step_3_bins_transformed/'
    pipeline_output_directory='../results/'+str(min_fold_change)+'/step_4_classes_transformed/'
    
    
    file_list=os.listdir(pipeline_input_panda_directory)
    file_list.remove('dummy.txt')
    
    for temp_file in file_list:


        initial_panda=pandas.read_pickle(pipeline_input_panda_directory+temp_file)
        #the classes are printed
        #and put in the proper file if necessary
        #plb 7-4-22 no reason to visit this

        #fill the class column 
        initial_panda=update_curated_classes_from_mapping(initial_panda, class_mapping_address)

        ###############################################
        ###############################################
        #later, we may add a class from a ML algorithm#
        #in this way, each class could be added to the compound hierarchy#
        #this was the original intent. not any more.
        ###############################################
        ###############################################
        temporary_file_integer=re.findall(r'\d+', temp_file)[0]
        initial_panda.to_pickle(pipeline_output_directory+'binvestigate_classes_transformed_'+str(temporary_file_integer)+'.bin',protocol=0)
