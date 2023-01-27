from singlepairresult import SinglePairResult
import pandas as pd
from itertools import product
import numpy as np
import time



class GroupQueryResult:

    def __init__(self,another_temp_panda,bin_type,connection):
        triplet_pair_list=self.form_triplet_pairs(another_temp_panda)
        # print(triplet_pair_list)
        for i, pair in enumerate(triplet_pair_list):
            if i==0:
                # get the result
                temp_SinglePairResult=SinglePairResult(pair[0],pair[1],bin_type,connection)                
                total_result_panda=temp_SinglePairResult.result_panda.copy()
            elif i>=1:
                temp_SinglePairResult=SinglePairResult(pair[0],pair[1],bin_type,connection)
                total_result_panda=pd.concat([total_result_panda,temp_SinglePairResult.result_panda],axis='index',ignore_index=True)
            #print(f'the total time to get or concatenate a single result is {end-start}')

            if i%5==0:
                #we do this again after all loopin gto clear up any remainder. we do it during the loopin gso that the memeory doesnt get insane
                total_result_panda_grouped=total_result_panda.groupby(
                    by=['compound_id', 'identifier', 'english_name', 'bin_type_dict'], axis='index',as_index=False,sort=False
                )

                rolling_result_significance=total_result_panda_grouped[['significance_welch']].max()
                rolling_result_fold_avg=self.one_df_transform_fold(total_result_panda_grouped['fold_change_average'],'fold_change_average')
                
                total_result_panda=rolling_result_significance.copy()
                total_result_panda['fold_change_average']=rolling_result_fold_avg

                self.result_panda=total_result_panda.copy()

        start=time.time()
        #we do this again after all loopin gto clear up any remainder. we do it during the loopin gso that the memeory doesnt get insane
        total_result_panda_grouped=total_result_panda.groupby(
            by=['compound_id', 'identifier', 'english_name', 'bin_type_dict'], axis='index',as_index=False,sort=False
        )

        #total_result_panda
        rolling_result_significance=total_result_panda_grouped[['significance_welch']].max()

        rolling_result_fold_avg=self.one_df_transform_fold(total_result_panda_grouped['fold_change_average'],'fold_change_average')

        total_result_panda=rolling_result_significance.copy()
        total_result_panda['fold_change_average']=rolling_result_fold_avg

        self.result_panda=total_result_panda.copy()
        

    def one_df_transform_fold(self,temp_groupby,temp_fold_type):
        '''
        given an numpy array, chooses what the aggregate value is
        '''

        my_groupby_min=temp_groupby.min()

        my_groupby_max=temp_groupby.max()       

        intermediate_min_df=my_groupby_min[temp_fold_type].where(
            np.sign(my_groupby_min[temp_fold_type])>0,
            other=0
        )
        intermediate_max_df=my_groupby_max[temp_fold_type].where(
            np.sign(my_groupby_max[temp_fold_type])<0,
            other=0
        )
        
        temp_output=intermediate_min_df.where(
            intermediate_min_df !=0,
            other=intermediate_max_df
        )
        return temp_output

    def form_triplet_pairs(self,temp_panda):
        from_trips=set(temp_panda.loc[temp_panda.from_or_to=='from']['triplet_id'].unique())
        to_trips=set(temp_panda.loc[temp_panda.from_or_to=='to']['triplet_id'].unique())
        return set(product(from_trips,to_trips))
