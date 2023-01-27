import networkx as nx
import re
import pandas as pd


######### HELPER FUNCTIONS ################
def create_compound_selection_labels(final_curations_address):
    '''
    reads in labels and transforms them by last-minute curations
    '''
    # return compound_dropdown_options
    final_curations=pd.read_pickle(final_curations_address)
    final_curations.loc[final_curations.bin_type=='known','english_name']='Known: '+final_curations.loc[final_curations.bin_type=='known']['english_name'].astype(str)
    final_curations.drop(['bin_type','identifier','integer_representation'],axis='columns',inplace=True)
    final_curations.rename(columns={'compound_identifier':'value','english_name':'label'},inplace=True)
    
    #compound_dropdown_options=
    return final_curations.to_dict(
        'records'
    )


def coerce_full_panda(df,value_column,column_list):
    '''
    the plotly express sunburst was not appropriate for our goals, so we coerce the data to the form
    requested by the go object
    '''
    pandas_list=list()
    for i in range(len(column_list),0,-1):
        pandas_list.append(
            pd.DataFrame(
                data={
                    'count':df.groupby(by=column_list[0:i]).size().to_list(),
                    'sum':df.groupby(by=column_list[0:i])[value_column].sum().to_list(),
                    'parent':['/'.join(group[0][:i-1]) for group in df.groupby(by=column_list[0:i])],
                    'id':['/'.join(group[0][:i]) for group in df.groupby(by=column_list[0:i])],
                    'name':df.groupby(by=column_list[0:i])[column_list[i-1]].unique().map(lambda x: x[0]).values
                }
            )
        )
    tree_panda=pd.concat(pandas_list,axis='index')
    tree_panda.reset_index(inplace=True,drop=True)
    tree_panda.at[len(tree_panda.index)-1,'id']='binvestigate'
    tree_panda['average']=tree_panda['sum']/tree_panda['count']

    #there is a known bug in the way that branch totals works
    #https://community.plotly.com/t/plotly-sunburst-returning-empty-chart-with-branchvalues-total/26582/8
    #no matter what i tried, i could not get the branch total thing to work for me
    #so we use a hack workaround for now - everything except for the lowest levels is set to 0 for valeus
    #now we can use remainder and it should work as intended
    first_parent_index=len(df.index)
    tree_panda.loc[first_parent_index:,'sum']=0

    if value_column=='intensity_average':
        tree_panda=tree_panda.round(decimals=0)
        tree_panda['average']=tree_panda['average'].astype(int)

    elif value_column=='percent_present':
        tree_panda=tree_panda.round(decimals=3)
    
    return tree_panda