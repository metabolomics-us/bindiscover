import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import io
import base64
from upsetplot import UpSet

#############Load pandas for data selection options ##########
def get_unique_sod_combinations():
    unique_sod_combinations_address = "../newer_datasets/unique_sod_combinations_common_and_database_names.bin"
    unique_sod_combinations_panda = pd.read_pickle(unique_sod_combinations_address)
    unique_sod_combinations_dict=dict(zip(unique_sod_combinations_panda['triplet_with_common_name'],unique_sod_combinations_panda['triplet']))
    return unique_sod_combinations_dict
##############################################################

########################upset###################################
def create_upset(temp_panda):
    '''
    Creates the upsetplot image.
    These matplotlib figures are saved in file and shown as images
    '''
    temp_MultiIndex=temp_panda.isnull()==False
    temp_panda.index=temp_MultiIndex
    temp_panda.index=pd.MultiIndex.from_tuples(temp_panda.index)
    temp_panda.index.set_names(names=temp_panda.columns,inplace=True)

    fig = plt.figure(figsize=(5, 5),dpi=200)
    my_UpSet = UpSet(temp_panda,subset_size='count',min_degree=1)
    #my_UpSet._plot_stacked_bars(title="hello")
    upset_subplot_dict=UpSet.plot(my_UpSet,fig=fig)
    upset_subplot_dict['intersections'].set_ylabel('Compound Count')
    buf = io.BytesIO() # in-memory files

    plt.savefig(buf, format = "png") # save to the above file object
    plt.close('all')
    data = base64.b64encode(buf.getbuffer()).decode("utf8") # encode to html elements
    plotly_fig="data:image/png;base64,{}".format(data)
    buf.close()
    return plotly_fig
####################################################################