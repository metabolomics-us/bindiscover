import plotly.graph_objects as go
import numpy as np
import dash_bootstrap_components as dbc
from dash import dcc, html, dash_table
import pandas as pd


def build_summary_statistics_div_children(x_hist,y_hist,summary_stats_table):
    output_children=list()

    output_children.append(
        dbc.Row(
            html.H2("Volcano Summary Statistics Triplets", className='text-center'),
        )
    )
    output_children.append(html.Br())

    output_children.append(
        # dbc.Row(
        #     ,
        # )
        dbc.Row(
            children=[
                dbc.Col(
                    children=[
                        dbc.Card(
                            children=[
                                html.H4('To update:',className='text-center'),
                                html.H4('• Lasso/Box select on Volcano Plot •',className='text-center'),
                                html.H4('• Filters in the DataTable •',className='text-center')
                            ],
                            color='#fff4e4',
                            style={
                                'textAlign':'center',
                                "box-shadow": "1px 2px 7px 0px red",
                                "border-radius": "10px"
                            }
                        ),
                        
                    ],
                    width=6
                ),
                dbc.Col(
                    children=[
                        html.H4("Compound Counts in Volcano Regions", className='text-center'),
                        summary_stats_table
                    ],
                    width=6
                ),
                # dbc.Col(width=6)
                # summary_stats_table        
            ]
        )
    )
    output_children.append(html.Br())
    output_children.append(
        dbc.Row(
            html.H2("Volcano Plot Compound Histograms", className='text-center'),
        )
    )
    output_children.append(
        dbc.Row(
            children=[
                dbc.Col(
                    children=[dcc.Graph(id='x_hist',figure=x_hist)],
                    width=6
                ),
                dbc.Col(
                    children=[dcc.Graph(id='y_hist',figure=y_hist)],
                    width=6
                ),
            ]
        )
    )
    output_children.append(html.Br())


    

    return output_children

def create_fold_hist(hist_panda):

    fold_bin_edges=list()


    if hist_panda['x'].min()<-1:
        if hist_panda['x'].max()>-1:
            upper_limit=-1
        else:
            upper_limit=hist_panda['x'].max()
        fold_bin_edges+=list(np.linspace(hist_panda['x'].min(),upper_limit,10))
    if ( len(hist_panda['x'].loc[(hist_panda['x']>-1) & (hist_panda['x']<1)].index) > 0 ) :
        if hist_panda['x'].max()>1:
            upper_limit=1
        else:
            upper_limit=hist_panda['x'].max()
        if hist_panda['x'].min()<-1:
            lower_limit=-1
        else:
            lower_limit=hist_panda['x'].min()


        fold_bin_edges+=list(np.linspace(lower_limit,upper_limit,9))
    if hist_panda['x'].max()>1:
        if hist_panda['x'].min()<1:
            lower_limit=1
        else:
            lower_limit=hist_panda['x'].min()


        fold_bin_edges+=list(np.linspace(lower_limit,hist_panda['x'].max(),10))

    #remove duplicates
    fold_bin_edges=sorted(list(set(fold_bin_edges)))
    widths=[(fold_bin_edges[i+1]-fold_bin_edges[i]) for i in range(len(fold_bin_edges)-1)]
    x_positions=[(0.5*fold_bin_edges[i]+0.5*fold_bin_edges[i+1]) for i in range(len(fold_bin_edges)-1)]


    marker_colors=[  "#FF0000" if ((temp_pos<-1) or (temp_pos>1)) else "#2884f4" for temp_pos in x_positions     ]

    ######################
    marker_colors_gfblp=[  "#2884f4" for temp_pos in x_positions     ]
    good_fold_but_large_p_hist=hist_panda.loc[
        (hist_panda['x'].abs()>1) & (hist_panda['y'].abs()<2)
    ]

    hist_panda=hist_panda.loc[
        ~hist_panda.index.isin(good_fold_but_large_p_hist.index.tolist())
    ]

    fold_counts_gfblp,_=np.histogram(good_fold_but_large_p_hist['x'],bins=fold_bin_edges)
    fold_counts_log_gfblp=np.log10(fold_counts_gfblp)#,where=[True if element!=0 else False for element in fold_counts])
    fold_counts_0_removed_gfblp=[0.01 if (element==0) else element for element in fold_counts_log_gfblp]
    fold_counts_log_neg_inf_removed_gfblp=[0 if (np.isneginf(element)==True) else element for element in fold_counts_0_removed_gfblp]
    # fold_counts_log_2=fold_counts_log
    
    hovertext_values_gfblp=list()
    for i in range(len(x_positions)):
        hovertext_values_gfblp.append(
            str(fold_bin_edges[i]).split('.')[0]+'.'+str(fold_bin_edges[i]).split('.')[1][:2]+' to '+str(fold_bin_edges[i+1]).split('.')[0]+'.'+str(fold_bin_edges[i+1]).split('.')[1][:2]+': '+str(fold_counts_gfblp[i])+' compounds'
        )
    #######################



    fold_counts,_=np.histogram(hist_panda['x'],bins=fold_bin_edges)
    fold_counts_log=np.log10(fold_counts)#,where=[True if element!=0 else False for element in fold_counts])
    fold_counts_0_removed=[0.01 if (element==0) else element for element in fold_counts_log]
    fold_counts_log_neg_inf_removed=[0 if (np.isneginf(element)==True) else element for element in fold_counts_0_removed]
    # fold_counts_log_2=fold_counts_log
    
    hovertext_values=list()
    for i in range(len(x_positions)):
        hovertext_values.append(
            str(fold_bin_edges[i]).split('.')[0]+'.'+str(fold_bin_edges[i]).split('.')[1][:2]+' to '+str(fold_bin_edges[i+1]).split('.')[0]+'.'+str(fold_bin_edges[i+1]).split('.')[1][:2]+': '+str(fold_counts[i])+' compounds'
        )



    x_hist=go.Figure()

    ###########################
    x_hist.add_trace(
        go.Bar(
            x=x_positions,
            y=fold_counts_log_neg_inf_removed_gfblp,
            width=widths,
            marker_color=marker_colors_gfblp,

            hovertext=hovertext_values_gfblp,
            hoverinfo='text'
        )

    )
    #############################



    x_hist.add_trace(
        go.Bar(
            x=x_positions,
            y=fold_counts_log_neg_inf_removed,
            width=widths,
            marker_color=marker_colors,

            hovertext=hovertext_values,
            hoverinfo='text'
        )

    )










    x_hist.update_layout(
        title_text='Fold Changes', 
        # title_x=0.5,
        showlegend=False,
        # plot_bgcolor='#FFFFFF',
        xaxis_title_text='log2(Fold Change)',
        yaxis_title_text='log10(Compound Count)',
        bargap=0,
        margin=dict(l=20, r=20, t=40, b=20),
        font=dict(size=20,color='black',family='Roboto'),
        barmode='stack'
    )

    return x_hist

def create_pvalue_hist(hist_panda):
    p_bin_edges=list()
    # print('%'*50)
    # print(hist_panda)
    # print(hist_panda['y'].min())
    # print(hist_panda['y'].max())

    #the "insignificant"
    if hist_panda['y'].min()<1:
        if hist_panda['y'].max()>=1:
            upper_limit=1
        else:
            upper_limit= hist_panda['y'].max()
        p_bin_edges+=list(np.linspace(hist_panda['y'].min(),upper_limit,9))
    #the "significant"
    if hist_panda['y'].max()>1:
        if hist_panda['y'].min()<1:
            lower_limit=1
        else:
            lower_limit=hist_panda['y'].min()
        p_bin_edges+=list(np.linspace(lower_limit,hist_panda['y'].max(),10))        
    


    p_bin_edges=sorted(list(set(p_bin_edges)))
    # print(p_bin_edges)


    widths=[(p_bin_edges[i+1]-p_bin_edges[i]) for i in range(len(p_bin_edges)-1)]
    # print('widths')
    # print(widths)
    y_positions=[(0.5*p_bin_edges[i]+0.5*p_bin_edges[i+1]) for i in range(len(p_bin_edges)-1)]
    # print(x_positions)

    ######################
    marker_colors_gpbsf=[  "#2884f4" for temp_pos in y_positions     ]
    good_p_but_small_fold_hist=hist_panda.loc[
        (hist_panda['x'].abs()<1) & (hist_panda['y'].abs()>2)
    ]

    hist_panda=hist_panda.loc[
        ~hist_panda.index.isin(good_p_but_small_fold_hist.index.tolist())
    ]

    p_counts_gpbsf,_=np.histogram(good_p_but_small_fold_hist['y'],bins=p_bin_edges)
    p_counts_log_gpbsf=np.log10(p_counts_gpbsf)#,where=[True if element!=0 else False for element in fold_counts])
    p_counts_0_removed_gpbsf=[0.01 if (element==0) else element for element in p_counts_log_gpbsf]
    p_counts_log_neg_inf_removed_gpbsf=[0 if (np.isneginf(element)==True) else element for element in p_counts_0_removed_gpbsf]
    # fold_counts_log_2=fold_counts_log
    
    hovertext_values_gpbsf=list()
    for i in range(len(y_positions)):
        hovertext_values_gpbsf.append(
            str(p_bin_edges[i]).split('.')[0]+'.'+str(p_bin_edges[i]).split('.')[1][:2]+' to '+str(p_bin_edges[i+1]).split('.')[0]+'.'+str(p_bin_edges[i+1]).split('.')[1][:2]+': '+str(p_counts_gpbsf[i])+' compounds'
        )
    #######################




    marker_colors=[  "#FF0000" if (temp_pos>1) else "#2884f4" for temp_pos in y_positions     ]

    p_counts,_=np.histogram(hist_panda['y'],bins=p_bin_edges)

    p_counts_log=np.log10(p_counts)#,where=[True if element!=0 else False for element in fold_counts])

    #ordering matters here
    #1 becomes 0, then -inf becomes 0
    p_counts_0_removed=[0.01 if (element==0) else element for element in p_counts_log]
    p_counts_log_neg_inf_removed=[0 if (np.isneginf(element)==True) else element for element in p_counts_0_removed]
    # fold_counts_log_2=fold_counts_log
    
    hovertext_values=list()
    for i in range(len(y_positions)):
        hovertext_values.append(
            str(p_bin_edges[i]).split('.')[0]+'.'+str(p_bin_edges[i]).split('.')[1][:2]+' to '+str(p_bin_edges[i+1]).split('.')[0]+'.'+str(p_bin_edges[i+1]).split('.')[1][:2]+': '+str(p_counts[i])+' compounds'
        )

    
    y_hist=go.Figure()


    y_hist.add_trace(
        go.Bar(
            y=y_positions,
            x=p_counts_log_neg_inf_removed_gpbsf,
            width=widths,
            marker_color=marker_colors_gpbsf,
            orientation='h',
            hovertext=hovertext_values_gpbsf,
            hoverinfo='text'
        )

    )


    y_hist.add_trace(
        go.Bar(
            y=y_positions,
            x=p_counts_log_neg_inf_removed,
            width=widths,
            marker_color=marker_colors,
            orientation='h',
            hovertext=hovertext_values,
            hoverinfo='text'
        )

    )
    y_hist.update_layout(
        title_text='p-Values', 
        # title_x=0.5,
        showlegend=False,
        # plot_bgcolor='#FFFFFF',
        xaxis_title_text='log10(Compound Count)',
        yaxis_title_text='log10(p-Value)',
        bargap=0,
        margin=dict(l=20, r=20, t=40, b=20),
        font=dict(size=20,color='black',family='Roboto'),
        barmode='stack'
    )

    return y_hist

def create_summary_stats(hist_panda):

    datatable_counts,x_,y_=np.histogram2d(
        x=hist_panda['x'],
        y=hist_panda['y'],
        bins=[
            [-9999999,-1,1,99999999],
            [0,1.30103,2,99999999999]
        ]
    )
    # print(datatable_counts)
    #datatable_counts[0]
    
    # datatable_counts=datatable_counts.T

    # datatable_counts[0],datatable_counts[2]=datatable_counts[2],datatable_counts[0]

    # print(datatable_counts)
    # print('++++++++++++++++++++++++++++++++++++++++++++++++++++++')

    summary_stats=pd.DataFrame.from_dict(
        {
            'decreased':np.flip(datatable_counts[0]),
            'neither':np.flip(datatable_counts[1]),
            'increased':np.flip(datatable_counts[2]),
            # 'sig_labels':['p<0.01','0.01<p<0.05','0.05<p']
            'sig_labels':['p<0.01','0.01<p<0.05','p>0.05']
        }
    )
    data=summary_stats.to_dict(orient='records')


    output_table=dash_table.DataTable(
        id='hgda_table_summary',
        columns=[
            {'name': '', 'id': 'sig_labels'},
            {'name': 'Decreased', 'id': 'decreased'}, 
            {'name': 'Neither', 'id': 'neither'},
            {'name': 'Increased', 'id': 'increased'}
        ],
        data=data,
        # page_current=0,
        # page_size=10,
        # page_action='native',
        # sort_action='native',
        # sort_mode='multi',
        # filter_action='native',
        
        style_cell={
            'fontSize': 17,
            'padding': '8px',
            'textAlign': 'center'
        },
        style_header={
            'font-family': 'arial',
            'fontSize': 15,
            'fontWeight': 'bold',
            'textAlign': 'center'
        },
        style_data={
            'textAlign': 'center',
            'fontWeight': 'bold',
            'font-family': 'Roboto',
            'fontSize': 15,
        },
        
        style_data_conditional=[
            {
                'if':{
                    'column_id':'sig_labels'
                },
                'font-family': 'arial',
                'fontSize': 15,
                'fontWeight': 'bold',
                'textAlign': 'center',

                'fontSize': 17,
                'padding': '8px',
                'textAlign': 'center',
                'backgroundColor':'#fafafa'
            }
        ]
    ) 

    return output_table

