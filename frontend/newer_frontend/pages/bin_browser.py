import dash
from dash import dcc, html, dash_table, callback
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.dash_table.Format import Format, Scheme, Group
import dash_bio as dashbio
from dash import callback_context
from dash.exceptions import PreventUpdate

import plotly.express as px
import requests
import pandas as pd
import plotly.graph_objects as go

from . import bin_browser_helper

dash.register_page(__name__,path_template="/bin-browser/<linked_compound>")

base_url_api = f"http://api_alias:4999/"
# base_url_api = "http://127.0.0.1:4999/"


########get things from helper script########
bins_dict=bin_browser_helper.generate_bin_dropdown_options()
compound_classes=bin_browser_helper.generate_compound_classes()
compound_translation_panda=pd.read_pickle('../newer_datasets/compound_translation_for_all_components.bin')
compound_translation_dict=dict(zip(compound_translation_panda.compound_identifier.tolist(),compound_translation_panda.english_name.tolist()))
final_curations=pd.read_pickle('../newer_datasets/compound_translation_for_all_components.bin')
final_curations.drop(['bin_type','english_name'],axis='columns',inplace=True)
final_curations.set_index('compound_identifier',drop=True,inplace=True)


layout=html.Div(
    children=[
        html.Br(),
        html.Br(),
        dcc.Location(id='url',refresh=False),
        dcc.Download(
            id='download_msp_known'
        ),
        dcc.Download(
            id='download_msp_unknown'
        ),
        dbc.Row(
            children=[
                dbc.Col(width=1),
                dbc.Col(
                    children=[
                        html.H2("Compound options", className='text-center'),
                        html.Br(),
                    ],
                    width={'size':5}
                ),
                dbc.Col(
                    children=[
                        html.H2("Download .msp Files", className='text-center'),
                        html.Br(),
                    ],
                    width={'size':4}
                ),
            ],
        ),
        dbc.Row(
            children=[
                dbc.Col(width=2),
                dbc.Col(
                    children=[
                        dcc.Dropdown(
                            id='dropdown_bin',
                            multi=False,
                            placeholder='Type compound name to search',
                            options=['Type substring to populate options.']
                        ),  
                        html.Br(),
                        html.Br(),
                        html.Div(
                            dbc.Button(
                                'Query and Visualize',
                                id='button_bin_visualize',
                            ),
                            className="d-grid gap-2 col-4 mx-auto",
                        ),
                    ],
                    width={'size':3}
                ),
                dbc.Col(width=1),
                dbc.Col(
                    children=[
                        html.Div(
                            dbc.Button(
                                'Download All Identified Compounds .msp',
                                id='button_download_msp_identified',
                            ),
                            className="d-grid gap-2 col-4 mx-auto",
                        ),
                        html.Br(),
                        html.Br(),
                        html.Div(
                            dbc.Button(
                                'Download All Unidentified Compounds .msp',
                                id='button_download_msp_unknown',
                            ),
                            className="d-grid gap-1 col-4 mx-auto",
                        )
                    ],
                    width=4
                ),
            ]
        ),
        dbc.Spinner(
            children=[
                html.Div(
                    id='bin_viewer_dif',
                    children=[]
                )
            ]
        ),
    ]
)

@callback(
    Output("dropdown_bin", "options"),
    Input("dropdown_bin", "search_value")
)
def update_options(search_value):
    if not search_value:
        raise PreventUpdate
    return [o for o in bins_dict if search_value in o["label"]]



@callback(
    [
        Output(component_id='bin_viewer_dif',component_property='children'),
        Output(component_id='url', component_property='pathname'),
        Output(component_id='dropdown_bin', component_property='value')
    ],
    [
        Input(component_id='button_bin_visualize', component_property='n_clicks'),
    ],
    [
        State(component_id='url', component_property='pathname'),
        State(component_id='dropdown_bin', component_property='value')
    ],
)
def query_figure(button_bin_visualize_n_clicks,url_pathname,dropdown_bin_value):
    '''
    '''
    if callback_context.triggered[0]['prop_id']=='.':
        bin_output={
            'bin_id':url_pathname.split('/')[-1]
        }
        output_url=url_pathname.split('/')[-1]
    else:
        bin_output={
            'bin_id':dropdown_bin_value
        }
        output_url=dropdown_bin_value


    response = requests.post(base_url_api + "/binresource/", json=bin_output)
    total_panda = pd.read_json(response.json(), orient="records")

    mzs=[float(x.split(':')[0]) for x in total_panda.at[0,'spectrum'].split(' ')]
    intensities=[float(x.split(':')[1]) for x in total_panda.at[0,'spectrum'].split(' ')]
    mzs=[0]+mzs
    intensities=[0]+intensities

    spectrum_figure=go.Figure(
        go.Bar(
            x=mzs,
            y=intensities,
            marker=dict(color="rgb(220, 53, 69)")
        )
    )
    #spectrum_figure.update_title(title="Relative Intensity")
    spectrum_figure.update_yaxes(title="Relative Intensity")
    spectrum_figure.update_xaxes(title="m/z")
    spectrum_figure.update_traces(width=1, hovertemplate="m/z: %{x}<br>Intensity: %{y}<br>")
    spectrum_figure.update_layout(showlegend=False,font=dict(size=18))
    compound_name=total_panda.at[0,'english_name']
    spectrum_figure.update_layout(title={'text':compound_name,'x':0.5})

    total_panda=total_panda[['english_name','compound_identifier','retentionIndex','kovats','spectrum','quantMass','uniqueMass','splash']]
    
    
    total_panda.rename(
        {
            'english_name':'Name',
            'compound_identifier':'InChIKey',
            'retentionIndex':'FAME RI',
            'kovats':'Kovats RI',
            'spectrum':'Spectrum',
            'quantMass':'Quant Mass',
            'uniqueMass':'Unique Mass',
            'splash':'SPLASH',
        },
        axis='columns',
        inplace=True
    )
    total_panda=total_panda.T
    total_panda.at['InChIKey',0]=final_curations.at[
        str(total_panda.at['InChIKey',0]),'identifier'
    ]
    try:
        total_panda.at['Superclass',0]=compound_classes.at[
            total_panda.at['InChIKey',0],'Superclass'
        ]
        total_panda.at['Class',0]=compound_classes.at[
            total_panda.at['InChIKey',0],'Class'
        ]
        total_panda.at['Subclass',0]=compound_classes.at[
            total_panda.at['InChIKey',0],'Subclass'
        ]
    except KeyError:
        total_panda.at['Superclass',0]=''
        total_panda.at['Class',0]=''      
        total_panda.at['Subclass',0]=''       

    total_panda.reset_index(inplace=True)
    total_panda.rename(
        {
            'index':'Attribute',
            0:'Value'
        },
        axis='columns',
        inplace=True
    )

    total_panda.at[0,'Value']=compound_translation_dict[
        str(bin_output['bin_id'])
    ]

    data = total_panda.to_dict(orient='records')
    column_list=[
        {'name': temp_column, 'id':temp_column} for temp_column in total_panda.columns
    ]

    bin_viewer_div_children=[
        html.Br(),
        html.Br(),
        dbc.Row(
            children=[
                dbc.Col(width=1),
                dbc.Col(width=1), 
                dbc.Col(
                    children=[
                        
                        html.Br(),
                        html.Br(),
                        html.H3('Chemical Metadata', className='text-center'),
                        html.Div(
                            dbc.Button(
                                'Download as .xlsx',
                                id='button_download_msp',
                            ),
                            className="d-grid gap-2 col-4 mx-auto",
                        ),
                        html.Br(),
                        dcc.Download(id="download_download_msp"),
                        dash_table.DataTable(
                            id='table_bin',
                            columns=column_list,
                            data=data,
                            style_table={'overflowX': 'scroll'},
                            style_cell={
                                'fontSize': 17,
                                'padding': '8px',
                            },
                            style_header={
                                'font-family': 'arial',
                                'fontSize': 15,
                                'fontWeight': 'bold',
                                'textAlign':'left'
                            },
                            style_data={
                                'textAlign': 'left',
                                'fontWeight': 'bold',
                                'font-family': 'Roboto',
                                'fontSize': 15,
                            },
                        )
                    ],
                    width=3
                ),        
                #dbc.Col(width=1),           
                dbc.Col(
                    children=[
                        html.Br(),
                        html.Br(),
                        html.H3('Mass Spectrum', className='text-center'),
                        dcc.Graph(
                            id='figure_bin',
                            figure=spectrum_figure
                        )                     
                    ]
                ),
                dbc.Col(width=1)
            ]
        ),
    ]
    
    return [bin_viewer_div_children,'/'+url_pathname.split('/')[1]+'/'+str(output_url),bin_output['bin_id']]


@callback(
    [
        Output(component_id="download_download_msp", component_property="data"),
    ],
    [
        Input(component_id="button_download_msp", component_property="n_clicks"),
    ],
    [
        State(component_id='table_bin', component_property='derived_virtual_data')
    ],
    prevent_initial_call=True
)
def download_bin_datatable(
    button_download_msp_n_clicks,
    table_bin_derived_virual_data
    ):
        """
        """

        return [dcc.send_data_frame(
            pd.DataFrame.from_records(table_bin_derived_virual_data).to_excel, "binvestigate_sunburst_datatable.xlsx", sheet_name="sheet_1"
        )]


@callback(
    [
        Output(component_id="download_msp_known", component_property="data"),
    ],
    [
        Input(component_id="button_download_msp_identified", component_property="n_clicks"),
    ],
    prevent_initial_call=True
)
def download_msp_known(    button_download_msp_identified_n_clicks    ):
    return [dcc.send_file(
        '../newer_datasets/GCBinbase_knowns_curated.msp'
    )]

@callback(
    [
        Output(component_id="download_msp_unknown", component_property="data"),
    ],
    [
        Input(component_id="button_download_msp_unknown", component_property="n_clicks"),
    ],
    prevent_initial_call=True
)
def download_msp_unknown(    button_download_msp_unknown_n_clicks    ):
    return [dcc.send_file(
        '../newer_datasets/GCBinbase_unknowns.msp'
    )]