import dash
from dash import dcc, html, dash_table, callback
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.dash_table.Format import Format, Scheme, Group
import dash_bio as dashbio
from dash.exceptions import PreventUpdate
import pandas as pd
import requests
import plotly.express as px
import time
from . import venn_helper

dash.register_page(__name__)

#when containerized, the url is not the local 127.0.0.1
base_url_api = f"http://api_alias:4999/"
#base_url_api = "http://127.0.0.1:4999/"

#populate constants for functionality#########
unique_sod_combinations_dict=venn_helper.get_unique_sod_combinations()
compound_translation_panda=pd.read_pickle('../newer_datasets/compound_translation_for_all_components.bin')
hyperlink_translation_dict=dict(zip(compound_translation_panda.integer_representation.tolist(),compound_translation_panda.compound_identifier.tolist()))
#############################################

#layout=dbc.Container(
layout=html.Div(
    children=[
       html.Br(),
       html.Br(),
       dbc.Row(
            children=[
                html.H2('Metadata Triplets and Compound Types'),
            ],
            style={'textAlign': 'center'}
        ),
        html.Br(),
        dbc.Row(
            children=[
                dbc.Col(width=1),
                dbc.Col(
                    children=[
                        dcc.Dropdown(
                            id='dropdown_triplet_selection_from',
                            options=sorted([
                                {'label': temp, 'value':unique_sod_combinations_dict[temp]} for temp in unique_sod_combinations_dict
                            ],key=lambda x:x['label']),
                            multi=True,
                            placeholder='Select Triplet'
                        ),  
                        html.Br(),
                    ],
                    width={'size':4}
                ),
                dbc.Col(
                    children=[
                        dcc.Dropdown(
                            id='dropdown_triplet_selection_to',
                            options=sorted([
                                {'label': temp, 'value':unique_sod_combinations_dict[temp]} for temp in unique_sod_combinations_dict
                            ],key=lambda x:x['label']),
                            multi=True,
                            placeholder='Select Triplet'
                        ),  
                        html.Br(),
                    ],
                    width={'size':4}
                ),
                dbc.Col(
                    children=[
                        html.Div(className="radio-group-container add-margin-top-1", children=[
                            html.Div(className="radio-group", children=[
                                dbc.RadioItems(
                                    id='radio_items_bin_type',
                                    options=[
                                        {'label': 'Knowns', 'value': 'known'},
                                        {'label': 'Classes', 'value': 'class'},
                                        {'label': 'Unknowns', 'value': 'unknown'},
                                    ],         
                                    value='known',
                                    className="btn-group",
                                    inputClassName="btn-check",
                                    labelClassName="btn btn-outline-primary",
                                    inputCheckedClassName="active",                                
                                ),
                            ])
                        ]),
                    ],
                    width={'size':3}
                ),
            ],
        ),
        dbc.Row(
            children=[
                dbc.Col(width=3),
                dbc.Col(
                    html.Div(
                        dbc.Button(
                            'Determine Valid Triplets',
                            id='metadata_query',
                        ),
                        className="d-grid gap-2 col-3 mx-auto",
                    ),
                ),
                dbc.Col(width=3),
            ]
        ),
        html.Br(),
        html.Br(),
        html.Br(),
        dbc.Spinner(
            children=[
                html.Div(
                    id='div_metadata_selection_dif',
                    children=[]
                )
            ]
        ),
        html.Br(),
        html.Br(),
        html.Br(),
        dbc.Spinner(
            children=[
                html.Div(
                    id='div_volcano_dif',
                    children=[]
                )
            ]
        ),
        dbc.Spinner(
            children=[
                html.Div(
                    id='div_datatable_dif',
                    children=[]
                )
            ]
        )
    ]
)


@callback(
    [
        Output(component_id="div_metadata_selection_dif", component_property="children")
    ],
    [
        Input(component_id='metadata_query', component_property='n_clicks'),
    ],
    [
        State(component_id='dropdown_triplet_selection_from',component_property='value'),
        State(component_id='dropdown_triplet_selection_to',component_property='value'),
    ],
    prevent_initial_call=True
)
def query_md_table(metadata_query_n_clicks,dropdown_triplet_selection_from_value,dropdown_triplet_selection_to_value):
    '''
    determines a set of metadata triplets for which we have data
    from the selected ontological nodes
    '''

    leaf_output={
        "triplet_from":dropdown_triplet_selection_from_value,
        "triplet_to":dropdown_triplet_selection_to_value
    }

    response = requests.post(base_url_api + "/leafmetadataresource/", json=leaf_output)
    total_panda = pd.read_json(response.json(), orient="records")
    data = total_panda.to_dict(orient='records')


    div_metadata_selection_dif_children=[
        dbc.Row(
            children=[
                dbc.Col(width=3),
                dbc.Col(
                    children=[
                        html.H2("Valid Triplets", className='text-center'),
                        dash_table.DataTable(
                            id='table_metadata',
                            columns=[
                                {'name': 'From or To', 'id': 'from_or_to'},
                                {'name': 'Triplet ID', 'id': 'triplet_id'}, 
                                {'name': 'Sample Count', 'id': 'sample_count'}
                            ],
                            data=data,
                            page_current=0,
                            page_size=10,
                            page_action='native',
                            sort_action='native',
                            sort_mode='multi',
                            filter_action='native',
                            row_deletable=True,
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
                            filter_options={
                                'case':'insensitive',
                                'placeholder_text':'Type here to filter'
                            }
                        ),
                        html.Br(),
                        html.Div(
                            dbc.Button(
                                'Perform Differential Analysis',
                                id='leaf_query',
                            ),
                            className="d-grid gap-2 col-3 mx-auto",
                        ),
                    ],
                ),
                dbc.Col(width=3),
            ]
        ),
    ]

    return [div_metadata_selection_dif_children]


@callback(
    [
        Output(component_id='div_volcano_dif', component_property='children'),
    ],
    [
        Input(component_id='leaf_table', component_property='derived_virtual_data'),
    ],
    [
        State(component_id='dropdown_triplet_selection_from',component_property='value'),
        State(component_id='dropdown_triplet_selection_to',component_property='value'),
        State(component_id='radio_items_bin_type',component_property='value'),
    ],
    prevent_initial_call=True
)
def query_figure(leaf_table_derived_virtual_data,dropdown_triplet_selection_from_value,dropdown_triplet_selection_to_value,radio_items_bin_type_value):
    '''
    create the differential analysis volcano plot.
    dot locations come from the datatable
    title comes from the ontological nodes in dropdowns and radio
    '''

    temp=pd.DataFrame.from_records(leaf_table_derived_virtual_data)
    
    if radio_items_bin_type_value!='class':
        temp['english_name']=temp['english_name'].str.extract('\[(.*)\]')

    p='significance_welch'
    effect_size='fold_change_average'


    #done so that we can have symmetric labels
    if len(dropdown_triplet_selection_from_value[0].title()) < len(dropdown_triplet_selection_to_value[0].title()):
        shorter_length=len(dropdown_triplet_selection_from_value[0].title())
        #title=dropdown_triplet_selection_from_value[0][:shorter_length].title()+'             vs.               '+dropdown_triplet_selection_to_value[0][:shorter_length].title()
        x_axis_label='log2 Fold Change<br>Increased in \"'+dropdown_triplet_selection_from_value[0].title()+'\"                                                                      Increased in \"'+dropdown_triplet_selection_to_value[0][:shorter_length].title()+'...\"<br><br>.'

    elif len(dropdown_triplet_selection_from_value[0].title()) > len(dropdown_triplet_selection_to_value[0].title()):
        shorter_length=len(dropdown_triplet_selection_to_value[0].title())   
        x_axis_label='log2 Fold Change<br>Increased in \"'+dropdown_triplet_selection_from_value[0][:shorter_length].title()+'...\"                                                                      Increased in \"'+dropdown_triplet_selection_to_value[0].title()+'\"<br><br>.' 

    elif len(dropdown_triplet_selection_from_value[0].title()) == len(dropdown_triplet_selection_to_value[0].title()):
        x_axis_label='log2 Fold Change<b    r>Increased in \"'+dropdown_triplet_selection_from_value[0].title()+'\"                                                                      Increased in \"'+dropdown_triplet_selection_to_value[0].title()+'\"<br><br>.' 

    volcano = dashbio.VolcanoPlot(
        dataframe=temp,
        snp="english_name",
        p=p,
        effect_size=effect_size,
        gene=None,
        xlabel=x_axis_label,
        genomewideline_value=2,
        title=dropdown_triplet_selection_from_value[0].title()+'             vs.               '+dropdown_triplet_selection_to_value[0].title(),
        title_x=0.5
    )
    volcano.update_layout(showlegend=False)

    div_volcano_dif_children=[
        dbc.Row(
            children=[
                dbc.Col(width={'size':2}),
                dbc.Col(
                    children=[
                        html.H2("Volcano Plot", className='text-center'),
                        dcc.Graph(
                            id='leaf_figure',
                            figure=volcano
                        ),
                    ],
                    width={'size':8}
                ),
                dbc.Col(width={'size':2}),
            ],
        ),        
    ]

    return [div_volcano_dif_children]


@callback(
    [
        Output(component_id="div_datatable_dif", component_property="children")
    ],
    [
        Input(component_id='leaf_query', component_property='n_clicks'),
    ],
    [
        State(component_id='radio_items_bin_type',component_property='value'),
        State(component_id='table_metadata', component_property='derived_virtual_data'),
    ],
    prevent_initial_call=True
)
def query_table(leaf_query_n_clicks,radio_items_bin_type_value,table_metadata_derived_virtual_data):
    '''
    retrieves the differential analysis results for the triplets in the metadata table
    '''

    if leaf_query_n_clicks==None:
        raise PreventUpdate

    leaf_output={
        "metadata_datatable":table_metadata_derived_virtual_data,
        "bin_type":radio_items_bin_type_value
    }

    response = requests.post(base_url_api + "/leafresource/", json=leaf_output)
    total_panda = pd.read_json(response.json(), orient="records")

    if radio_items_bin_type_value!='class':
        total_panda['compound_id']=total_panda['compound_id'].map(hyperlink_translation_dict.get)
        total_panda['english_name']='['+total_panda['english_name']+'](/sunburst/'+total_panda['compound_id'].astype(str)+')'
        total_panda['identifier']='['+total_panda['identifier']+'](/bin-browser/'+total_panda['compound_id'].astype(str)+')'

    data = total_panda.to_dict(orient='records')

    div_datatable_dif_children=[
        dbc.Row(
            children=[
                dbc.Col(width=2),
                dbc.Col(
                        children=[
                            html.Br(),
                            html.Br(),
                            html.H2("Result Datatable", className='text-center'),
                            html.Div(
                                dbc.Button(
                                    'Download Datatable as .xlsx',
                                    id='button_download',
                                ),
                                className="d-grid gap-2 col-3 mx-auto",
                            ),
                            dcc.Download(id="download_leaf_datatable"),
                            dash_table.DataTable(
                                id='leaf_table',
                                columns=[
                                    {"name": "English Name", "id": "english_name",'presentation':'markdown'},
                                    {"name": "Identifier", "id": "identifier",'presentation':'markdown'},
                                    {"name": "log2 Fold Change", "id": "fold_change_average","type": "numeric","format": Format(group=Group.yes, precision=2)},#, scheme=Scheme.exponent)},
                                    {"name": "Significance Welch", "id": "significance_welch","type": "numeric","format": Format(group=Group.yes, precision=2)},#, scheme=Scheme.exponent)},
                                ],
                                markdown_options={"link_target": "_blank"},
                                data=data,
                                page_current=0,
                                page_size=50,
                                page_action='native',
                                sort_action='native',
                                sort_mode='multi',
                                filter_action='native',
                                filter_options={
                                    'case':'insensitive',
                                    'placeholder_text':'Type here to filter'
                                },
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
                            )
                        ],
                        width={'size':8}
                    ),
                dbc.Col(width=2),
            ],
        )
    ]

    return [div_datatable_dif_children]


@callback(
    [
        Output(component_id="download_leaf_datatable", component_property="data"),
    ],
    [
        Input(component_id="button_download", component_property="n_clicks"),
    ],
    [
        State(component_id="leaf_table",component_property="data"),
        State(component_id='radio_items_bin_type',component_property='value')
    ],
    prevent_initial_call=True
)
def download_leaf_datatable(
    download_click,
    table_data,
    radio_items_bin_type_value
    ):
        '''
        download the information in the datatable
        '''

        downloaded_panda=pd.DataFrame.from_records(table_data)

        if radio_items_bin_type_value!='class':
            downloaded_panda['english_name']=downloaded_panda['english_name'].str.extract('\[(.*)\]')
            downloaded_panda['identifier']=downloaded_panda['identifier'].str.extract('\[(.*)\]')

        return [dcc.send_data_frame(
            downloaded_panda.to_excel, "binvestigate_differential_datatable.xlsx", sheet_name="sheet_1"
        )]

