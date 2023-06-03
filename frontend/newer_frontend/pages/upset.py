import dash
from dash import dcc, html, dash_table, callback
import plotly.express as px
import dash_bootstrap_components as dbc
import requests
from . import venn_helper
import pandas as pd
from dash.dependencies import Input, Output, State
from dash.dash_table.Format import Format, Scheme, Group
import xlsxwriter

#when containerized, the url is not the local 127.0.0.1
base_url_api = f"http://api_alias:4999/"
# base_url_api = "http://127.0.0.1:4999/"

dash.register_page(__name__)

#populate constants for functionality#########
unique_sod_combinations_dict=venn_helper.get_unique_sod_combinations()
final_curations=pd.read_pickle('../newer_datasets/compound_translation_for_all_components.bin')
compound_bin_translator_dict=dict(zip(final_curations.loc[final_curations.bin_type=='known']['compound_identifier'].astype(int).tolist(),final_curations.loc[final_curations.bin_type=='known']['english_name'].tolist()))
#############################################

layout=html.Div(
    children=[
        html.Br(),
        dbc.Row(
            children=[
                dbc.Col(
                    children=[
                        html.Br(),
                    ],
                    width={'size':6}
                )
            ],
        ),
        dbc.Row(
            children=[
                dbc.Col(width=2),
                dbc.Col(
                    children=[
                        html.H2("Metadata Combinations", className='text-center'),
                        dcc.Dropdown(
                            id='dropdown_triplet_selection',
                            options=sorted([
                                {'label': temp.title(), 'value':unique_sod_combinations_dict[temp]} for temp in unique_sod_combinations_dict
                            ],key=lambda x:x['label']),
                            multi=True,
                        ),  
                        html.Br(),
                    ],
                    width={'size':4}
                ),
                dbc.Col(
                    children=[
                        html.H2("Options", className='text-center'),
                        html.H6("Display All Compounds or Only Compounds In-common", className='text-center'),
                        html.Div(className="radio-group-container add-margin-top-1", children=[
                            html.Div(className="radio-group", children=[
                                dbc.RadioItems(
                                    id='radio_items_filter',
                                    options=[
                                        {'label': 'No Filter', 'value': 'no_filter'},
                                        {'label': 'Common', 'value': 'common'},
                                    ],         
                                    value='no_filter',
                                    className="btn-group",
                                    inputClassName="btn-check",
                                    labelClassName="btn btn-outline-primary",
                                    inputCheckedClassName="active",                                
                                ),
                            ])
                        ]),
                        html.Br(),
                        html.H6("Display Knowns or Unknowns", className='text-center'),
                        html.Div(className="radio-group-container add-margin-top-1", children=[
                            html.Div(className="radio-group", children=[
                                dbc.RadioItems(
                                    id='radio_items_bin_type_upset',
                                    options=[
                                        {'label': 'Knowns', 'value': 'known'},
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
                        html.Br(),
                        html.H6("What Percent Observed Constitutes Present", className='text-center'),
                        dcc.Slider(
                            id='slider_percent_present',
                            min=0,
                            max=100,
                            step=1,
                            value=50,   
                            marks=None,
                            tooltip={"placement": "bottom", "always_visible": True}       
                        ),
                        html.Br(),
                    ],
                    width={'size':3}
                ),
                dbc.Col(width=2),
            ],
        ),
        html.Br(),
        html.Br(),
        html.Br(),
        dbc.Row(
            children=[
                dbc.Col(width=3),
                dbc.Col(
                    children=[
                        html.Div(
                            dbc.Button(
                                'Search Metadata Combinations',
                                id='button_query',
                            ),
                            className="d-grid gap-2 col-6 mx-auto",
                        ),
                        html.Br(),
                    ],
                    width={'size':6}
                ),
                dbc.Col(width=3),
            ],
        ),
        html.Br(),
        html.Br(),
        dbc.Spinner(
            children=[
                html.Div(
                    id='div_image_upset',
                    children=[]
                )
            ]
        ),
        html.Br(),
        html.Br(),
        html.Br(),
        html.Br(),
        html.Br(),
        dbc.Spinner(
            children=[
                html.Div(
                    id='div_datatable_upset',
                    children=[]
                )
            ]
        ),        
    ]
)


@callback(
    [
        Output(component_id='div_datatable_upset', component_property='children')
    ],
    [
        Input(component_id="button_query", component_property="n_clicks"),
    ],
    [
        State(component_id="dropdown_triplet_selection",component_property="value"),
        State(component_id="slider_percent_present", component_property="value"),
        State(component_id="radio_items_filter",component_property="value"),
        State(component_id='radio_items_bin_type_upset',component_property='value')
    ],
    prevent_initial_call=True
)
def perform_query_table(
    button_query,
    dropdown_triplet_selection_value,
    slider_percent_present_value,
    radio_items_filter_value,
    radio_items_bin_type_upset_value
    ):

        '''
        '''

        ##################volcano query######################
        #prepare json for api
        venn_data_table_output={
            "dropdown_triplet_selection_value":dropdown_triplet_selection_value,
            "slider_percent_present_value":slider_percent_present_value,
            "toggle_average_true_value":True,
            "radio_items_filter_value":radio_items_filter_value,
        }

        response = requests.post(base_url_api + "/venntableresource/", json=venn_data_table_output)
        
        total_panda = pd.read_json(response.json(), orient="records")

        total_panda=total_panda.loc[
            total_panda['bin_type']==radio_items_bin_type_upset_value,
            :
        ]

        if radio_items_bin_type_upset_value=='known':
            total_panda['english_name']=total_panda['bin'].map(compound_bin_translator_dict.get)
            total_panda=total_panda.loc[
                total_panda['english_name'].isnull()==False
            ]

        total_panda.drop(['compound','bin_type','compound_identifier'],axis='columns',inplace=True)

        #prepare columns and data for the table
        column_list = [
            {"name": "Bin", "id": "bin"},
            {"name": "English Name", "id":"english_name",'presentation':'markdown'},
            {"name": "InChIKey","id":"identifier",'presentation':'markdown'}

        ]
        sod_column_list=[
            {"name": temp_column, "id": temp_column,"type": "numeric","format": Format(group=Group.yes, precision=2, scheme=Scheme.exponent)} for temp_column in total_panda.columns 
            if (temp_column != "bin" and temp_column!="english_name" and temp_column!="identifier")
        ]

        column_list+=sod_column_list
        temp_column_names=[element['id'] for element in column_list]
        total_panda=total_panda.loc[
            :,
            temp_column_names
        ].copy()

        total_panda['english_name']='['+total_panda['english_name']+'](/sunburst/'+total_panda['bin'].astype(str)+')'
        total_panda['identifier']='['+total_panda['identifier']+'](/bin-browser/'+total_panda['bin'].astype(str)+')'

        data = total_panda.to_dict(orient='records')

        del total_panda

        div_datatable_upset_children=[
            dbc.Row(
                children=[
                    dbc.Col(width=2),
                    dbc.Col(
                        children=[
                            html.H2("Result Datatable", className='text-center'),
                            html.Div(
                                dbc.Button(
                                    'Download Datatable as .xlsx',
                                    id='button_download',
                                ),
                                className="d-grid gap-2 col-3 mx-auto",
                            ),
                            dcc.Download(id="download_datatable"),
                            dash_table.DataTable(
                                id='table',
                                columns=column_list,
                                data=data,
                                markdown_options={"link_target": "_blank"},
                                page_current=0,
                                page_size=50,
                                page_action='native',
                                sort_action='native',
                                sort_mode='multi',
                                filter_action='native',
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
                        ] 
                    ),
                    dbc.Col(width=2),
                ]
            )
        ]
        return [div_datatable_upset_children]


@callback(
    [
        Output(component_id='div_image_upset',component_property='children')
    ],
    [
        Input(component_id="table",component_property="derived_virtual_data")
    ],
    prevent_initial_call=True
)
def perform_query_diagram(
    table_derived_virtual_data
    ):
        '''
        '''
        temp_img=venn_helper.create_upset(pd.DataFrame.from_records(table_derived_virtual_data).drop(['english_name','identifier','bin'],axis='columns'))
        #return [temp_img,temp_img]

        div_image_upset_children=[
            dbc.Row(
                children=[
                    dbc.Col(width=2),
                    dbc.Col(
                        children=[
                            dbc.Row(html.H2("Upset Plot"),style={'textAlign': 'center'}),
                            dbc.Row(
                                html.Div(className="venn-thumbnail-container",
                                    children=[
                                        html.Img(
                                            id='Img_venn',
                                            height=200,
                                            width=200,
                                            src=temp_img
                                        ),
                                    ]
                                ),
                                style={'textAlign': 'center'}
                            ),
                            dbc.Modal(
                                children=[
                                    dbc.ModalHeader(dbc.ModalTitle("Right Click + Save for High-Res"),close_button=True),
                                    dbc.ModalBody(
                                        html.Img(
                                            id='modal_Img_venn',
                                            style={"height": "40vh"},
                                            src=temp_img
                                        )
                                    ),
                                ],
                                className="modal-overarching",
                                id='modal',
                                centered=True,
                                size='xl',
                                is_open=False,
                                style={"max-width": "none", "width": "90%"}
                            ),
                        ]
                    ),
                    dbc.Col(width=2),
                ]
            )
        ]
        return [div_image_upset_children]




@callback(
    [
        Output(component_id='modal', component_property='is_open'),
    ],
    [
        Input(component_id='Img_venn', component_property='n_clicks'),
    ],
    prevent_initial_call=True
)
def open_modal(Img_venn_n_clicks):
    return [True]

@callback(
    [
        Output(component_id="download_datatable", component_property="data"),
    ],
    [
        Input(component_id="button_download", component_property="n_clicks"),
    ],
    [
        State(component_id="table",component_property="data"),
        State(component_id='radio_items_bin_type_upset',component_property='value')
    ],
    prevent_initial_call=True
)
def download_datatable(
    download_click,
    table_data,
    radio_items_bin_type_upset_value
    ):
        '''
        '''

        downloaded_panda=pd.DataFrame.from_records(table_data)

        if radio_items_bin_type_upset_value!='class':
            downloaded_panda['english_name']=downloaded_panda['english_name'].str.extract('\[(.*)\]')
            downloaded_panda['identifier']=downloaded_panda['identifier'].str.extract('\[(.*)\]')

        return [dcc.send_data_frame(
            downloaded_panda.to_excel, "binvestigate_upset_datatable.xlsx", sheet_name="sheet_1"
        )]

