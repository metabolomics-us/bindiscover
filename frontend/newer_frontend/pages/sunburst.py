import dash
from dash import dcc, html, dash_table, callback
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.dash_table.Format import Format, Scheme, Group
from dash import callback_context
from dash.exceptions import PreventUpdate

import plotly.express as px
import plotly.graph_objects as go
import requests
import pandas as pd

from . import sunburst_helper

#when containerized, the url is not the local 127.0.0.1
base_url_api = f"http://api_alias:4999/"
# base_url_api = "http://127.0.0.1:4999/"
dash.register_page(__name__,path_template="/sunburst/<linked_compound>")

#populate constants for functionality#########
compound_dropdown_options=sunburst_helper.create_compound_selection_labels("../newer_datasets/compound_translation_for_all_components.bin")
compound_dropdown_options_sorted=sorted(
                                compound_dropdown_options,
                                key=lambda x:x['label']
                            )
sunburst_compound_translation_dict={
    element['value']:element['label'] for element in compound_dropdown_options_sorted
}
##############################################

layout=html.Div(
    children=[
        #start the app on a particular compound for simplicity
        dcc.Location(id='url2',refresh=False),
        dbc.Col(width=3),
        dbc.Row(
            children=[
                dbc.Col(
                    children=[
                        html.Br(),
                    ],
                    width={'size':6}
                )
            ],
            justify='center'
        ),
        dbc.Row(
            children=[
                dbc.Col(width=3),
                dbc.Col(
                    children=[
                        html.H2("Compounds", className='text-center'),
                        dcc.Dropdown(
                            id='compound_selection',
                            multi=False,
                            placeholder='Type compound name to search',
                            options=['Type substring to populate options.']
                        ),
                        html.Br(),
                    ],
                    width={'size':3}
                ),
                dbc.Col(
                    children=[
                        html.H2("Display Options", className='text-center'),
                        html.H6("Metric for Sunburst", className='text-center'),
                        html.Div(className="radio-group-container add-margin-top-1", children=[
                            html.Div(className="radio-group", children=[
                                dbc.RadioItems(
                                    id='radio_items_sunburst_value',
                                    options=[
                                        {'label': 'Average', 'value': 'intensity_average'},
                                        {'label': 'Percent Present', 'value':'percent_present'}
                                    ],         
                                    value='intensity_average',
                                    className="btn-group",
                                    inputClassName="btn-check",
                                    labelClassName="btn btn-outline-primary",
                                    inputCheckedClassName="active"                               
                                ),
                            ])
                        ]),
                        html.Br(),
                        html.H6("Order for Sunburst", className='text-center'),
                        html.Div(className="radio-group-container add-margin-top-1", children=[
                            html.Div(className="radio-group", children=[
                                dbc.RadioItems(
                                    id='radio_items_sod_order',
                                    options=[
                                        {'label': 'Species, Organ, Disease', 'value': 'binvestigate,species,organ,disease'},
                                        {'label': 'Species, Disease, Organ', 'value': 'binvestigate,species,disease,organ'},
                                        {'label': 'Organ, Species, Disease', 'value': 'binvestigate,organ,species,disease'},
                                        {'label': 'Organ, Disease, Species', 'value': 'binvestigate,organ,disease,species'},
                                        {'label': 'Disease, Species, Organ', 'value': 'binvestigate,disease,species,organ'},
                                        {'label': 'Disease, Organ, Species', 'value': 'binvestigate,disease,organ,species'},
                                    ],         
                                    value='binvestigate,species,organ,disease',
                                    className="btn-group",
                                    inputClassName="btn-check",
                                    labelClassName="btn btn-outline-primary",
                                    inputCheckedClassName="active"                               
                                ),
                            ])
                        ]),
                        html.Br(),
                    ],
                    width={'size':3}
                ),
                dbc.Col(width=3),
            ],
        ),
        html.Br(),
        html.Br(),
        dbc.Row(
            children=[
                dbc.Col(
                    children=[
                        #html.H2("Execute or Update Query", className='text-center'),
                        
                        html.Div(
                            dbc.Alert(
                                "WARNING: Comparisons involving multiple types of organs are intrinsically semi-quantitative. Each organ's extraction method has a different sample amount, e.g., 20 µL of plasma or 4 mg of liver.",
                                color='primary'
                            ),
                            className="d-grid gap-2 col-6 mx-auto",
                        ),
                        
                        html.Br(),
                        html.Div(
                            dbc.Button(
                                'Search Compound',
                                id='button_query',
                            ),
                            className="d-grid gap-2 col-6 mx-auto",
                        ),
                        html.Br(),
                    ],
                    width={'size':6}
                )
            ],
            justify='center'
        ),
        html.Br(),
        dbc.Spinner(
            children=[
                dbc.Row(
                    children=[
                        dbc.Col(
                            children=[
                                html.H3(
                                    children=' ',
                                    id='compound_selected_sunburst',
                                    style={'textAlign': 'center'}
                                )
                            ],
                        ),
                    ]
                ),
            ]
        ),
        dbc.Row(
            children=[
                dbc.Col(width=2),
                dbc.Col(
                    dbc.Spinner(
                        children=[
                            html.Div(className="sunburst-container",
                                children=[
                                    dcc.Graph(
                                        id='sunburst_figure'
                                    )
                                ]
                            ),
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
                            dcc.Download(id="download_sunburst_datatable"),
                            dash_table.DataTable(
                                id='sunburst_table',
                                columns=[
                                    {"name": "Species", "id": "bin_id"},
                                    {"name": "Organ", "id": "compound_name"},
                                    {"name": "Disease", "id": "group_1"},
                                    {"name": "Metric", "id": "metric"}
                                ],
                                data=[],
                                page_current=0,  
                                #for whatever reason, case insensitive breaks the magnitude filter                              
                                # filter_options={
                                #     'case':'insensitive',
                                #     'placeholder_text':'Type here to filter'
                                # },
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
                        ],
                    ),
                    width=8
                ),
                dbc.Col(width=2)
            ],
        ),
        html.Br(),
        html.Br(),
    ]
)

@callback(
    Output(component_id='compound_selection',component_property='options'),
    Input("compound_selection", "search_value")
)
def update_options(search_value):
    if not search_value:
        raise PreventUpdate
    return [o for o in compound_dropdown_options_sorted if search_value in o["label"]]


@callback(
    [
        Output(component_id="sunburst_table", component_property="columns"),
        Output(component_id="sunburst_table", component_property="data"),
        Output(component_id='url2', component_property='pathname'),
        Output(component_id='compound_selected_sunburst',component_property='children')
    ],
    [
        Input(component_id='button_query', component_property='n_clicks'),
    ],
    [
        State(component_id='url2', component_property='pathname'),
        State(component_id='compound_selection',component_property='value'),
        State(component_id='radio_items_sunburst_value',component_property='value')
    ],
)
def query_table(button_query_n_clicks,url_pathname,compound_selection_value,radio_items_sunburst_value_value):
    '''
    '''

    if callback_context.triggered[0]['prop_id']=='.' :
        sunburst_output={
            "compound":url_pathname.split('/')[-1]
        }
        output_url=url_pathname.split('/')[-1]
    elif callback_context.triggered[0]['prop_id']=='button_query.n_clicks' and compound_selection_value==None:
        sunburst_output={
            "compound":url_pathname.split('/')[-1]
        }
        output_url=url_pathname.split('/')[-1]
    else:
        sunburst_output={
            "compound":compound_selection_value
        }
        output_url=compound_selection_value

    response = requests.post(base_url_api + "/sunburstresource/", json=sunburst_output)
    total_panda = pd.read_json(response.json(), orient="records")
    total_panda['binvestigate']='binvestigate'

    if radio_items_sunburst_value_value=='intensity_average':
        last_column={'name': 'Average Intensity','id':'intensity_average',"type": "numeric","format": Format(group=Group.yes, precision=2, scheme=Scheme.exponent)}
    elif radio_items_sunburst_value_value=='percent_present':
        last_column={'name': 'Percent Present','id':'percent_present',"type": "numeric","format": Format(group=Group.yes, precision=2, scheme=Scheme.exponent)}
    column_list=[
        {'name':'Species','id':'species'},
        {'name':'Organ','id':'organ'},
        {'name':'Disease','id':'disease'},
        last_column
    ]
    total_panda=total_panda[['binvestigate','species','organ','disease',last_column['id']]]

    data = total_panda.to_dict(orient='records')

    output_compound_name=sunburst_compound_translation_dict[sunburst_output['compound']].split(': ')[1]

    return [column_list,data,'/'+url_pathname.split('/')[1]+'/'+str(output_url),'Compound: ' +output_compound_name]

@callback(
    [
        Output(component_id='sunburst_figure', component_property='figure'),
    ],
    [
        Input(component_id='sunburst_table', component_property='derived_virtual_data'),
        Input(component_id='radio_items_sod_order',component_property='value')
    ],
    [
        State(component_id='radio_items_sunburst_value',component_property='value')
    ],
    prevent_initial_call=True
)
def query_figure(sunburst_table_derived_virtual_data,radio_items_sod_order_value,radio_items_sunburst_value_value):
    '''
    '''
    #get dataframe from derived data
    temp=pd.DataFrame.from_records(sunburst_table_derived_virtual_data)
    #coerce it into sunburst form with helper function
    temp_in_sunburst_form=sunburst_helper.coerce_full_panda(temp,radio_items_sunburst_value_value,radio_items_sod_order_value.split(','))

    if radio_items_sunburst_value_value=='percent_present':
        my_hovertext_values=(temp_in_sunburst_form['id'].str.split('/').str[1:].str.join(' - '))+': '+((temp_in_sunburst_form['average']).astype(str))
    else:
        my_hovertext_values=(temp_in_sunburst_form['id'].str.split('/').str[1:].str.join(' - '))+': '+(temp_in_sunburst_form['average'].map('{:,.0f}'.format))
    current_figure=go.Figure(
        go.Sunburst(
            parents=temp_in_sunburst_form['parent'].to_list(),
            labels=[' ' for x in temp_in_sunburst_form['name'].to_list()],
            values=temp_in_sunburst_form['sum'].to_list(),
            ids=temp_in_sunburst_form['id'].to_list(),
            hovertext=my_hovertext_values,
            hoverinfo='text',
            # scene={'bgcolor':'green'}
            # bgcolor='#fff4e4'
        ),
        layout=go.Layout(height=1000,width=1000)
    )
    current_figure.update_layout(
        margin = dict(t=0, l=0, r=0, b=0),
        hoverlabel=dict(font_size=24),
        #values=dict(separatethousands=True)
    )
    return [current_figure]

@callback(
    [
        Output(component_id="download_sunburst_datatable", component_property="data"),
    ],
    [
        Input(component_id="button_download", component_property="n_clicks"),
    ],
    [
        State(component_id="sunburst_table",component_property="data")
    ],
    prevent_initial_call=True
)
def download_sunburst_datatable(
    download_click,
    table_data
    ):
    '''
    '''
    return [dcc.send_data_frame(
        pd.DataFrame.from_records(table_data).to_excel, "binvestigate_sunburst_datatable.xlsx", sheet_name="sheet_1"
    )]