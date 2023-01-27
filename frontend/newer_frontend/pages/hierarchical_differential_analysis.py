import dash
import dash_bio as dashbio
import dash_bootstrap_components as dbc
from dash import dcc, html, dash_table, callback
from dash.dash_table.Format import Format, Scheme, Group
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import plotly.express as px
import requests
import pandas as pd
import networkx as nx
from . import hierarchical_differential_analysis_helper
import datetime

dash.register_page(__name__)
#when containerized, the url is not the local 127.0.0.1
base_url_api = f"http://api_alias:4999/"
#base_url_api = "http://127.0.0.1:4999/"

#populate constants for functionality#########
#used for dropdown values and dropdown filtering logic (pick a species, get a reduced set of organ options)
species_networkx,species_node_dict=hierarchical_differential_analysis_helper.extract_networkx_selections_species()
organ_networkx,organ_node_dict=hierarchical_differential_analysis_helper.extract_networkx_selections_organ()
disease_networkx,disease_node_dict=hierarchical_differential_analysis_helper.extract_networkx_selections_disease()
#translate species id/english-names
index_panda=pd.read_pickle('../newer_datasets/index_panda.bin')
index_panda=index_panda.sort_index()
index_panda['species']=index_panda['species'].astype(str)
#translate compound names
compound_translation_panda=pd.read_pickle('../newer_datasets/compound_translation_for_all_components.bin')
hyperlink_translation_dict=dict(zip(compound_translation_panda.integer_representation.tolist(),compound_translation_panda.compound_identifier.tolist()))
##############################################

layout=html.Div(
    children=[
        html.Br(),
        html.Br(),
        dbc.Row(
            children=[
                html.H2('Ontological Nodes and Compound Types'),
            ],
            style={'textAlign': 'center'}
        ),
        html.Br(),
        dbc.Row(
            children=[
                #blank columns to buffer whitespace
                dbc.Col(width=1),
                dbc.Col(
                    children=[
                        dcc.Dropdown(
                            id='dropdown_from_species',
                            options=sorted([
                                {'label':species_node_dict[temp], 'value':temp.title()} for temp in species_node_dict
                            ],key=lambda x:x['label']),
                            multi=False,
                            placeholder='Select species ontology node'
                        ),  
                        dcc.Dropdown(
                            id='dropdown_from_organ',
                            options=sorted([
                                {'label':organ_node_dict[temp], 'value':temp} for temp in organ_node_dict
                            ],key=lambda x:x['label']),
                            multi=False,
                            placeholder='Select organ ontology node'
                        ), 
                        dcc.Dropdown(
                            id='dropdown_from_disease',
                            options=sorted([
                                {'label':disease_node_dict[temp], 'value':temp} for temp in disease_node_dict
                            ],key=lambda x:x['label']),
                            multi=False,
                            placeholder='Select disease ontology node',
                            value='No Disease'
                        ), 
                        html.Br(),
                    ],
                    width={'size':3}
                ),
                dbc.Col(
                    children=[
                        dcc.Dropdown(
                            id='dropdown_to_species',
                            options=sorted([
                                {'label':species_node_dict[temp], 'value':temp.title()} for temp in species_node_dict
                            ],key=lambda x:x['label']),
                            multi=False,
                            placeholder='Select species ontology node'
                        ),  
                        dcc.Dropdown(
                            id='dropdown_to_organ',
                            options=sorted([
                                {'label':organ_node_dict[temp], 'value':temp} for temp in organ_node_dict
                            ],key=lambda x:x['label']),
                            multi=False,
                            placeholder='Select organ ontology node'
                        ), 
                        dcc.Dropdown(
                            id='dropdown_to_disease',
                            options=sorted([
                                {'label':disease_node_dict[temp], 'value':temp} for temp in disease_node_dict
                            ],key=lambda x:x['label']),
                            multi=False,
                            placeholder='Select disease ontology node',
                            value='No Disease'
                        ), 
                        html.Br(),
                    ],
                    width={'size':3}
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
                            id='hgda_metadata_query',
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
                    id='div_metadata_selection_onto',
                    children=[]
                )
            ]
        ),
        html.Div(
            id='div_metadata_time_estimator',
            children=[]
        ),
        html.Br(),
        html.Br(),
        html.Br(),
        dbc.Spinner(
            children=[
                html.Div(
                    id='div_volcano_onto',
                    children=[]
                )
            ]
        ),
        dbc.Spinner(
            children=[
                html.Div(
                    id='div_datatable_onto',
                    children=[]
                )
            ]
        )
    ]
)


@callback(
    [
        Output(component_id="dropdown_from_species", component_property="options"),
        Output(component_id="dropdown_from_organ", component_property="options"),
        Output(component_id="dropdown_from_disease", component_property="options"),
    ],
    [
        Input(component_id="dropdown_from_species", component_property="value"),
        Input(component_id="dropdown_from_organ", component_property="value"),
        Input(component_id="dropdown_from_disease", component_property="value"),
    ],
    prevent_initial_call=True
)
def update_input_options_from(
    from_species_value_input,
    from_organ_value_input,
    from_disease_value_input,
):
    '''
    Filter species, organs, diseases if one or more of those option types is specified
    eg, selecting bacteria narrows organ options to things like "cells"
    the left-hand-side option set
    '''

    temp_view=index_panda.copy()
    if from_species_value_input!=None:
        temp_set=nx.algorithms.dag.descendants(species_networkx,from_species_value_input)
        temp_set.add(from_species_value_input)
        temp_view=temp_view.loc[
            temp_view.species.isin(temp_set)
        ]

    if from_organ_value_input!=None:
        temp_set=nx.algorithms.dag.descendants(organ_networkx,from_organ_value_input)
        temp_set.add(from_organ_value_input)
        temp_view=temp_view.loc[
            temp_view.organ.isin(temp_set)
        ]

    if from_disease_value_input!=None:
        temp_set=nx.algorithms.dag.descendants(disease_networkx,from_disease_value_input)
        temp_set.add(from_disease_value_input)
        temp_view=temp_view.loc[
            temp_view.disease.isin(temp_set)
        ]

    all_basic_species_options=set(temp_view.species.values)
    all_valid_species_options=all_basic_species_options
    [all_valid_species_options:=all_valid_species_options.union(nx.ancestors(species_networkx,temp_option)) for temp_option in all_basic_species_options]

    all_basic_organ_options=set(temp_view.organ.values)
    all_valid_organ_options=all_basic_organ_options
    [all_valid_organ_options:=all_valid_organ_options.union(nx.ancestors(organ_networkx,temp_option)) for temp_option in all_basic_organ_options]

    all_basic_disease_options=set(temp_view.disease.values)
    all_valid_disease_options=all_basic_disease_options
    [all_valid_disease_options:=all_valid_disease_options.union(nx.ancestors(disease_networkx,temp_option)) for temp_option in all_basic_disease_options]

    species_options=sorted([
        {'label':species_node_dict[temp], 'value':temp.title()} for temp in species_node_dict if temp in all_valid_species_options
    ],key=lambda x:x['label'])

    organ_options=sorted([
        {'label':organ_node_dict[temp], 'value':temp} for temp in organ_node_dict if temp in all_valid_organ_options
    ],key=lambda x:x['label'])

    disease_options=sorted([
        {'label':disease_node_dict[temp], 'value':temp} for temp in disease_node_dict if temp in all_valid_disease_options
    ],key=lambda x:x['label'])

    return species_options,organ_options,disease_options


@callback(
    [
        Output(component_id="dropdown_to_species", component_property="options"),
        Output(component_id="dropdown_to_organ", component_property="options"),
        Output(component_id="dropdown_to_disease", component_property="options"),
    ],
    [
        Input(component_id="dropdown_to_species", component_property="value"),
        Input(component_id="dropdown_to_organ", component_property="value"),
        Input(component_id="dropdown_to_disease", component_property="value"),
    ],
    prevent_initial_call=True
)
def update_input_options_to(
    to_species_value_input,
    to_organ_value_input,
    to_disease_value_input,
):
    '''
    Filter species, organs, diseases if one or more of those option types is specified
    eg, selecting bacteria narrows organ options to things like "cells"
    the right hand side option set
    '''

    temp_view=index_panda.copy()
    if to_species_value_input!=None:
        temp_set=nx.algorithms.dag.descendants(species_networkx,to_species_value_input)
        temp_set.add(to_species_value_input)
        temp_view=temp_view.loc[
            temp_view.species.isin(temp_set)
        ]

    if to_organ_value_input!=None:
        temp_set=nx.algorithms.dag.descendants(organ_networkx,to_organ_value_input)
        temp_set.add(to_organ_value_input)
        temp_view=temp_view.loc[
            temp_view.organ.isin(temp_set)
        ]

    if to_disease_value_input!=None:
        temp_set=nx.algorithms.dag.descendants(disease_networkx,to_disease_value_input)
        temp_set.add(to_disease_value_input)
        temp_view=temp_view.loc[
            temp_view.disease.isin(temp_set)
        ]

    all_basic_species_options=set(temp_view.species.values)
    all_valid_species_options=all_basic_species_options
    [all_valid_species_options:=all_valid_species_options.union(nx.ancestors(species_networkx,temp_option)) for temp_option in all_basic_species_options]

    all_basic_organ_options=set(temp_view.organ.values)
    all_valid_organ_options=all_basic_organ_options
    [all_valid_organ_options:=all_valid_organ_options.union(nx.ancestors(organ_networkx,temp_option)) for temp_option in all_basic_organ_options]

    all_basic_disease_options=set(temp_view.disease.values)
    all_valid_disease_options=all_basic_disease_options
    [all_valid_disease_options:=all_valid_disease_options.union(nx.ancestors(disease_networkx,temp_option)) for temp_option in all_basic_disease_options]

    species_options=sorted([
        {'label':species_node_dict[temp], 'value':temp.title()} for temp in species_node_dict if temp in all_valid_species_options
    ],key=lambda x:x['label'])

    organ_options=sorted([
        {'label':organ_node_dict[temp], 'value':temp} for temp in organ_node_dict if temp in all_valid_organ_options
    ],key=lambda x:x['label'])

    disease_options=sorted([
        {'label':disease_node_dict[temp], 'value':temp} for temp in disease_node_dict if temp in all_valid_disease_options
    ],key=lambda x:x['label'])

    return species_options,organ_options,disease_options


@callback(
    [
        Output(component_id="div_metadata_selection_onto", component_property="children")
    ],
    [
        Input(component_id="hgda_metadata_query", component_property="n_clicks"),
    ],
    [
        State(component_id="dropdown_from_species", component_property="value"),
        State(component_id="dropdown_from_organ", component_property="value"),
        State(component_id="dropdown_from_disease", component_property="value"),
        State(component_id="dropdown_to_species", component_property="value"),
        State(component_id="dropdown_to_organ", component_property="value"),
        State(component_id="dropdown_to_disease", component_property="value"),
    ],
    prevent_initial_call=True
)
def perform_metadata_query(
    query,
    from_species_value,
    from_organ_value,
    from_disease_value,
    to_species_value,
    to_organ_value,
    to_disease_value,
):
    '''
    determines a set of metadata triplets for which we have data
    from the selected ontological nodes
    '''

    metadata_json_output = {
        "from_species": from_species_value,
        "from_organ": from_organ_value,
        "from_disease": from_disease_value,
        "to_species": to_species_value,
        "to_organ": to_organ_value,
        "to_disease": to_disease_value,
    }

    #obtain results from api
    response = requests.post(base_url_api + "/hgdametadataresource/", json=metadata_json_output)
    total_panda = pd.read_json(response.json(), orient="records")
    data = total_panda.to_dict(orient='records')

    div_metadata_selection_onto_children=dbc.Row(
        children=[
            dbc.Col(width=3),
            dbc.Col(
                children=[
                    html.H2("Valid Triplets", className='text-center'),
                    dash_table.DataTable(
                        id='hgda_table_metadata',
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
                        row_deletable=True,
                    ),
                    html.Br(),
                    html.Div(
                        dbc.Button(
                            'Perform Differential Analysis',
                            id='hgda_query',
                        ),
                        className="d-grid gap-2 col-3 mx-auto",
                    ),
                ],
            ),
            dbc.Col(width=3),
        ]
    )

    return [div_metadata_selection_onto_children]


@callback(
    [
        Output(component_id="div_metadata_time_estimator", component_property="children")
    ],
    [
        Input(component_id='radio_items_bin_type',component_property='value'),
        Input(component_id='hgda_table_metadata', component_property='derived_virtual_data')
    ],
    prevent_initial_call=True
)
def update_time_requirement_estimate(radio_items_bin_type_value,hgda_table_metadata_derived_virtual_data):
    
    multiple_metadata_panda=pd.DataFrame.from_dict(hgda_table_metadata_derived_virtual_data)
    number_from=multiple_metadata_panda.from_or_to.value_counts()['from']
    number_to=multiple_metadata_panda.from_or_to.value_counts()['to']

    if radio_items_bin_type_value=='known':
        seconds_per=1
    elif radio_items_bin_type_value=='class':
        seconds_per=0.3
    elif radio_items_bin_type_value=='unknown':
        seconds_per=7
    total_time=number_from*number_to*seconds_per
    total_time/=60
    total_time=str(datetime.timedelta(minutes=total_time))
    total_string=f'The estimated duration of this query is {total_time} hours, minutes, and seconds.'

    div_metadata_time_estimator_children=[
        dbc.Row(
            children=[
                dbc.Col(width={'size':2}),
                dbc.Col(
                    children=[
                        html.Br(),
                        html.H4(total_string, className='text-center'),
                    ],
                    width={'size':8}
                ),
                dbc.Col(width={'size':2}),
            ],
        ),
    ]    

    return [div_metadata_time_estimator_children]


@callback(
    [
        Output(component_id="div_datatable_onto", component_property="children")
    ],
    [
        Input(component_id='hgda_query', component_property='n_clicks'),
    ],
    [
        State(component_id='radio_items_bin_type',component_property='value'),
        State(component_id='hgda_table_metadata', component_property='derived_virtual_data')
    ],
    prevent_initial_call=True
)
def query_table(
    hgda_query_n_clicks,
    radio_items_bin_type_value,
    hgda_table_metadata_derived_virtual_data
):
    '''
    retrieves the differential analysis results for the triplets in the metadata table
    '''

    if hgda_query_n_clicks==None:
        raise PreventUpdate

    json_output = {
        'metadata_datatable':hgda_table_metadata_derived_virtual_data,
        "bin_type":radio_items_bin_type_value
    }

    #obtain results from api
    response = requests.post(base_url_api + "/hgdaresource/", json=json_output)
    total_panda = pd.read_json(response.json(), orient="records")

    #compounds in datatable stored as integers. translate integers to english names
    if radio_items_bin_type_value!='class':
        total_panda['compound_id']=total_panda['compound_id'].map(hyperlink_translation_dict.get)
        total_panda['english_name']='['+total_panda['english_name']+'](/sunburst/'+total_panda['compound_id'].astype(str)+')'
        total_panda['identifier']='['+total_panda['identifier']+'](/bin-browser/'+total_panda['compound_id'].astype(str)+')'

    data = total_panda.to_dict(orient='records')
    
    div_datatable_onto_children=[
        dbc.Row(
            children=[
                dbc.Col(width={'size':2}),
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
                        dcc.Download(id="download_hgda_datatable"),
                        dash_table.DataTable(
                            id='hgda_table',
                            columns=[
                                {"name": "English Name", "id": "english_name",'presentation':'markdown'},
                                {"name": "Identifier", "id": "identifier",'presentation':'markdown'},
                                {"name": "log2 Fold Change", "id": "fold_change_average","type": "numeric","format": Format(group=Group.yes, precision=2, scheme=Scheme.exponent)},
                                {"name": "Significance Welch", "id": "significance_welch","type": "numeric","format": Format(group=Group.yes, precision=2, scheme=Scheme.exponent)},
                            ],
                            markdown_options={"link_target": "_blank"},
                            data=data,
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
                    ],
                    width={'size':8}
                ),
                dbc.Col(width={'size':2}),
            ],
        ),
    ]
    return [div_datatable_onto_children]
        


@callback(
    [
        Output(component_id='div_volcano_onto', component_property='children'),
    ],
    [
        Input(component_id='hgda_table', component_property='derived_virtual_data'),
    ],
    [
        State(component_id='dropdown_from_species',component_property='value'),
        State(component_id='dropdown_from_organ',component_property='value'),
        State(component_id='dropdown_from_disease',component_property='value'),
        State(component_id='dropdown_to_species',component_property='value'),
        State(component_id='dropdown_to_organ',component_property='value'),
        State(component_id='dropdown_to_disease',component_property='value'),
        State(component_id='radio_items_bin_type',component_property='value'),

    ],
    prevent_initial_call=True
)
def query_figure(
    hgda_table_derived_virtual_data,
    dropdown_from_species_value,
    dropdown_from_organ_value,
    dropdown_from_disease_value,
    dropdown_to_species_value,
    dropdown_to_organ_value,
    dropdown_to_disease_value,
    radio_items_bin_type_value
):
    '''
    create the differential analysis volcano plot.
    dot locations come from the datatable
    title comes from the ontological nodes in dropdowns and radio
    '''

    temp=pd.DataFrame.from_records(hgda_table_derived_virtual_data)

    if radio_items_bin_type_value!='class':
        temp['english_name']=temp['english_name'].str.extract('\[(.*)\]')

    p='significance_welch'
    effect_size='fold_change_average'
        
    title_string_from=' - '.join([species_node_dict[dropdown_from_species_value],organ_node_dict[dropdown_from_organ_value].split(' - ')[0],disease_node_dict[dropdown_from_disease_value].split(' - ')[0]])
    title_string_to=' - '.join([species_node_dict[dropdown_to_species_value],organ_node_dict[dropdown_to_organ_value].split(' - ')[0],disease_node_dict[dropdown_to_disease_value].split(' - ')[0]])

    #done so that we can have symmetric labels
    if len(title_string_from) < len(title_string_to):
        shorter_length=len(title_string_from)
        #title=dropdown_triplet_selection_from_value[0][:shorter_length].title()+'             vs.               '+dropdown_triplet_selection_to_value[0][:shorter_length].title()
        x_axis_label='log2 Fold Change<br>Increased in \"'+title_string_from.title()+'\"                                                                      Increased in \"'+title_string_to[:shorter_length].title()+'...\"<br><br>.'

    elif len(title_string_from) > len(title_string_to):
        shorter_length=len(title_string_to.title())   
        x_axis_label='log2 Fold Change<br>Increased in \"'+title_string_from[:shorter_length].title()+'...\"                                                                      Increased in \"'+title_string_to.title()+'\"<br><br>.' 

    elif len(title_string_from) == len(title_string_to):
        x_axis_label='log2 Fold Change<br>Increased in \"'+title_string_from.title()+'\"                                                                      Increased in \"'+title_string_to.title()+'\"<br><br>.' 



    volcano = dashbio.VolcanoPlot(
        dataframe=temp,
        snp="english_name",
        p=p,
        effect_size=effect_size,
        gene=None,
        xlabel=x_axis_label,
        genomewideline_value=2,
        title=title_string_from+'        vs.       '+title_string_to,
        title_x=0.5
    )
    volcano.update_layout(showlegend=False)

    div_volcano_onto_children=[
        dbc.Row(
            children=[
                dbc.Col(width={'size':2}),
                dbc.Col(
                    children=[
                        html.H2("Volcano Plot", className='text-center'),
                    ],
                    width={'size':8}
                ),
                dbc.Col(width={'size':2}),
            ],
        ),
        dbc.Row(
            children=[
                dbc.Col(width={'size':2}),
                dbc.Col(
                    children=[
                        dcc.Graph(
                            id='hgda_figure',
                            figure=volcano
                        ),
                    ],
                    width={'size':8}
                ),
                dbc.Col(width={'size':2}),
            ],
        ),
    ]
    return [div_volcano_onto_children]


@callback(
    [
        Output(component_id="download_hgda_datatable", component_property="data"),
    ],
    [
        Input(component_id="button_download", component_property="n_clicks"),
    ],
    [
        State(component_id="hgda_table",component_property="data"),
        State(component_id='radio_items_bin_type',component_property='value')
    ],
    prevent_initial_call=True
)
def download_hgda_datatable(
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
