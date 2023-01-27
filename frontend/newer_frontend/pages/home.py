

import dash
from dash import dcc, html
import plotly.express as px
import dash_bootstrap_components as dbc

dash.register_page(__name__, path='/')

layout = html.Div(
    children=[
        html.Br(),
        html.Br(),
        dbc.Row(
            children=[
                dbc.Col(width=2),
                dbc.Col(
                    children=[
                        html.H3('Publications:'),
                        html.H5('Forthecoming'),
                        html.Br(),
                        html.Br(),
                        html.H3('What is this tool?'),
                        html.H5('BinDiscover enables the easy exploration of data from twenty years of Gas-Chromatograph/Mass-Spectrometry metabolomics performed at the West Coast Metabolomics Center at the  University of California, Davis. We hope that this tool enables rapid hypothesis generation from these data.'),
                        html.Br(),
                        html.Br(), 
                        html.H3('What are the different components?'),
                        html.H5('Ontological differential analysis and Phylometabolomic trees effectively visualize all compounds for many metadata combinations simultaneously.'),
                        html.H5('Differential analysis and upset plots are good for visualizing all compounds for several metadata combinations.'),
                        html.H5('The sunburst diagrams are good for focusing on single compound against all metadata combinations.'),
                        html.H5('The BinBrowser tool is good for exploring chemical metadata and downloading BinDiscover .msp files.'),
                        #html.H3('In the past decade, the field of metabolomics has transformed from an obscure specialty into a major “-omics” platform for studying metabolic processes and biomolecular characterization. However, as a whole the field is still very fractured, as the nature of the instrumentation and of the information produced by the platform essentially creates incompatible “islands” of datasets. This lack of data coherency results in the inability to accumulate a critical mass of metabolomics data that has enabled other –omics platforms to make impactful discoveries and meaningful advances.\n-Titus Mak'),
                        html.Br(),
                        html.Br(),
                        html.H3('What is a Bin?'),
                        html.H5('For the user, a bin is roughly a compound. The word derives from the computational-science algorithm \"Bin Sorting\". This algorithm places detected chromatographic peaks and their corresponding spectra into Bins, each of which have been previously identified as a particular compound. In this way, previous identifications can be propagated on newer studies, which allows the West Coast Metabolomics Center to very efficiently annotate studies.'),
                        html.Br(),
                        html.Br(),    
                        html.H3('Acknowledgements:'),
                        html.H5('The authors would like to thank everyone who has contributed samples to the West Coast Metabolomics Center, which has enabled the synchronization of such diverse data. We would also like to thank the Core lab members who contributed to the instrument collection of this dataset.')
                    ],
                    width=8
                ),
                dbc.Col(width=2)
            ]
        )
    ],
)