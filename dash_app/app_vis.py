from math import sqrt

import dash
from dash.dependencies import Output, Input  # State
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import dash_html_components as html
import networkx as nx
import pandas as pd

import bio_networks.network_generator as ng
import bio_networks.helpers as h
import dash_app.cytoscape_stylesheets as stylesheets

print(stylesheets.combined)

def make_cyto_elements(network):
    """Takes a BioNetwork and outputs Cytoscape elements that can be visualized with Dash. Also creates selector
    classes according to the attributes and the layout coordinates."""

    # Get nodes of the largest component and generate sub-network
    main_component_nodes = [component for component in nx.connected_components(temp_network)][0]
    network = network.subgraph(main_component_nodes)

    # Convert nx network to cytoscape JSON
    json_elements = nx.readwrite.json_graph.cytoscape_data(network)['elements']
    layout = nx.spring_layout(network, k=10 / sqrt(len(network)), scale=1000, center=[500, 500])
    # layout = nx.spring_layout(network)

    nodes = json_elements['nodes']
    for node in nodes:
        node['data']['label'] = node['data']['shortName']  # Use short name as node label

        pos = layout[node['data']['id']]  # Extract positions from layout dict
        node['position'] = {'x': pos[0], 'y': pos[1]}
        #node['position'] = {'x': pos[0] * 1000 + 500, 'y': pos[1] * 1000 + 500}

    elements = nodes + json_elements['edges']

    return elements


# temp_bionetwork = ng.CombinedNetwork(h.get_genes('/home/javier/Documents/Corrie/TnSeq_CB_Manuscript-master/RNASeq/results/mediarpmi.treatmentazm.csv'),
#                                      pd.read_csv('/home/javier/Documents/Corrie/TnSeq_CB_Manuscript-master/RNASeq/results/mediarpmi.treatmentazm.csv'),
#                                      h.get_genes('/home/javier/Documents/Corrie/TnSeq_CB_Manuscript-master/TnSeq/in-vitro/essential/finalEss_noTrue_RPMI_AZMvsunt_20190715.csv'),
#                                      'PAO1',
#                                      0,
#                                      3,
#                                      False)
#
# temp_network = temp_bionetwork.network
# nx.write_graphml(temp_network, 'corries_AZM_combined.graphml')
temp_network = nx.read_graphml('corries_AZM_combined.graphml')
cyto_elements = make_cyto_elements(temp_network)


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = html.Div(
    dbc.Container(
        [
            dbc.Row(
                [
                    dbc.Col(
                        html.Div(
                            [
                                dbc.Label('Color Mapping'),
                                dbc.RadioItems(
                                    options=[
                                        {'label': 'Significance Source', 'value': 'ss'},
                                        {'label': 'Differential Expression', 'value': 'de'}
                                    ],
                                    value='ss',
                                    id='color-map',
                                    ),
                                html.H5('Select nodes:'),
                                html.Details(
                                    [
                                        html.Summary('By source of interest: '),
                                        dbc.RadioItems(
                                            id='source-selection',
                                            options=[
                                                {'label': 'RNASeq', 'value': 'transcriptional'},
                                                {'label': 'TnSeq', 'value': 'phenotypical'},
                                                {'label': 'Both', 'value': 'both'}
                                            ],
                                            value='transcriptional'
                                        )
                                    ],
                                ),
                                html.Br(),
                                dbc.Button(id='reset-view',
                                           children='Reset View')
                            ],
                            style={'margin': '2vh'}
                        ),
                        width=3,
                        style={'backgroundColor': '#7FDBFF', 'padding': '5px'}
                    ),
                    dbc.Col(
                        [
                            cyto.Cytoscape(
                                id='cyto-network',
                                style={'width': '100%',
                                       'height': '80vh',
                                       'line-color': 'red'
                                       },
                                elements=cyto_elements,
                                zoom=1,
                                maxZoom=5,
                                minZoom=0.3
                            ),
                            html.Div(
                                [
                                    html.H5('Selected Node Details')
                                ],
                                style={'height': '20vh'}
                            )
                        ],
                        width=8
                    ),
                ]
            )
        ],
        fluid=True,
    ),
)


@app.callback(
    [Output('cyto-network', 'zoom'),
     Output('cyto-network', 'elements')],
    [Input('reset-view', 'n_clicks')],
)
def reset_layout(n_clicks):
    return [1, cyto_elements]


@app.callback(
    [Output('cyto-network', 'stylesheet')],
    [Input('color-map', 'value')]
)
def change_color_map(value):
    if value == 'ss':
        return [stylesheets.combined]
    else:
        return [stylesheets.fold_change]


if __name__ == "__main__":
    app.run_server(debug=True)