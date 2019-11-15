import re

import dash
from dash.dependencies import Output, Input, State
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import dash_html_components as html
import networkx as nx


def make_cyto_elements(network):
    """Takes a networkx network and outputs cyto elements that can be visualized with Dash."""
    cytoscape_elements = []

    for node in network.nodes():
        node_info = network.nodes(data=True)[node]
        cyto_node = {'data': {
            'id': node,
            'label': node_info['shortName']},
            'classes': '{} {}'.format(node_info['type'], node_info['localization'])
        }
        cytoscape_elements.append(cyto_node)

    for edge in nx.generate_edgelist(network):
        edge = re.sub(r'\s\{.*?\}', '', edge)  # Only keep node names
        node1, node2 = edge.split(' ')
        cyto_edge = {'data': {'source': node1,
                              'target': node2,
                              'label': '{} to {}'.format(node1, node2)}}
        cytoscape_elements.append(cyto_edge)

    return cytoscape_elements


temp_network = nx.read_graphml('/home/javier/PycharmProjects/PaIntDB/temp_data/azm_vs_control_network.graphml')
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
                                #dbc.Label('By mniguebo'),
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
                                layout={'name': 'grid'},
                                style={'width': '100%',
                                       'height': '80vh',
                                       'line-color': 'red',
                                       'backgroundColor': '#7FDBFF',
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


if __name__ == "__main__":
    app.run_server(debug=True)