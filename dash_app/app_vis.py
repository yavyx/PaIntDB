import re

import dash
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import dash_html_components as html
import networkx as nx


def make_cyto_elements(graph):
    """Takes a networkx network and outputs cyto elements that can be visualized with Dash."""
    cytoscape_elements = []

    for node in graph.nodes():
        node_info = graph.nodes(data=True)[node]
        cyto_node = {'data': {
            'id': node,
            'label': node_info['shortName']}
        }
        cytoscape_elements.append(cyto_node)

    for edge in nx.generate_edgelist(graph):
        edge = re.sub(r'\s\{.*?\}', '', edge)  # keep only node names
        node1, node2 = edge.split(' ')
        cyto_edge = {'data': {'source': node1,
                              'target': node2,
                              'label': '{} to {}'.format(node1, node2)}}
        cytoscape_elements.append(cyto_edge)

    return cytoscape_elements


temp_network = nx.read_graphml('/home/javier/PycharmProjects/PaIntDB/temp_data/azm_vs_control_network.graphml')
#temp_network = nx.powerlaw_cluster_graph(1000, 3, 0.5)
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
                                html.H4('Select nodes:', style={'margin': '30px', 'backgroundColor': '#111111'}),
                                html.Br(),
                                dbc.Label('By source of interest: '),
                                dbc.RadioItems(
                                    id='source-selection',
                                    options=[
                                        {'label': 'RNASeq', 'value': 'transcriptional'},
                                        {'label': 'TnSeq', 'value': 'phenotypical'},
                                        {'label': 'Both', 'value': 'both'}
                                    ],
                                    value='transcriptional'
                                ),
                                html.Br(),
                                dbc.Label('By mniguebo')

                            ],
                            style={'backgroundColor': '#7FDBFF', 'margin': '0px'}
                        ),
                        width=2,
                    ),
                    dbc.Col(
                        [
                            cyto.Cytoscape(
                                id='network',
                                layout={'name': 'grid'},
                                style={'width': '100%', 'height': '100vh', 'line-color': 'red',
                                       'backgroundColor': '#7FDBFF'},
                                elements=cyto_elements
                            )
                        ],
                        width=7
                    ),
                    dbc.Col([html.H4("Subcluster Details"),
                             html.Div(id='node_genomes')],
                            width=3)
                ]
            )
        ],
        fluid=True,
        style={'height': '50vh'}
    ),
    style={'height': '50vh'}
)

if __name__ == "__main__":
    app.run_server(debug=True)