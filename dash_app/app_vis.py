from math import sqrt

import dash
from dash.dependencies import Output, Input, State
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_cytoscape as cyto
import dash_html_components as html
import dash_table
import networkx as nx
import pandas as pd

import dash_app.cytoscape_stylesheets as stylesheets


def make_network_df(network):
    """Takes a network and outputs a DataFrame of every node with its attributes."""
    data = [i[1] for i in network.nodes(data=True)]  # Get node attributes
    index = [i[0] for i in network.nodes(data=True)]  # Use node ids as index
    df = pd.DataFrame(data, index)

    # format fields
    df['log2FoldChange'] = df['log2FoldChange'].round(2)
    df['padj'] = df['padj'].map('{:.3g}'.format)
    df['padj'] = df['padj'].astype(float)
    return df


def make_cyto_elements(network):
    """Takes a network and outputs Cytoscape elements that can be visualized with Dash. Also creates selector
    classes according to the attributes and the layout coordinates."""

    # Get nodes of the largest component and generate sub-network
    main_component_nodes = [component for component in nx.connected_components(temp_network)][0]
    network = network.subgraph(main_component_nodes)

    # Convert nx network to cytoscape JSON
    json_elements = nx.readwrite.json_graph.cytoscape_data(network)['elements']

    # Make layout (much faster than default Cytoscape layouts)
    layout = nx.spring_layout(network, k=10 / sqrt(len(network)), scale=1000, center=[450, 450])
    nodes = json_elements['nodes']
    for node in nodes:
        node['data']['label'] = node['data']['shortName']  # Use short name as node label

        pos = layout[node['data']['id']]  # Extract positions from layout dict
        node['position'] = {'x': pos[0], 'y': pos[1]}

    edges = json_elements['edges']
    elements = nodes + edges

    return elements, nodes, edges


temp_network = nx.read_graphml('corries_AZM_combined.graphml')
network_df = make_network_df(temp_network)
cyto_elements, cyto_nodes, cyto_edges = make_cyto_elements(temp_network)


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = html.Div(
    dbc.Container(
        [
            dbc.Row(
                [
                    dbc.Col(
                        html.Div(
                            [
                                html.H5('Color Mapping'),
                                dbc.RadioItems(
                                    options=[
                                        {'label': 'Significance Source', 'value': 'ss'},
                                        {'label': 'Differential Expression', 'value': 'de'}
                                    ],
                                    value='ss',
                                    id='color-map',
                                    ),
                                html.Hr(),
                                html.H5('Select nodes'),
                                html.Details(
                                    [
                                        html.Summary('By source of interest '),
                                        dbc.Checklist(
                                            id='source-selection',
                                            options=[
                                                {'label': 'RNASeq', 'value': 'RNASeq'},
                                                {'label': 'TnSeq', 'value': 'TnSeq'},
                                                {'label': 'Both', 'value': 'both'}
                                            ],
                                            value=['RNASeq', 'TnSeq', 'both']
                                        )
                                    ],
                                ),
                                html.Br(),
                            ],
                            style={'margin': '2vh'}
                        ),
                        width=2,
                        style={'backgroundColor': '#7FDBFF', 'padding': '5px'}
                    ),
                    dbc.Col(
                        [
                            cyto.Cytoscape(
                                id='cytoscape',
                                style={
                                    'width': '100%',
                                    'height': '75vh'
                                },
                                elements=cyto_elements,
                                maxZoom=5,
                                minZoom=0.3,
                                zoom=1,
                                layout={'name': 'preset'},
                            ),
                            html.Div(
                                [
                                    html.H5('Selected Node Details'),
                                    html.Div(id='node-details')
                                ],
                                style={'height': '25vh'}
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
    Output('cytoscape', 'stylesheet'),
    [Input('color-map', 'value'),
     Input('source-selection', 'value')]
)
def change_color_map(value, elements):
    if value == 'ss':
        return stylesheets.combined
    else:
        return stylesheets.fold_change


@app.callback(
    Output('cytoscape', 'elements'),
     #Output('cytoscape', 'stylesheet')],
    [Input('source-selection', 'value')],
    #[State('cytoscape', 'stylesheet')]
)
def select_by_source(value):
    """Select nodes according to significance source."""
    nodes = cyto_nodes  # Make new variable to avoid modifying the global cyto_nodes
    if value:
        for node in nodes:
            if node['data']['significanceSource'] in value:
                node['data']['selected'] = True
            else:
                node['data']['selected'] = False
        print(value)
        print(nodes[-1]['data'], '\n')
        return nodes + cyto_edges


@app.callback(
    Output('node-details', 'children'),
    [Input('cytoscape', 'selectedNodeData')]
)
def show_node_details(node_data):
    """Filters the network DataFrame with the user-selected nodes and returns a DataTable."""
    if node_data:
        # Columns to display
        cols = ['shortName', 'description', 'log2FoldChange', 'padj', 'ncbi', 'uniprotkb']
        node_ids = [node['label'] for node in node_data]
        filtered_df = network_df.loc[network_df.shortName.isin(node_ids), cols]

        table = dash_table.DataTable(
            data=filtered_df.to_dict('records'),
            columns=[{"name": i, "id": i} for i in filtered_df.columns],
            fixed_rows={'headers': True, 'data': 0},
            style_table={
                'maxHeight': '20vh',
                'overflowY': 'auto'
            },
            style_data={'whiteSpace': 'normal',
                        'height': 'auto'
                        },
            style_cell={
                'height': 'auto',
                # all three widths are needed
                'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
                'whiteSpace': 'normal',
                'textAlign': 'left'
            }
        )
        return table


if __name__ == "__main__":
    app.run_server(debug=True)
