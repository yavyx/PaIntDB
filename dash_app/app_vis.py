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
    layout = nx.spring_layout(network, k=5 / sqrt(len(network)), scale=1000, center=[450, 450])
    nodes = json_elements['nodes']
    for node in nodes:
        node['data']['label'] = node['data']['shortName']  # Use short name as node label

        pos = layout[node['data']['id']]  # Extract positions from layout dict
        node['position'] = {'x': pos[0], 'y': pos[1]}

    edges = json_elements['edges']
    elements = nodes + edges

    return elements, nodes, edges


# temp_network = nx.read_graphml('temp_data/Biofilm_DJK6_vs_Biofilm_network.graphml')
# temp_network = nx.read_graphml('temp_data/Biofilm_DJK5_vs_Biofilm_network (1).graphml')
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
                                            value=['RNASeq', 'TnSeq', 'both'],
                                        )
                                    ],
                                ),
                                html.Details(
                                    [
                                        html.Summary('By localization '),
                                        dbc.Checklist(
                                            id='location-selection',
                                            options=[
                                                dict(label=location, value=location)
                                                for location in network_df['localization'].unique()
                                            ],
                                            value=[],
                                        )
                                    ],
                                ),
                                html.Br(),
                                html.Div(
                                    html.Button('Make Selection',
                                                id='make-selection',
                                                style={'display': 'inline-block'}),
                                    style={'text-align': 'center'}),
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
                        width=10
                    ),
                ]
            )
        ],
        fluid=True,
    ),
)


@app.callback(
    Output('cytoscape', 'stylesheet'),
    [Input('color-map', 'value')]
)
def change_color_map(value):
    if value == 'ss':
        return stylesheets.combined
    else:
        return stylesheets.fold_change


@app.callback(
    Output('cytoscape', 'elements'),
    [Input('make-selection', 'n_clicks')],
    [State('source-selection', 'value'),
     State('location-selection', 'value')]
)
def select_nodes(n_clicks, significance_source, location):
    """Select nodes according to significance source."""
    nodes = cyto_nodes  # Make new variable to avoid modifying the global cyto_nodes

    print('significance', significance_source)
    print('location', location)

    query = []  # query to filter nodes
    print(query)
    if significance_source:
        query.append('significanceSource in @significance_source')
    if location:
        query.append('localization in @location')
    print(query)
    query_str = ' & '.join(query)
    if query_str:
        queried_nodes = list(network_df.query(query_str).index)
    else:
        queried_nodes = nodes
    print('# nodes', len(queried_nodes))

    for node in nodes:
        if node['data']['id'] in queried_nodes:
            node['selected'] = True
        else:
            node['selected'] = False

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
        filtered_df = (network_df.loc[network_df.shortName.isin(node_ids), cols]
                       .reset_index()
                       .rename(columns={'index': 'Locus Tag'})
                       )

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
