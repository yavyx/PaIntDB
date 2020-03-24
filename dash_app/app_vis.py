from datetime import datetime
from math import sqrt
from itertools import chain
import os

import dash
from dash.dependencies import Output, Input, State
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_cytoscape as cyto
import dash_html_components as html
import dash_table
import flask
import networkx as nx
from networkx.algorithms import approximation
from networkx.utils import pairwise
import pandas as pd

import dash_app.cytoscape_stylesheets as stylesheets
import go_enrichment.go_enrichment as goe


def make_network_df(network):
    """Takes a network and outputs a DataFrame of every node with its attributes."""
    data = [i[1] for i in network.nodes(data=True)]  # Get node attributes
    index = [i[0] for i in network.nodes(data=True)]  # Use node ids as index
    df = pd.DataFrame(data, index)
    # format fields
    df['log2FoldChange'] = df['log2FoldChange'].round(2)
    df['padj'] = df['padj'].map('{:.3g}'.format)
    df['padj'] = df['padj'].astype(float)
    df['regulation'] = ['up' if change > 0 else 'down' for change in df['log2FoldChange']]
    if 'significanceSource' in df.columns:
        df['regulation'] = [None if sig == 'TnSeq' else reg for sig, reg in zip(df['significanceSource'], df['regulation'])]
    return df


def make_cyto_elements(network):
    """Takes a network and outputs Cytoscape elements that can be visualized with Dash. Also creates selector
    classes according to the attributes and the layout coordinates."""

    # Get nodes of the largest component and generate sub-network
    main_component_nodes = [component for component in nx.connected_components(temp_network)][0]
    network = network.subgraph(main_component_nodes)
    nx.set_node_attributes(network, dict(network.degree()), 'degree')

    # Convert nx network to cytoscape JSON
    json_elements = nx.readwrite.json_graph.cytoscape_data(network)['elements']

    # Make layout (much faster than default Cytoscape layouts)
    layout = nx.spring_layout(network, k=2 / sqrt(len(network)), scale=1000)
    nodes = json_elements['nodes']
    for node in nodes:
        node['data']['label'] = node['data']['shortName']  # Use short name as node label

        pos = layout[node['data']['id']]  # Extract positions from layout dict
        node['position'] = {'x': pos[0], 'y': pos[1]}  # Set positions

    edges = json_elements['edges']
    elements = nodes + edges

    return elements, nodes, edges, network


network_path = os.path.join('temp_data', 'Biofilm_vs_planktonic_PA14_combined_network.graphml')
temp_network = nx.read_graphml(network_path)
all_nodes = len(temp_network.nodes())

# Extract and use just main component
nodes = [component for component in nx.connected_components(temp_network)][0]
temp_network = temp_network.subgraph(nodes)
new_nodes = len(temp_network.nodes())
print('Lost {} nodes'.format(all_nodes - new_nodes))

print('Calculating metric closure')
t_start = datetime.now()
metric_closure = approximation.metric_closure(temp_network)
print(datetime.now() - t_start)

strain = 'PA14'

cyto_elements, cyto_nodes, cyto_edges, network_main_comp = make_cyto_elements(temp_network)
network_df = make_network_df(network_main_comp)

enrichment_results, goea_results = goe.run_go_enrichment(strain, list(network_df.index))


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = html.Div(
    [html.Div(
        style={
            'width': '23vw',
            'backgroundColor': '#7FDBFF',
            'padding': '10px',
            'display': 'inline-block',
            'height': '100vh',
            'vertical-align': 'top'
        },
        children=[
            html.H5('Color Mapping'),
            dbc.RadioItems(
                options=[
                    {'label': 'Experiment', 'value': 'ss'},
                    {'label': 'Differential Expression', 'value': 'de'}
                ],
                value='ss',
                id='color-map',
            ),
            html.Div(style={'padding': '10px'},
                     children=html.Img(
                         id='legend',
                         width=100)
                     ),
            html.Hr(),
            html.H5('Select nodes'),
            html.Details(
                [
                    html.Summary('By name'),
                    dcc.Dropdown(
                        id='name-selection',
                        # Search by both short name and locus tag (index)
                        options=[{'label': term, 'value': term}
                                 for term in network_df['shortName']] +
                                [{'label': term, 'value': term}
                                 for term in network_df.index],
                        multi=True,
                        optionHeight=50
                    )
                ],
            ),
            html.Br(),
            html.Details(
                children=[
                    html.Summary('By experiment '),
                    dbc.Checklist(
                        id='source-selection',
                        options=[
                            {'label': 'RNASeq', 'value': 'RNASeq'},
                            {'label': 'TnSeq', 'value': 'TnSeq'}
                        ],
                        value=[],
                    )
                ],
            ),
            html.Br(),
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
            html.Details(
                [
                    html.Summary('By differential expression'),
                    dbc.Checklist(
                        id='diff-exp-selection',
                        options=[
                            {'label': 'Up-Regulated', 'value': 'up'},
                            {'label': 'Down-Regulated', 'value': 'down'}
                        ],
                        value=[],
                    )
                ],
            ),
            html.Br(),
            html.Details(
                [
                    html.Summary('By enriched GO term'),
                    dcc.Dropdown(
                        id='enrichment-selection',
                        options=[{'label': term, 'value': term}
                                 for term in enrichment_results['name']],
                        multi=True,
                        optionHeight=50
                    )
                ],
            ),
            html.Br(),
            html.Div(
                [
                    html.P(id='num-selected-nodes'),
                    html.Button('Make Sub-Network',
                                id='make-subnetwork',
                                style={'display': 'inline-block'}),
                    html.Button(id='download-table-button',
                                children=html.A('Download Table',
                                                id='table-download-link'
                                                )
                                ),
                    html.Div(id='node-details-download', style={'display': 'none'}),
                    html.Div(id='hidden-div', style={'display': 'none'})
                ]
            )
        ],
    ),
        html.Div(
            style={
                'display': 'inline-block',
                'width': '74vw'
            },
            children=dbc.Container(
                fluid=True,
                children=[
                    dbc.Row(
                        [
                            dbc.Col(
                                id='main-view',
                                width=6,
                                children=[
                                    html.H5('Full Network View'),
                                    cyto.Cytoscape(
                                        id='cytoscape',
                                        style={
                                            'width': '100%',
                                            'height': '65vh'
                                        },
                                        elements=cyto_elements,
                                        maxZoom=5,
                                        minZoom=0.3,
                                        zoom=1,
                                        layout={'name': 'preset'},
                                    )
                                ],
                            ),
                            dbc.Col(
                                id='sub-view',
                                width=6,
                                children=[
                                    html.H5('Sub-network View'),
                                    cyto.Cytoscape(
                                        id='sub_cytoscape',
                                        style={
                                            'width': '100%',
                                            'height': '65vh'
                                        },
                                        elements=cyto_elements,
                                        maxZoom=5,
                                        minZoom=0.3,
                                        zoom=1,
                                        layout={'name': 'preset'},
                                        boxSelectionEnabled=True
                                    )
                                ]
                            )
                        ]
                    ),
                    dbc.Row(
                        dbc.Col(
                            html.Div(
                                [
                                    html.H5('Selected Node Details'),
                                    html.Div(id='node-details-table')
                                ],
                                style={'height': '30vh'}
                            )
                        )
                    )
                ]
            ))
    ]
)


@app.callback(
    [Output('cytoscape', 'stylesheet'),
     Output('sub_cytoscape', 'stylesheet'),
     Output('legend', 'src')],
    [Input('color-map', 'value')]
)
def change_color_map(value):
    if value == 'ss':
        source = app.get_asset_url('sig_source_legend.svg')
        return stylesheets.combined, stylesheets.combined, source
    else:
        source = app.get_asset_url('de_legend.svg')
        return stylesheets.fold_change, stylesheets.fold_change, source


@app.callback(
    [Output('cytoscape', 'elements'),
     Output('num-selected-nodes', 'children')],
    [Input('name-selection', 'value'),
     Input('source-selection', 'value'),
     Input('location-selection', 'value'),
     Input('diff-exp-selection', 'value'),
     Input('enrichment-selection', 'value')]
)
def select_nodes(short_name, significance_source, location, regulation, enriched_terms):
    """Select nodes according to significance source."""
    nodes = cyto_nodes  # Make new variable to avoid modifying the global cyto_nodes

    query = []  # Query to filter nodes

    # Add queries depending on user menu selections
    if short_name:
        query.append('shortName in @short_name | index in @short_name')

    if significance_source:
        significance_source.append('both')
        query.append('significanceSource in @significance_source')

    if location:
        query.append('localization in @location')

    if regulation:
        query.append('regulation in @regulation')

    if enriched_terms:
        # Get genes associated with selected GO term(s)
        genes_in_term = enrichment_results.loc[enrichment_results['name'].isin(enriched_terms), 'study_items']
        total_genes = [gene for term in genes_in_term for gene in term.split(', ')]
        if strain == 'PA14':
            total_genes = goe.map_pao1_genes(total_genes)
        query.append('index in @total_genes')

    query_str = ' & '.join(query)
    print(query_str)

    # Use query to select nodes
    if query_str:
        queried_nodes = network_df.query(query_str).index.tolist()
        selected_msg = 'Selected {} out of {} nodes'.format(len(queried_nodes), len(cyto_nodes))
    else:
        queried_nodes = nodes
        selected_msg = ''

    for node in nodes:
        node['selected'] = True if node['data']['id'] in queried_nodes else False

    return [nodes + cyto_edges, selected_msg]


@app.callback(
    [Output('node-details-table', 'children'),
     Output('node-details-download', 'children')],
    [Input('cytoscape', 'selectedNodeData'),
     Input('sub_cytoscape', 'selectedNodeData')]
)
def show_node_details(node_data, sub_node_data):
    """Filters the network DataFrame with the user-selected nodes and returns a DataTable."""
    if node_data:
        # Columns to display
        cols = ['shortName', 'description', 'log2FoldChange', 'padj']
        node_ids = [node['label'] for node in node_data]
        filtered_df = (network_df.loc[network_df.shortName.isin(node_ids), cols]
                       .reset_index()
                       .rename(columns={'index': 'Locus Tag',
                                        'shortName': 'Short Name',
                                        'description': 'Descripton',
                                        'log2FoldChange': 'Log2 Fold Change',
                                        'padj': 'Adjusted p-value'})
                       )

        if sub_node_data:
            node_ids = [node['label'] for node in sub_node_data]
            filtered_df = (network_df.loc[network_df.shortName.isin(node_ids), cols]
                           .reset_index()
                           .rename(columns={'index': 'Locus Tag',
                                            'shortName': 'Short Name',
                                            'description': 'Descripton',
                                            'log2FoldChange': 'Log2 Fold Change',
                                            'padj': 'Adjusted p-value'})
                           )

        nodes_table = dash_table.DataTable(
            data=filtered_df.to_dict('records'),
            columns=[{"name": i, "id": i} for i in filtered_df.columns],
            fixed_rows={'headers': True, 'data': 0},
            style_table={
                'maxHeight': '25vh',
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
        return nodes_table, filtered_df.to_json()
    else:
        return None, None


@app.callback(
    Output('table-download-link', 'href'),
    [Input('download-table-button', 'n_clicks')],
    [State('node-details-download', 'children')]
)
def update_download_link(n_clicks, node_details):
    if n_clicks:
        rel_filename = os.path.join('downloads', 'node_table.csv')
        abs_filename = os.path.join(os.getcwd(), rel_filename)
        pd.read_json(node_details).to_csv(abs_filename)
        return '/{}'.format(rel_filename)


@app.server.route('/downloads/<path:path>')
def download_csv(path):
    root_dir = os.getcwd()
    return flask.send_from_directory(
        os.path.join(root_dir, 'downloads'),
        path,
        cache_timeout=-1  # Prevent browser from caching previous file
    )


@app.callback(
    Output('hidden-div', 'children'),
    [Input('sub_cytoscape', 'selectedNodeData')],
    [State('sub_cytoscape', 'selectedNodeData')]
)
def print_nodes(n_clicks, node_data):
    if n_clicks:
        print('miguebo  ', [node['id'] for node in node_data])
    return None


@app.callback(
    Output('sub_cytoscape', 'elements'),
    [Input('make-subnetwork', 'n_clicks')],
    [State('cytoscape', 'selectedNodeData')]
)
def make_subnetwork(n_clicks, node_data):
    if n_clicks and node_data:
        node_ids = [node['id'] for node in node_data]  # Get selected node ids.
        H = metric_closure.subgraph(node_ids)
        mst_edges = nx.minimum_spanning_edges(H, weight='distance', data=True)
        edges = chain.from_iterable(pairwise(d['path']) for u, v, d in mst_edges)
        T = temp_network.edge_subgraph(edges)
        sub_network = T

        cyto_sub_network, sub_cyto_nodes, sub_cyto_edges, subnetwork = make_cyto_elements(sub_network)
        return cyto_sub_network
    else:
        return []


# Might Change this back later (for now, running enrichment when app starts)
# @app.callback(
#     Output('enrichment-selection', 'options'),
#     [Input('enrichment-results', 'children')]
# )
# def create_enrichment_options(results):
#     results_df = pd.read_json(results)
#     options = [{'label': term, 'value': term} for term in results_df['name']]
#     return options


# Might Change this back later (for now, running enrichment when app starts)
# @app.callback(
#     Output('enrichment-results', 'children'),
#     [Input('run-enrichment', 'n_clicks')]
# )
# def run_enrichment(n_clicks):
#     if n_clicks:
#         enrichment_results, goea_results = goe.run_go_enrichment('PAO1', list(network_df.index))
#         print(enrichment_results.head())
#         table = enrichment_results.to_json()
#     else:
#         table = ""
#     return table


if __name__ == "__main__":
    app.run_server(debug=True)
