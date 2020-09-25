from datetime import datetime
from math import sqrt
from itertools import chain
import json
import os

import dash
from dash.dependencies import Output, Input, State, ALL
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_cytoscape as cyto
import dash_html_components as html
import dash_table
import flask
import networkx as nx
from networkx.utils import pairwise
import pandas as pd

from dash_app import app
import dash_app.vis_stylesheets as stylesheets
from dash_app.app import app  # Loads app variable from app script
import go_enrichment.go_enrichment as goe

# For development  # TODO: delete
pd.set_option('display.max_columns', None)


def make_cyto_elements(network, k, scale):
    """Takes a network and outputs Cytoscape elements that can be visualized with Dash. Also creates selector
    classes according to the attributes and the layout coordinates."""

    # Get nodes of the largest component and generate sub-network
    main_component_nodes = max(nx.connected_components(network), key=len)
    network = network.subgraph(main_component_nodes)
    nx.set_node_attributes(network, dict(network.degree()), 'degree')

    # Convert nx network to cytoscape JSON
    json_elements = nx.readwrite.json_graph.cytoscape_data(network)['elements']

    # Make layout (much faster than default Cytoscape layouts)
    layout = nx.spring_layout(network, k=k / sqrt(len(network)), scale=scale)
    nodes = json_elements['nodes']
    for node in nodes:
        node['data']['label'] = node['data']['shortName']  # Use short name as node label

        pos = layout[node['data']['id']]  # Extract positions from layout dict
        node['position'] = {'x': pos[0], 'y': pos[1]}  # Set positions

    edges = json_elements['edges']
    elements = nodes + edges

    return elements, nodes, edges, network


def make_vis_layout(network_df, enrichment_results, cyto_network, network_params):
    network_params = json.loads(network_params)
    # This filter is only used with RNASeq/Combined networks
    regulation_filter = html.Details(
        id='diff-exp-details',
        open=True,
        children=[
            html.Summary('By differential expression'),
            dbc.Checklist(
                id={
                    'type': 'filter',
                    'index': 3
                },
                options=[
                    {'label': 'Up-Regulated', 'value': 'up'},
                    {'label': 'Down-Regulated', 'value': 'down'}
                ],
                value=[],
            )
        ],
    )
    # This filter is only used with Combined networks
    source_filter = html.Details(
        id='source-details',
        open=True,
        children=[
            html.Summary('By experiment '),
            dbc.Checklist(
                id={
                    'type': 'filter',
                    'index': 4
                },
                options=[
                    {'label': 'RNASeq', 'value': 'RNASeq'},
                    {'label': 'TnSeq', 'value': 'TnSeq'}
                ],
                value=[],
            )
        ],
    )
    # These filters are used for all networks
    sidebar_filters = [
        html.H5('Select nodes'),
        html.Details(
            open=True,
            children=[
                html.Summary('By name'),
                dcc.Dropdown(
                    id={
                        'type': 'filter',
                        'index': 0
                    },
                    # Allow search by both short name and locus tag (index)
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
            open=True,
            children=[
                html.Summary('By localization '),
                dcc.Dropdown(
                    id={
                        'type': 'filter',
                        'index': 1
                    },
                    options=[
                        dict(label=location, value=location)
                        for location in network_df['localization'].unique()
                    ],
                    multi=True,
                )
            ],
        ),
        html.Br(),
        html.Details(
            open=True,
            children=[
                html.Summary('By enriched GO term'),
                dcc.Dropdown(
                    id={
                        'type': 'filter',
                        'index': 2
                    },
                    # Use enriched GO terms as options.
                    options=[{'label': term, 'value': term}
                             for term in enrichment_results['name']],
                    multi=True,
                    optionHeight=50
                )
            ],
        ),
        html.Br()
    ]
    if network_params['type'] == 'rna_seq' or network_params['type'] == 'combined':
        sidebar_filters.extend([
            regulation_filter,
            html.Br()
        ])
        if network_params['type'] == 'combined':
            sidebar_filters.extend([
                source_filter,
            ])

    return html.Div(
        [
            html.Div(
                style={
                    'width': '23vw',
                    'backgroundColor': '#7FDBFF',
                    'padding': '10px',
                    'display': 'inline-block',
                    'height': '95vh',
                    'vertical-align': 'top'
                },
                children=[
                    html.Div(id='color-map-div',
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
                                          )
                             ]
                             ),
                    html.Hr(),
                    html.Div(
                        id='select-nodes',
                        children=sidebar_filters
                    ),
                    html.Br(),
                    html.Div(
                        [
                            html.P(id='num-selected-nodes'),
                            dbc.Button('Make Sub-Network',
                                       id='make-subnetwork',
                                       style={'display': 'inline-block'}),
                            dbc.Button('Reset Network',
                                       id='reset-network',
                                       style={'display': 'inline-block'}),
                            dbc.Button('Download Table',
                                       id='download-table-button',
                                       ),
                            dbc.Button('Download Subnetwork',
                                       id='download-subnetwork'),
                            # Hidden Div to store node details table download file
                            html.Div(id='node-details-download', style={'display': 'none'}),
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
                                    children=[
                                        cyto.Cytoscape(
                                            id='main-view',
                                            style={
                                                'width': '100%',
                                                'height': '60vh',
                                                # 'position': 'absolute'
                                            },
                                            stylesheet=stylesheets.default,
                                            maxZoom=5,
                                            minZoom=0.3,
                                            zoom=1,
                                            layout={'name': 'preset'},
                                            elements=cyto_network['elements']
                                        ),
                                        cyto.Cytoscape(
                                            id='sub-view',
                                            style={
                                                'width': '100%',
                                                'height': '60vh',
                                                # 'display': 'none'
                                                # 'position': 'absolute'
                                            },
                                            stylesheet=stylesheets.default,
                                            maxZoom=5,
                                            minZoom=0.3,
                                            zoom=1,
                                            layout={'name': 'preset'},
                                            boxSelectionEnabled=True
                                        )

                                    ],
                                ),
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
                )
            )
        ]
    )


@app.callback(
    [Output('main-view', 'stylesheet'),
     Output('sub-view', 'stylesheet'),
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
    [Output('main-view', 'elements'),
     Output('num-selected-nodes', 'children')],
    [Input({'type': 'filter', 'index': ALL}, 'value')],  # Pattern-matching all callbacks with filter type
    [State('node-details-df', 'children'),
     State('enrichment-results', 'children'),
     State('network-parameters', 'children'),
     State('cyto-network', 'children')]
)
def select_nodes(values, node_details, enrichment_results, network_params, cyto_network):
    """Select nodes according to user selected filters."""
    enrichment_results = pd.read_json(enrichment_results)
    cyto_network = json.loads(cyto_network)
    nodes = cyto_network['nodes']
    edges = cyto_network['edges']
    network_df = pd.read_json(node_details)
    network_params = json.loads(network_params)
    strain = network_params['strain']
    network_type = network_params['type']
    query = []  # Query to filter nodes

    # Filter output values
    short_name = values[0]
    location = values[1]
    enriched_terms = values[2]
    if network_type == 'rna_seq' or network_params['type'] == 'combined':
        regulation = values[3]
        if network_params['type'] == 'combined':
            significance_source = values[4]
    else:
        regulation, significance_source = [], []

    # Add queries depending on user filter selections
    if short_name:
        query.append('shortName in @short_name | index in @short_name')
    if location:
        query.append('localization in @location')
    if enriched_terms:
        # Get genes associated with selected GO term(s)
        genes_in_term = enrichment_results.loc[enrichment_results['name'].isin(enriched_terms), 'study_items']
        total_genes = [gene for term in genes_in_term for gene in term.split(', ')]  # Unnest genes in term
        if strain == 'PA14':
            total_genes = goe.map_pao1_genes(total_genes)
        query.append('index in @total_genes')
    if significance_source:
        significance_source.append('both')
        query.append('significanceSource in @significance_source')
    if regulation:
        query.append('regulation in @regulation')
    query_str = ' & '.join(query)
    print(query_str)

    # Use query to select nodes
    if query_str:
        queried_nodes = network_df.query(query_str).index.tolist()
        selected_msg = 'Selected {} out of {} nodes'.format(len(queried_nodes), len(nodes))
    else:
        queried_nodes = nodes
        selected_msg = ''

    for node in nodes:
        node['selected'] = True if node['data']['id'] in queried_nodes else False

    return nodes + edges, selected_msg


@app.callback(
    [Output('node-details-table', 'children'),
     Output('node-details-download', 'children')],
    [Input('main-view', 'selectedNodeData'),
     Input('sub-view', 'selectedNodeData')],
    [State('node-details-df', 'children'),
     State('network-parameters', 'children')]
)
def show_node_details(node_data, sub_node_data, node_details, network_params):
    """Filters the network DataFrame with the user-selected nodes and returns a DataTable."""
    if node_data:
        # Columns to display
        print('miguebo')
        cols = ['shortName', 'description']
        network_params = json.loads(network_params)
        if network_params['type'] == 'rna_seq' or network_params['type'] == 'combined':
            cols.extend(['log2FoldChange', 'padj'])
        node_ids = [node['label'] for node in node_data]
        print(node_ids)
        network_df = pd.read_json(node_details)
        filtered_df = (network_df.loc[network_df.shortName.isin(node_ids), cols]
                       .reset_index()
                       .rename(columns={'index': 'Locus Tag',
                                        'shortName': 'Short Name',
                                        'description': 'Description',
                                        'log2FoldChange': 'Log2 Fold Change',
                                        'padj': 'Adjusted p-value'
                                        }
                               )
                       )
        if sub_node_data:
            node_ids = [node['label'] for node in sub_node_data]
            filtered_df = (network_df.loc[network_df.shortName.isin(node_ids), cols]
                           .reset_index()
                           .rename(columns={'index': 'Locus Tag',
                                            'shortName': 'Short Name',
                                            'description': 'Descripton',
                                            'log2FoldChange': 'Log2 Fold Change',
                                            'padj': 'Adjusted p-value'
                                            }
                                   )
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
   [Output('sub-view', 'elements'),
    Output('main-view', 'style'),
    Output('sub-view', 'style')],
   [Input('make-subnetwork', 'n_clicks'),
    Input('reset-network', 'n_clicks')],
   [State('main-view', 'selectedNodeData'),
    State('metric-closure', 'children'),
    State('hidden-bionetwork', 'children')]
)
def make_subnetwork(n_clicks, reset_clicks, node_data, json_metric_closure, json_str_network):
    if n_clicks:
        print('make subnetwork')
        metric_closure = nx.node_link_graph(json.loads(json_metric_closure))
        network = nx.node_link_graph(json.loads(json_str_network))
        node_ids = [node['id'] for node in node_data]  # Get selected nodes in main network
        H = metric_closure.subgraph(node_ids)
        mst_edges = nx.minimum_spanning_edges(H, weight='distance', data=True)
        edges = chain.from_iterable(pairwise(d['path']) for u, v, d in mst_edges)
        T = network.edge_subgraph(edges)
        sub_network = T
        print(sub_network)
        cyto_sub_network, sub_cyto_nodes, sub_cyto_edges, subnetwork = make_cyto_elements(sub_network, 2, 300)
        return cyto_sub_network, {'display': 'none'}, {'display': 'block'}

    if reset_clicks:
        print('resetting')
        return [], {'display': 'block'}, {'display': 'none'}


@app.callback(
    Output('download-table-button', 'href'),
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
def download_csv():
    root_dir = os.getcwd()
    # downloads_dir = os.path.join(root_dir, 'downloads')
    # if not os.path.exists(downloads_dir):
    #    os.mkdir(downloads_dir)

    return flask.send_file(
        os.path.join(root_dir, 'downloads', 'table.csv'),
        mimetype='text/csv',
        attachment_filename='node_table.csv',
        as_attachment=True,
        cache_timeout=-1  # Prevent browser from caching previous file
    )


@app.callback(
    Output('download-subnetwork', 'href'),
    [Input('download-subnetwork', 'n_clicks')],
    [State('hidden-bionetwork', 'children'),
     State('sub-view', 'selectedNodeData')]
)
def update_sub_download_link(n_clicks, json_str_network, node_data):
    if n_clicks:
        rel_filename = os.path.join('downloads', 'sub_network.graphml')
        abs_filename = os.path.join(os.getcwd(), rel_filename)
        network = nx.node_link_graph(json.loads(json_str_network))
        sub_nodes = [node['id'] for node in node_data]
        sub_network = nx.subgraph(network, sub_nodes)
        nx.write_graphml(sub_network, abs_filename)

        return '/{}'.format(rel_filename)


@app.server.route('/downloads/<path:path>')
def download_sub_graphml(path):
    root_dir = os.getcwd()
    downloads_dir = os.path.join(root_dir, 'downloads')
    if not os.path.exists(downloads_dir):
        os.mkdir(downloads_dir)

    return flask.send_from_directory(
        os.path.join(downloads_dir),
        path,
        cache_timeout=-1  # Prevent browser from caching previous file
    )


@app.callback(
    Output('hidden-div', 'children'),
    [Input('sub-view', 'selectedNodeData')],
   # [State('cytoscape', 'selectedNodeData')]
)
def print_nodes(node_data):
    if node_data:
        print('miguebo  ', [node['id'] for node in node_data])
    return None


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
