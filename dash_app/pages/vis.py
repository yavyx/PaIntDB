import json
import os
from math import sqrt

from dash.dependencies import Output, Input, State, ALL
from dash.dash import no_update
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_cytoscape as cyto
from dash_extensions import Download
from dash_extensions.snippets import send_file, send_data_frame
import dash_html_components as html
import dash_table
import networkx as nx
from OmicsIntegrator import Graph
import pandas as pd

import dash_app.vis_stylesheets as stylesheets
from dash_app.app import app  # Loads app variable from app script
import go_enrichment.go_enrichment as goe

# For development  # TODO: delete
pd.set_option('display.max_columns', None)


def make_cyto_elements(network):
    """Takes a networkx network and outputs Cytoscape elements that can be visualized with Dash. Also creates selector
    classes according to the attributes and the layout coordinates."""
    # Get node degrees
    nx.set_node_attributes(network, dict(network.degree()), 'degree')

    # Convert nx network to cytoscape JSON
    json_elements = nx.readwrite.json_graph.cytoscape_data(network)['elements']

    # Make layout (much faster than default Cytoscape layouts)
    # layout = nx.nx_agraph.graphviz_layout(network)
    # layout = nx.drawing.kamada_kawai_layout(network, scale=1000)
    layout = nx.drawing.spring_layout(network, scale=1000, seed=1, k=3/sqrt(len(network)))
    nodes = json_elements['nodes']
    for node in nodes:
        node['data']['label'] = node['data']['shortName']  # Use short name as node label

        pos = layout[node['data']['id']]  # Extract positions from layout dict
        node['position'] = {'x': pos[0], 'y': pos[1]}  # Set positions

    edges = json_elements['edges']
    elements = nodes + edges

    return elements, nodes, edges


def make_vis_layout(network_df, enrichment_results, cyto_network, network_params):
    """Generates a custom layout depending on the network type."""
    network_params = json.loads(network_params)

    # This filter is only used with RNASeq/Combined networks
    regulation_filter = html.Details(
        id='diff-exp-details',
        open=True,
        children=[
            html.Summary('By differential expression'),
            dbc.Checklist(
                # ID type for pattern matching callback
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
                # ID type for pattern matching callback
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
                    # ID type for pattern-matching callback
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
                    # ID type for pattern-matching callback
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
                    # ID type for pattern-matching callback
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

    # Add extra filters for DE/Combined networks
    if network_params['type'] == 'rna_seq' or network_params['type'] == 'combined':
        sidebar_filters.extend([regulation_filter, html.Br()])
        if network_params['type'] == 'combined':
            sidebar_filters.extend([source_filter])

    # Add color mapping functionality for DE/Combined networks
    color_mapping = None
    stylesheet = stylesheets.default
    # if network_params['type'] == 'gene_list':
    #    stylesheet = stylesheets.default
    if network_params['type'] == 'rna_seq':
        color_mapping = [
            html.H5('Color Mapping'),
            html.Div(style={'padding': '10px'},
                     children=html.Img(
                         src=app.get_asset_url('de_legend.svg'),
                         id='legend',
                         width=100)
                     )
        ]
        stylesheet = stylesheets.fold_change
    elif network_params['type'] == 'combined':
        color_mapping = [
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

    # Layout begins here
    return html.Div(
        [
            html.Div(
                style={
                    'width': '25vw',
                    'backgroundColor': '#a6edff',
                    'padding': '10px',
                    'display': 'inline-block',
                    'height': '95vh',
                    'vertical-align': 'top',
                    'overflow': 'auto'
                },
                children=[
                    html.Div(id='color-map-div',
                             children=color_mapping
                             ),
                    html.Hr(),
                    html.Div(
                        id='full-network-panel',
                        children=[
                            html.Div(
                                id='node-filters',
                                children=sidebar_filters
                            ),
                            html.Br(),
                            html.P(id='num-selected-nodes'),
                            dbc.Button('Make Sub-Network', id='make-subnetwork', color='primary')
                        ]
                    ),
                    html.Div(
                        id='subnetwork-btns',
                        style={'display': 'none'},
                        children=[
                            dbc.Button('Return to selection', id='reset-network', color='primary'),
                            dbc.DropdownMenu(
                                id='download-dropdown',
                                color='primary',
                                style={'padding-top': '5px'},
                                label='Download',
                                direction='right',
                                children=[
                                    dbc.DropdownMenuItem('Download Table (.csv)',
                                                         id='download-table'),
                                    dbc.DropdownMenuItem('Network (.graphml)',
                                                         id='download-network'),
                                    dbc.DropdownMenuItem('Network (.svg)',
                                                         id='download-network-img')
                                ]
                            ),
                        ]
                    ),
                    Download(id='graphml-download2'),
                    Download(id='csv-download'),
                    # Hidden Divs to store node details and subnetwork for download
                    html.Div(id='filtered-node-details', style={'display': 'none'}),
                    html.Div(id='hidden-subnetwork', style={'display': 'none'})
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
                                    children=cyto.Cytoscape(
                                        id='main-view',
                                        style={
                                            'width': '100%',
                                            'height': '90vh',
                                        },
                                        stylesheet=stylesheet,
                                        maxZoom=5,
                                        minZoom=0.3,
                                        zoom=1,
                                        layout={'name': 'preset'},
                                        elements=cyto_network['elements'],
                                        boxSelectionEnabled=True
                                    )
                                ),
                            ]
                        ),
                        dbc.Row(
                            dbc.Col(
                                html.Div(id='node-details-table',
                                         style={'margin-top': '-30vh'}
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
     Output('legend', 'src')],
    [Input('color-map', 'value')]
)
def change_color_map(value):
    if value == 'ss':
        source = app.get_asset_url('sig_source_legend.svg')
        return stylesheets.combined, source
    else:
        source = app.get_asset_url('de_legend.svg')
        return stylesheets.fold_change, source


@app.callback(
    [Output('full-network-panel', 'style'),
     Output('subnetwork-btns', 'style'),
     Output('main-view', 'elements'),
     Output('hidden-subnetwork', 'children'),
     Output('num-selected-nodes', 'children')],
    [Input({'type': 'filter', 'index': ALL}, 'value'),  # Pattern-matching all callbacks with filter type
     Input('make-subnetwork', 'n_clicks')],
    [State('main-view', 'selectedNodeData'),
     State('node-details-df', 'children'),
     State('enrichment-results', 'children'),
     State('network-parameters', 'children'),
     State('cyto-network', 'children'),
     State('hidden-bionetwork', 'children')]
)
def select_nodes(values, subnetwork_clicks, node_data, node_details, enrichment_results,
                 network_params, cyto_network, bio_network):
    """Select nodes according to user selected filters. Creates subnetwork with selected nodes."""
    cyto_network = json.loads(cyto_network)
    enrichment_results = pd.read_json(enrichment_results)
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
            significance_source = []
    else:
        regulation, significance_source = [], []

    # Add queries depending on GUI filter selections
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

    # Use query to select nodes
    if query_str:
        queried_nodes = network_df.query(query_str).index.tolist()
    else:
        queried_nodes = nodes

    n_selected = 0  # Selected nodes counter
    for node in nodes:
        if node['data']['id'] in queried_nodes:
            node['selected'] = True
            n_selected += 1
        else:
            node['selected'] = False

    selected_msg = 'Selected {} out of {} nodes'.format(n_selected, len(nodes))

    # Generate subnetwork when button is clicked.
    if subnetwork_clicks:
        cyto_sub_network, json_sub_network = make_subnetwork(node_data, network_df, bio_network, strain)
        # Throws warning if subnetwork solution is empty.
        if json_sub_network is None:
            selected_msg = dbc.Alert('Could not compute subnetwork using the selected nodes. Try selecting more nodes.',
                                     color='danger')
            return {'display': 'none'}, {'display': 'none'}, no_update, no_update, selected_msg
        # Return subnetwork
        return {'display': 'none'}, {'display': 'block'}, cyto_sub_network, json_sub_network, ''
    # Return full network
    return {'display': 'block'}, {'display': 'none'}, nodes + edges, no_update, selected_msg


@app.callback(
    Output('make-subnetwork', 'n_clicks'),
    [Input('reset-network', 'n_clicks')])
def reset_subnetwork_clicks(n_clicks):
    """Reset subnetwork clicks to cycle through full network/subnetwork view."""
    return 0


def make_subnetwork(node_data, network_df, json_str_network, strain):
    """Returns a subnetwork using the PCSF algorithm, using the user-selected nodes as terminals."""

    def make_prize_file(network_df, node_data):
        """Generates .tsv file with node prizes for use with OmicsIntegrator."""
        # User-selected nodes are terminals
        terminals = [node['id'] for node in node_data]
        print(terminals)
        terminal_prizes = network_df.loc[network_df.index.isin(terminals), ['log2FoldChange']]
        # The bigger the fold change, the bigger the prize
        terminal_prizes.log2FoldChange = abs(terminal_prizes.log2FoldChange)
        # Set TnSeq prizes to the max log2FoldChange
        terminal_prizes.loc[terminal_prizes['log2FoldChange'].isna(), :] = max(network_df['log2FoldChange'])
        terminal_prizes = terminal_prizes.rename(columns={'log2FoldChange': 'prize'})
        terminal_prizes.to_csv(os.path.join('temp_data', 'node_prizes.tsv'), sep='\t')

    network = nx.node_link_graph(json.loads(json_str_network))
    print(len(network.nodes()))
    # Make Graph object for prize-collecting Steiner forest (PCSF)
    graph = Graph(os.path.join('data', '{}_interactome.tsv'.format(strain)),
                  {'b': 5})  # b > 1 results in more terminal nodes in solution
    make_prize_file(network_df, node_data)
    graph.prepare_prizes(os.path.join('temp_data', 'node_prizes.tsv'))
    os.remove(os.path.join('temp_data', 'node_prizes.tsv'))
    vertex_indices, edge_indices = graph.pcsf()
    forest, augmented_forest = graph.output_forest_as_networkx(vertex_indices, edge_indices)
    # If solution is empty, warning is shown
    if len(forest.nodes) == 0:
        return None, None
    sub_network = network.edge_subgraph(forest.edges())
    cyto_sub_network, sub_cyto_nodes, sub_cyto_edges = make_cyto_elements(sub_network)
    json_sub_network = json.dumps(nx.node_link_data(sub_network))  # For downloading
    return cyto_sub_network, json_sub_network


@app.callback(
    [Output('node-details-table', 'children'),
     Output('filtered-node-details', 'children')],
    [Input('main-view', 'selectedNodeData')],
    [State('node-details-df', 'children'),
     State('network-parameters', 'children')]
)
def show_node_details(node_data, node_details, network_params):
    """Filters the network DataFrame with the user-selected nodes and returns a DataTable."""
    if node_data:
        # Columns to display
        cols = ['shortName', 'description']
        network_params = json.loads(network_params)
        if network_params['type'] == 'rna_seq' or network_params['type'] == 'combined':
            cols.extend(['log2FoldChange', 'padj'])
        # Get selected nodes
        node_ids = [node['label'] for node in node_data]
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

        nodes_table = [
            html.H5('Selected Node(s) Details'),
            dash_table.DataTable(
                data=filtered_df.to_dict('records'),
                columns=[{'name': i, 'id': i} for i in filtered_df.columns],
                fixed_rows={'headers': True},
                css=[{'selector': '.row', 'rule': 'margin: 0'}],  # Fixes left margin crop
                page_action='none',
                sort_action='native',
                style_as_list_view=True,
                style_table={
                    'maxHeight': '30vh',
                    'overflowY': 'auto'
                },
                style_cell={'textAlign': 'left',
                            'minWidth': '150px',
                            'width': '150px',
                            'maxWidth': '150px',
                            'font-family': 'sans-serif'
                            },
                style_header={'backgroundColor': 'rgb(166, 237, 255)',
                              'fontWeight': 'bold'
                              },
                style_data={'whiteSpace': 'normal',
                            'table-layout': 'fixed'
                            }
            )
        ]
        return nodes_table, filtered_df.to_json()
    else:
        return None, None


@app.callback(
    Output('csv-download', 'data'),
    [Input('download-table', 'n_clicks')],
    [State('filtered-node-details', 'children')]
)
def download_csv(n_clicks, json_df):
    if n_clicks:
        # downloads_dir = os.path.join(os.getcwd(), 'downloads')
        # if not os.path.exists(downloads_dir):
        #     os.mkdir(downloads_dir)
        node_details = pd.read_json(json_df)
        # abs_filename = os.path.join(downloads_dir, 'node_details.csv')
        print(node_details.head())
        # TODO: Fix index issue
        return send_data_frame(node_details.to_csv, 'node_details.csv')


@app.callback(
    Output('graphml-download2', 'data'),
    Input('download-network', 'n_clicks'),
    State('hidden-subnetwork', 'children')
)
def download_sub_graphml(n_clicks, json_str_sub_network):
    if n_clicks:
        downloads_dir = os.path.join(os.getcwd(), 'downloads')
        if not os.path.exists(downloads_dir):
            os.mkdir(downloads_dir)
        rel_filename = os.path.join('downloads', 'subnetwork.graphml')
        abs_filename = os.path.join(os.getcwd(), rel_filename)
        sub_network = nx.node_link_graph(json.loads(json_str_sub_network))
        nx.write_graphml(sub_network, path=abs_filename)
        return send_file(abs_filename)


# @app.callback(
#     Output('hidden-div', 'children'),
#     [Input('main-view', 'selectedNodeData')],
# )
# def print_nodes(node_data):
#     if node_data:
#         print('miguebo  ', [node['id'] for node in node_data])
#     return None

