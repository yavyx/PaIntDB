import json
import os
import sqlite3

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_cytoscape as cyto
import dash_html_components as html
import dash_table
import networkx as nx
import pandas as pd
from OmicsIntegrator import Graph
from dash.dash import no_update
from dash.dependencies import Output, Input, State, ALL
from dash_extensions import Download
from dash_extensions.snippets import send_file

import dash_app.vis_stylesheets as stylesheets
import go_enrichment.go_enrichment as goe
from dash_app.app import app  # Loads app variable from app script


def make_cyto_elements(network):
    """Takes a networkx network and outputs Cytoscape elements that can be visualized with Dash. Also creates selector
    classes according to the attributes and the layout coordinates."""
    # Get node degrees
    nx.set_node_attributes(network, dict(network.degree()), 'degree')

    # Convert nx network to cytoscape JSON
    json_elements = nx.readwrite.json_graph.cytoscape_data(network)['elements']

    # Make layout (much faster than default Cytoscape layouts)
    layout = nx.nx_agraph.graphviz_layout(network)
    nodes = json_elements['nodes']
    for node in nodes:
        node['data']['label'] = node['data']['shortName']  # Use short name as node label

        pos = layout[node['data']['id']]  # Extract positions from layout dict
        node['position'] = {'x': pos[0], 'y': pos[1]}  # Set positions

    elements = dict()
    elements['nodes'] = nodes
    elements['edges'] = json_elements['edges']
    return elements


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
                    {'label': 'RNA-Seq', 'value': 'RNASeq'},
                    {'label': 'Tn-Seq', 'value': 'TnSeq'}
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

    color_mapping = [html.Div(id='color-map'), html.Div(id='legend')]
    stylesheet = stylesheets.default
    # Add color mapping functionality for DE/Combined networks
    if network_params['type'] == 'gene_list':
        stylesheet = stylesheets.default
    elif network_params['type'] == 'rna_seq':
        color_mapping = [
            html.H5('Color Mapping'),
            html.Div(id='color-map'),
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
                    {'label': 'Experiment', 'value': 'experiment'},
                    {'label': 'Differential Expression', 'value': 'regulation'}
                ],
                value='experiment',
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
                    'width': '24vw',
                    'backgroundColor': '#a6edff',
                    'padding': '10px',
                    'display': 'inline-block',
                    'height': 'calc(100vh - 65px)',
                    'vertical-align': 'top',
                    'overflow': 'auto'
                },
                children=[
                    html.Div(color_mapping),
                    # legend,
                    dbc.Checklist(
                        id='show-labels',
                        options=[
                            {'label': 'Show node labels', 'value': 1}
                        ],
                        switch=True,
                        value=[]
                    ),
                    html.Hr(),
                    html.Div(
                        id='full-network-panel',
                        children=[
                            html.Div(
                                id='node-filters',
                                children=sidebar_filters
                            ),
                            html.P(id='num-selected-nodes', style={'padding-top': '5px'}),
                            dbc.Button('Make Sub-Network', id='make-subnetwork', color='primary',
                                       style={'display': 'none'})
                        ]
                    ),
                    html.Div(
                        id='subnetwork-btns',
                        style={'display': 'none'},
                        children=[
                            dbc.Checklist(id='include-extra-genes',
                                          options=[
                                              {'label': 'Include additional genes', 'value': 1}
                                          ],
                                          switch=True,
                                          value=[]
                                          ),
                            html.Abbr('Help',
                                      title=(('Include additional genes, called Steiner nodes, that are are not '
                                              'included in the original data, but help connect other genes that are. '
                                              'Useful to connect subnetworks with '
                                              'many smaller components.')),
                                      ),
                            html.Br(),
                            dbc.Checklist(id='include-low-confidence',
                                          options=[
                                              {'label': 'Include low-confidence interactions', 'value': 1}
                                          ],
                                          switch=True,
                                          value=[]
                                          ),
                            html.Abbr('Help',
                                      title=(('Include all interactions in the subnetwork, instead of prioritizing '
                                              'experimental interactions.')),
                                      ),
                            html.Br(),
                            dbc.Button(
                                'Return to selection',
                                id='reset-network',
                                color='primary',
                            )
                        ]
                    ),
                    html.Hr(),
                    dbc.DropdownMenu(
                        id='download-dropdown',
                        color='primary',
                        style={'padding-top': '5px'},
                        label='Download',
                        direction='right',
                        children=[
                            dbc.DropdownMenuItem('Network (.graphml)',
                                                 id='download-network'),
                            dbc.DropdownMenuItem('Network Image (.png)',
                                                 id='download-network-img'),
                            dbc.DropdownMenuItem('Table (.csv)',
                                                 id='download-table',
                                                 style={'display': 'none'})
                        ]
                    ),
                    Download(id='graphml-download'),
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
                                    children=cyto.Cytoscape(
                                        id='main-view',
                                        style={
                                            'width': '100%',
                                            'height': 'calc(100vh - 80px)'
                                        },
                                        stylesheet=stylesheet,
                                        maxZoom=5,
                                        minZoom=0.3,
                                        zoom=1,
                                        layout={'name': 'preset'},
                                        elements=cyto_network,
                                        boxSelectionEnabled=True
                                    )
                                ),
                            ]
                        ),
                        dbc.Row(
                            dbc.Col(
                                html.Div(id='node-details-table',
                                         style={'margin-top': '-25vh'}
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
    [Input('color-map', 'value'),
     Input('show-labels', 'value')]
)
def change_stylesheet(color_map, show_labels):
    if color_map == 'experiment':
        legend = app.get_asset_url('sig_source_legend.svg')
        stylesheet = stylesheets.combined
    else:
        legend = app.get_asset_url('de_legend.svg')
        stylesheet = stylesheets.fold_change
    if show_labels:
        return stylesheets.add_labels(stylesheet), legend
    else:
        return stylesheet, legend


@app.callback(
    [Output('full-network-panel', 'style'),
     Output('subnetwork-btns', 'style'),
     Output('main-view', 'elements'),
     Output('hidden-subnetwork', 'children'),
     Output('num-selected-nodes', 'children'),
     Output('make-subnetwork', 'style')],
    [Input({'type': 'filter', 'index': ALL}, 'value'),  # Pattern-matching all callbacks with filter type
     Input('make-subnetwork', 'n_clicks'),
     Input('include-low-confidence', 'value'),
     Input('include-extra-genes', 'value')],
    [State('node-details-df', 'children'),
     State('enrichment-results', 'children'),
     State('network-parameters', 'children'),
     State('cyto-network', 'children'),
     State('hidden-bionetwork', 'children')]
)
def select_nodes(values, subnetwork_clicks, low_confidence, extra_genes, node_details, enrichment_results,
                 network_params, cyto_network, bio_network):
    """Select nodes according to user selected filters. Creates subnetwork with selected nodes."""
    cyto_network = json.loads(cyto_network)
    enrichment_results = pd.read_json(enrichment_results)
    nodes = cyto_network['nodes']
    # edges = cyto_network['edges']
    network_df = pd.read_json(node_details)
    network_params = json.loads(network_params)
    strain = network_params['strain']
    network_type = network_params['type']
    query = []  # Query to filter nodes

    # Get selected filter input values
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
    query_str = ' & '.join(query)  # Join all queries
    if short_name:
        if not query:
            query_str += 'shortName in @short_name | index in @short_name'
        else:
            query_str += '| shortName in @short_name | index in @short_name'

    # Use query to select nodes
    if query_str:
        queried_nodes = network_df.query(query_str).index.tolist()
    else:
        queried_nodes = []

    # Select nodes
    for node in nodes:
        node['selected'] = True if node['data']['id'] in queried_nodes else False

    selected_msg = 'Selected {} out of {} nodes'.format(len(queried_nodes), len(nodes))

    # Display make network button after selecting nodes
    btn_display = {'display': 'block'} if len(queried_nodes) != 0 else {'display': 'none'}

    # Generate subnetwork when button is clicked.
    if subnetwork_clicks:
        cyto_sub_network, json_sub_network = make_subnetwork(queried_nodes, network_df, bio_network, strain,
                                                             network_type, low_confidence, extra_genes)
        # Throws warning if subnetwork sub_network is empty.
        if json_sub_network is None:
            selected_msg = dbc.Alert('Could not compute subnetwork using the selected nodes. Try selecting more nodes.',
                                     color='warning')
            cyto_sub_network = no_update
            json_sub_network = no_update
        # Return subnetwork
        else:
            selected_msg = ''
        return {'display': 'none'}, {'display': 'block'}, cyto_sub_network, json_sub_network, selected_msg, \
               btn_display
    # Return full network
    return {'display': 'block'}, {'display': 'none'}, cyto_network, no_update, selected_msg, btn_display


@app.callback(
    Output('make-subnetwork', 'n_clicks'),
    [Input('reset-network', 'n_clicks')])
def reset_subnetwork_clicks(n_clicks):
    """Reset subnetwork clicks to cycle through full network/subnetwork view."""
    return 0


def make_subnetwork(queried_nodes, network_df, json_str_network, strain, network_type, low_confidence, extra_genes):
    """Returns a subnetwork using the PCSF algorithm, using the user-selected nodes as terminals."""

    def make_prize_file(network_df, queried_nodes, network_type):
        """Generates .tsv file with node prizes for use with OmicsIntegrator."""
        if network_type == 'gene_list':
            # If there is no expression data, all prizes = 1
            network_df['prize'] = 1
            terminal_prizes = network_df.loc[network_df.index.isin(queried_nodes), 'prize']
        elif network_type == 'rna_seq' or network_type == 'combined':
            # Set prizes to expression values
            terminal_prizes = network_df.loc[network_df.index.isin(queried_nodes), ['log2FoldChange']]
            # The bigger the fold change, the bigger the prize
            terminal_prizes.log2FoldChange = abs(terminal_prizes.log2FoldChange)
            terminal_prizes = terminal_prizes.rename(columns={'log2FoldChange': 'prize'})
            if network_type == 'combined':
                # Set TnSeq prizes to the max prize
                terminal_prizes.loc[network_df['significanceSource'] == 'TnSeq', :] = terminal_prizes['prize'].max()
        terminal_prizes.to_csv('node_prizes.tsv', sep='\t')

    network = nx.node_link_graph(json.loads(json_str_network))
    # Make Graph object for prize-collecting Steiner forest (PCSF)
    graph = Graph(os.path.join('data', '{}_interactome.tsv'.format(strain)),  # Get interactome with costs
                  {'b': 10,  # b > 1 results in more terminal nodes in sub_network
                   'g': 0}  # g = 0 = disable degree cost correction
                  )
    make_prize_file(network_df, queried_nodes, network_type)
    graph.prepare_prizes('node_prizes.tsv')
    os.remove('node_prizes.tsv')  # Delete prize file (not needed anymore after running PCSF)
    vertex_indices, edge_indices = graph.pcsf()
    forest, augmented_forest = graph.output_forest_as_networkx(vertex_indices, edge_indices)
    # Include low confidence edges if selected by the user
    sub_network = augmented_forest if low_confidence else forest
    # If sub-network is empty, warning is shown
    if len(sub_network.nodes) == 0:
        return None, None

    # Sub-network includes extra genes (not in the input genes)
    if extra_genes:
        nodes = [node for node in sub_network.nodes]

        # Get extra gene information from database
        with sqlite3.connect('PaIntDB.db') as db_connection:
            descriptions = pd.read_sql_query("""SELECT id, product_name
                                                FROM protein
                                                WHERE id IN (%s)""" % ', '.join('?' * len(nodes)),
                                             con=db_connection, params=nodes)
            short_names = pd.read_sql_query("""SELECT id, name
                                               FROM interactor
                                               WHERE id IN (%s)""" % ', '.join('?' * len(nodes)),
                                            con=db_connection, params=nodes)
        # Format results to use as node attributes
        descriptions = descriptions.set_index('id').to_dict(orient='index')
        short_names = short_names.set_index('id').to_dict(orient='index')
        description_attr = dict()
        short_name_attr = dict()
        for key, value in descriptions.items():
            description_attr[key] = dict(description=value['product_name'])
        for key, value in short_names.items():
            short_name_attr[key] = dict(shortName=value['name'])
        nx.set_node_attributes(sub_network, description_attr)
        nx.set_node_attributes(sub_network, short_name_attr)

        # Set locus tags as short names for new genes
        for node in sub_network.nodes:
            if sub_network.nodes[node]['shortName'] is None:
                sub_network.nodes[node]['shortName'] = node

        sub_network.remove_edges_from(nx.selfloop_edges(sub_network))
    # Sub-network only includes genes in input genes
    else:
        sub_network = network.edge_subgraph(sub_network.edges())

    unfrozen_sub = nx.Graph(sub_network)  # Copy needed to remove orphan nodes
    unfrozen_sub.remove_nodes_from(list(nx.isolates(unfrozen_sub)))
    cyto_sub_network = make_cyto_elements(unfrozen_sub)
    json_sub_network = json.dumps(nx.node_link_data(unfrozen_sub))  # For downloading
    return cyto_sub_network, json_sub_network


@app.callback(
    [Output('node-details-table', 'children'),
     Output('filtered-node-details', 'children'),
     Output('download-table', 'style')],
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
        # Get selected nodes
        node_ids = [node['label'] for node in node_data]
        network_df = pd.read_json(node_details)
        if network_params['type'] == 'rna_seq' or network_params['type'] == 'combined':
            cols.extend(['log2FoldChange', 'padj'])
            network_df['log2FoldChange'] = network_df['log2FoldChange'].round(2)
            # network_df['padj'] = [sigfig.round(n, sigfigs=3) for n in network_df['padj']]

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
            dash_table.DataTable(
                data=filtered_df.to_dict('records'),
                columns=[{'name': i, 'id': i} for i in filtered_df.columns],
                fixed_rows={'headers': True},
                css=[{'selector': '.row', 'rule': 'margin: 0'}],  # Fixes left margin crop
                page_action='none',
                sort_action='native',
                style_as_list_view=True,
                style_table={
                    'maxHeight': '25vh',
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
        return nodes_table, filtered_df.to_json(), {'display': 'block'}
    else:
        return None, None, {'display': 'none'}


@app.callback(
    Output('csv-download', 'data'),
    [Input('download-table', 'n_clicks')],
    [State('filtered-node-details', 'children')]
)
def download_nodes_csv(n_clicks, json_df):
    if n_clicks:
        downloads_dir = os.path.join(os.getcwd(), 'downloads')
        if not os.path.exists(downloads_dir):
            os.mkdir(downloads_dir)
        nodes_df = pd.read_json(json_df)
        abs_filename = os.path.join(downloads_dir, 'node_details.csv')
        nodes_df.to_csv(abs_filename, index=False)
        return send_file(abs_filename)


@app.callback(
    Output('graphml-download', 'data'),
    Input('download-network', 'n_clicks'),
    State('hidden-subnetwork', 'children')
)
def download_sub_graphml(n_clicks, json_str_sub_network):
    if n_clicks:
        downloads_dir = os.path.join(os.getcwd(), 'downloads')
        if not os.path.exists(downloads_dir):
            os.mkdir(downloads_dir)
        rel_filename = os.path.join('downloads', 'network.graphml')
        abs_filename = os.path.join(os.getcwd(), rel_filename)
        sub_network = nx.node_link_graph(json.loads(json_str_sub_network))
        nx.write_graphml(sub_network, path=abs_filename)
        return send_file(abs_filename)


@app.callback(
    Output('main-view', 'generateImage'),
    [Input('download-network-img', 'n_clicks')]
)
def download_png(n_clicks):
    file_type = 'png'
    action = 'store'
    if n_clicks:
        file_type = 'png'
        action = 'download'
    return {'type': file_type, 'action': action, 'filename': 'subnetwork'}
