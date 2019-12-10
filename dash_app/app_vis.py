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
    df['regulation'] = [None if sig == 'TnSeq'else reg for sig, reg in zip(df['significanceSource'], df['regulation'])]
    print(df.head())
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
    layout = nx.spring_layout(network, k=1 / sqrt(len(network)), scale=1000, center=[450, 450])
    nodes = json_elements['nodes']
    for node in nodes:
        node['data']['label'] = node['data']['shortName']  # Use short name as node label

        pos = layout[node['data']['id']]  # Extract positions from layout dict
        node['position'] = {'x': pos[0], 'y': pos[1]}

    edges = json_elements['edges']
    elements = nodes + edges

    return elements, nodes, edges, network


# temp_network = nx.read_graphml('temp_data/Biofilm_DJK6_vs_Biofilm_network.graphml')
# temp_network = nx.read_graphml('temp_data/Biofilm_DJK5_vs_Biofilm_network (1).graphml')
# network_path = 'corries_AZM_combined.graphml'
network_path = '/home/javier/PycharmProjects/PaIntDB/temp_data/mediarpmi.treatmentazm_exp_network.graphml'
temp_network = nx.read_graphml(network_path)
cyto_elements, cyto_nodes, cyto_edges, network = make_cyto_elements(temp_network)
network_df = make_network_df(network)
enrichment_results, goea_results = goe.run_go_enrichment('PAO1', list(network_df.index))


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
                                        html.Button('Make Selection',
                                                    id='make-selection',
                                                    style={'display': 'inline-block'}),
                                        html.P(id='num-selected-nodes'),
                                        html.Br(),
                                        #html.Button('Run GO Term Enrichment',
                                        #            id='run-enrichment',
                                        #            style={'display': 'inline-block'}),
                                        #html.Div(id='enrichment-results', style={'display': 'none'})
                                    ])
                                    #style={'text-align': 'center'})
                            ],
                            style={'margin': '2vh'}
                        ),
                        width=3,
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
                        width=9
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
    [Output('cytoscape', 'elements'),
     Output('num-selected-nodes', 'children')],
    [Input('make-selection', 'n_clicks')],
    [State('source-selection', 'value'),
     State('location-selection', 'value'),
     State('diff-exp-selection', 'value'),
     State('enrichment-selection', 'value')]
)
def select_nodes(n_clicks, significance_source, location, regulation, enriched_terms):
    """Select nodes according to significance source."""
    nodes = cyto_nodes  # Make new variable to avoid modifying the global cyto_nodes

    query = []  # query to filter nodes
    if significance_source:
        significance_source.append('both')
        query.append('significanceSource in @significance_source')

    if location:
        query.append('localization in @location')

    if regulation:
        query.append('regulation in @regulation')

    if enriched_terms:
        # Get genes associated with selected GO term(s)
        print(enriched_terms)
        genes_in_term = enrichment_results.loc[enrichment_results['name'].isin(enriched_terms), 'study_items'].iloc[0]
        print(len(genes_in_term))
        query.append('index in @genes_in_term')

    query_str = ' & '.join(query)
    print(query_str)
    if query_str:
        queried_nodes = list(network_df.query(query_str).index)
        selected_msg = 'Selected {} out of {} nodes'.format(len(queried_nodes), len(cyto_nodes))
    else:
        queried_nodes = nodes
        selected_msg = ''

    for node in nodes:
        node['selected'] = True if node['data']['id'] in queried_nodes else False

    return [nodes + cyto_edges, selected_msg]


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
                       .rename(columns={'index': 'Locus Tag',
                                        'shortName': 'Short Name',
                                        'description': 'Descripton',
                                        'log2FoldChange': 'Log2 Fold Change',
                                        'padj': 'Adjusted p-value',
                                        'ncbi': 'NCBI Accession #',
                                        'uniprotkb': 'UniProtKB Accession #'})
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
