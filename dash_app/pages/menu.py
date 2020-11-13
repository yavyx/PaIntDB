import base64
from datetime import datetime
import io
import json
import os

from dash.dependencies import Output, Input, State
from dash.exceptions import PreventUpdate
from dash_extensions import Download
from dash_extensions.snippets import send_file
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import networkx as nx
import pandas as pd


from bio_networks.network_generator import BioNetwork, DENetwork, CombinedNetwork
from go_enrichment.go_enrichment import run_go_enrichment
from dash_app.app import app  # Loads app variable from app script


layout = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Br(),
                        'Select your input data: ',
                        dbc.RadioItems(
                            id='network-type',
                            options=[
                                {'label': 'Gene List', 'value': 'basic'},
                                {'label': 'Differentially-Expressed (DE) Gene List', 'value': 'DE'},
                                {'label': 'Combined (DE genes and TnSeq genes)', 'value': 'combined'}
                            ],
                            value='basic'),
                        html.Br(),
                        html.P('Your genes must be in the first column of a CSV file.'),
                        dcc.Upload(
                            id='gene-list-upload',
                            children=dbc.Button(
                                id='upload-button',
                                children='1. Upload Gene List',
                                color='primary')),
                        html.Div([
                            dcc.Upload(
                                id='tnseq-gene-list-upload',
                                children=dbc.Button(
                                    id='tnseq_upload-button',
                                    children='Upload TnSeq Gene List',
                                    color='primary'))
                        ], style={'marginTop': '1vh'}
                        )

                    ],
                    width=4
                ),
                dbc.Col(
                    [
                        html.Br(),
                        html.Div(id='data-upload-output', style={'height': '30vh'})
                    ],
                    width=5
                ),
            ]
        ),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Div(
                            ['Select the strain:',
                             dbc.RadioItems(
                                 id='strain',
                                 options=[
                                     {'label': 'PAO1', 'value': 'PAO1'},
                                     {'label': 'PA14', 'value': 'PA14'},
                                 ],
                                 value='PAO1'
                             )]
                        ),
                    ],
                    width=2
                ),
                dbc.Col(
                    [
                        html.Div(
                            ['Select the network order: ',
                             html.Abbr('?',
                                       title=(('Zero-order: Maps direct interactions between your queried genes. '
                                               '\nRecommended for long lists (>200 genes).'
                                               '\n\nFirst-order: Uses your queried genes as "seed" genes and finds any '
                                               'interaction between them and the other genes in the database. '
                                               '\nRecommended for short lists (<200 genes).'))
                                       ),
                             dbc.RadioItems(
                                 id='order',
                                 options=[
                                     {'label': 'Zero-order', 'value': 0},
                                     {'label': 'First-order', 'value': 1},
                                 ],
                                 value=0
                             )]
                        )
                    ],
                    width=2
                ),
                dbc.Col(
                    [
                        html.Div(
                            [
                                'Select the interaction detection method: ',
                                html.Abbr('?',
                                          title=(('Choose which interactions you want to use to generate the network.'
                                                  '\n\n"All" includes interactions predicted computationally.'
                                                  '\n\nExperimentally-verified interactions have the highest '
                                                  'confidence, but result in smaller networks.'))
                                          ),
                                dbc.RadioItems(
                                    id='detection-method',
                                    options=[
                                        {'label': 'All', 'value': 3},
                                        {'label': 'Experimental', 'value': 2},
                                    ],
                                    value=3
                                )
                            ])
                    ],
                    width=3
                )
            ]
        ),
        dbc.Row(
            dbc.Col(
                html.Div(
                    [
                        # html.Br(),
                        dbc.Checklist(
                            id='metabolites',
                            options=[
                                {'label': 'Include metabolites', 'value': 1}
                            ],
                            style={'display': 'none'}  # Removed for now
                        )
                    ]
                )
            )
        ),
        html.Br(),
        dbc.Button('2. Make Network', id='make-network-btn', color='primary'),
        html.Br(),
        dcc.Loading(id='loading',
                    children=html.Div(id='make-network-message'),
                    type='dot'),
        html.Div(
            [
                dbc.Button('Download (.graphml)', id='download-btn', color='link', style={'margin-top': '1vh'}),
                Download(id='graphml-download1')
            ],
            id='download',
            style={'display': 'none'},
        ),
        html.Div(
            [
                html.Hr(),
                'Select the genes for enrichment: ',
                dbc.RadioItems(id='enrichment-options'),
                html.Br(),
                dbc.Button('3. Run GO Term Enrichment', id='run-enrichment', color='primary')
            ],
            id='enrichment-btns',
            style={'display': 'none'}
        ),
        html.Br(),
        dcc.Loading(id='loading',
                    children=html.Div(id='enrichment-loading'),
                    type='dot'),
        html.Br(),
        dbc.Button('4. Explore Network', id='explore-btn', href='/vis', color='primary', block=False,
                   style={'display': 'none'}),
    ],
    fluid=True,
    style={'background-color': '#ededed',
           'height': '95vh',
           'overflow': 'auto'}
)


def parse_gene_list(contents, network_type):
    """Parses the uploaded gene list, returns a Bootstrap table and a Pandas DataFrame."""
    # TODO: Add checks for table headers
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    genes_df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
    cols = [col for col in genes_df.columns if col not in ['pvalue', 'padj']]  # select all columns except p-values
    genes_df.loc[:, cols] = genes_df[cols].round(2)
    if network_type == 'DE' or network_type == 'combined':
        # Check RNASeq headers are there
        if not {'log2FoldChange', 'padj'}.issubset(genes_df.columns):
            return dbc.Alert('Check your header names. They should include "log2FoldChange" and "padj"',
                             color='danger', style={'display': 'inline-block'}), []

    small_df = genes_df.head()  # smaller df to display on app
    table = dash_table.DataTable(
        data=small_df.to_dict('records'),
        columns=[{"name": i, "id": i} for i in small_df.columns],
        style_table={
            'maxHeight': '20vh',
        },
        style_cell={
            'font-family': 'sans-serif',
            'textAlign': 'left'
        }
    )

    upload_contents = html.Div(
        [
            dbc.Alert('Your list was uploaded successfully!', color='primary', dismissable=True,
                      style={'display': 'inline-block'}),
            table
        ]
    )
    return upload_contents, genes_df


def parse_tnseq_list(contents, filename):
    """Parses the uploaded TnSeq gene list"""
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    tnseq_df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
    tnseq_genes = tnseq_df.iloc[:, 0].tolist()
    upload_msg = html.Div(['Your TnSeq genes were uploaded successfully!',
                           filename])
    return upload_msg, tnseq_genes


@app.callback(
    Output('tnseq-gene-list-upload', 'style'),
    [Input('network-type', 'value')]
)
def show_tnseq_upload_btn(network_type):
    """Show TnSeq upload button when combined networks are selected."""
    return {'display': 'block'} if network_type == 'combined' else {'display': 'none'}


@app.callback(
    Output('make-network-btn', 'style'),
    [Input('gene-list-upload', 'contents'),
     Input('tnseq-gene-list-upload', 'contents'),
     Input('network-type', 'value')]
)
def show_make_network_btn(contents, tnseq_contents, network_type):
    """Shows make network button upload after upload is read."""
    if network_type == 'combined':
        if contents and tnseq_contents:
            return {'display': 'block'}
        else:
            return {'display': 'none'}
    elif network_type == 'basic' or network_type == 'DE':
        if contents:
            return {'display': 'block'}
        else:
            return {'display': 'none'}


@app.callback(
    Output('data-upload-output', 'children'),
    [Input('gene-list-upload', 'contents')],
    [State('network-type', 'value')]
)
def upload_message(contents, network_type):
    """Returns a successful message after file was uploaded."""
    if contents:
        try:
            small_table, genes_df = parse_gene_list(contents, network_type)
            return small_table
        except ValueError:
            return dbc.Alert('There was a problem uploading your file. Check that it is the correct format.',
                             color='danger', style={'display': 'inline-block'})


@app.callback(
    [Output('enrichment-btns', 'style'),
     Output('download', 'style'),
     Output('enrichment-options', 'options'),
     Output('enrichment-options', 'value'),
     Output('make-network-message', 'children'),
     Output('hidden-bionetwork', 'children'),
     Output('node-details-df', 'children'),
     Output('network-parameters', 'children'),
     Output('genes-of-interest', 'children')],
    [Input('make-network-btn', 'n_clicks')],
    [State('network-type', 'value'),
     State('strain', 'value'),
     State('order', 'value'),
     State('detection-method', 'value'),
     State('metabolites', 'value'),
     State('gene-list-upload', 'contents'),
     State('tnseq-gene-list-upload', 'contents'),
     State('gene-list-upload', 'filename'),
     State('tnseq-gene-list-upload', 'filename')]
)
def build_network(n_clicks, network_type, strain, order, detection_method, metabolites, rnaseq_contents,
                  tnseq_contents, rnaseq_filename, tnseq_filename):
    """Generates a network every time the make network button is clicked. Serializes results to JSON
    and stores them in hidden divs. Shows download and explore network button."""
    if n_clicks is None:
        raise PreventUpdate

    start_time = datetime.now()
    upload_msg, genes_df = parse_gene_list(rnaseq_contents, rnaseq_filename)
    genes_df.rename(columns={genes_df.columns[0]: 'gene'}, inplace=True)
    gene_list = genes_df.gene.tolist()
    if network_type == 'basic':
        bio_network = BioNetwork(gene_list=gene_list, strain=strain, order=order, detection_method=detection_method,
                                 metabolites=metabolites)
    elif network_type == 'DE':
        bio_network = DENetwork(gene_list=gene_list, strain=strain, order=order, detection_method=detection_method,
                                metabolites=metabolites, de_genes_df=genes_df)
    elif network_type == 'combined':
        upload_msg, tnseq_genes = parse_tnseq_list(tnseq_contents, tnseq_filename)
        bio_network = CombinedNetwork(gene_list=gene_list, strain=strain, order=order,
                                      detection_method=detection_method, metabolites=metabolites,
                                      de_genes_df=genes_df, tnseq_gene_list=tnseq_genes)
    else:
        bio_network = None

    if metabolites:
        mapping_msg = html.Div('''{} genes were mapped to the network out of {} genes in your list.\n{} 
                                       metabolites were mapped to these genes.'''
                               .format(len(bio_network.mapped_genes),
                                       len(bio_network.genes_of_interest),
                                       len(bio_network.mapped_metabolites)
                                       )
                               )
    else:
        mapping_msg = html.Div('{} genes were mapped to the network out of {} genes in your list.'
                               .format(len(bio_network.mapped_genes),
                                       len(bio_network.genes_of_interest))
                               )
    end_time = datetime.now()

    print("order = {}, detection_method = {}, metabolites = {}".format(order, detection_method, str(metabolites)))
    print(end_time - start_time)

    enrichment_options = [
            {'label': 'Full gene list ({} genes)'.format(len(bio_network.genes_of_interest)), 'value': 'all'},
            {'label': 'Genes mapped to network ({} genes)'.format(len(bio_network.mapped_genes)), 'value': 'network'},
        ]

    json_network = json.dumps(nx.node_link_data(bio_network.network))
    network_df = bio_network.network_df
    genes_of_interest = bio_network.genes_of_interest
    network_params = {'strain': bio_network.strain, 'type': bio_network.network_type}

    return {'display': 'block'}, {'display': 'block'}, enrichment_options, 'all', mapping_msg, json_network, \
           network_df.to_json(), json.dumps(network_params), json.dumps(genes_of_interest)


@app.callback(
    Output('graphml-download1', 'data'),
    [Input('download-btn', 'n_clicks')],
    [State('gene-list-upload', 'filename'),
     State('hidden-bionetwork', 'children')]
)
def download_graphml(n_clicks, filename, json_str_network):
    """Generates and sends a graphml file for download."""
    if n_clicks:
        # Create downloads directory if there isn't one
        downloads_dir = os.path.join(os.getcwd(), 'downloads')
        if not os.path.exists(downloads_dir):
            os.mkdir(downloads_dir)

        rel_filename = os.path.join('downloads', '{}_network.graphml'.format(filename[:-4]))
        abs_filename = os.path.join(os.getcwd(), rel_filename)
        network = nx.node_link_graph(json.loads(json_str_network))
        nx.write_graphml(network, path=abs_filename)
        return send_file(abs_filename)


@app.callback(
    [Output('enrichment-results', 'children'),
     Output('enrichment-loading', 'children')],
    [Input('run-enrichment', 'n_clicks')],
    [State('network-parameters', 'children'),
     State('genes-of-interest', 'children'),
     State('enrichment-options', 'value'),
     State('hidden-bionetwork', 'children')]
)
def run_enrichment(n_clicks, network_params, genes_of_interest, gene_list, json_network):
    if n_clicks is None:
        raise PreventUpdate
    # Load JSON data
    network_params = json.loads(network_params)
    strain = network_params['strain']
    network = nx.node_link_graph(json.loads(json_network))
    # Use full gene list or genes mapped to network depending on user selection
    enrichment_genes = json.loads(genes_of_interest) if gene_list == 'all' else network.nodes
    # Run enrichment
    enrichment_results, goea_results = run_go_enrichment(strain, enrichment_genes)
    # Keep only overrepresented terms (remove underrepresented)
    enrichment_results = enrichment_results.loc[enrichment_results['enrichment'] == 'e', :]
    return enrichment_results.to_json(), 'Found {} enriched GO terms.'.format(len(enrichment_results))


@app.callback(
    Output('explore-btn', 'style'),
    [Input('enrichment-results', 'children')]
)
def show_explore_network_btn(enrichment_results):
    """Shows explore network button after enrichment is done."""
    return {'display': 'inline-block'} if enrichment_results else {'display': 'none'}


