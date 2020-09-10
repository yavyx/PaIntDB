import base64
from datetime import datetime
import io
import json
import os

from dash.dependencies import Output, Input, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import flask
import networkx as nx
from networkx.readwrite import json_graph
import pandas as pd

import bio_networks.network_generator as ng
from dash_app.app import app  # Loads app variable from app script
import dash_app.apps.app_vis as app_vis


layout = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(
                    [
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
                                children='Upload Gene List')),
                        html.Div([
                            dcc.Upload(
                                id='tnseq-gene-list-upload',
                                children=dbc.Button(
                                    id='tnseq_upload-button',
                                    children='Upload TnSeq Gene List'))
                        ], style={'marginTop': '1vh'}
                        )

                    ],
                    width=4
                ),
                dbc.Col(
                    html.Div(id='data-upload-output', style={'height': '30vh'}),
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
                             html.Abbr("?",
                                       title=('Zero-order: Maps direct interactions between your queried genes.\n'
                                              'Recommended for long lists (>200 genes).\n\n'
                                              'First-order: Uses your queried genes as "seed" genes and finds any\n'
                                              'interaction between them and the other genes in the database.\n'
                                              'Recommended for short lists (<200 genes).')
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
                                html.Abbr("?",
                                          title=('Choose which interactions you want to use to generate the '
                                                 'network.\n\nAll includes interactions predicted computationally.\n\n'
                                                 'Experimentally-verified interactions have the highest confidence,\n '
                                                 'but result in smaller networks.\n\n ')
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
                        html.Br(),
                        dbc.Checklist(
                            id='metabolites',
                            options=[
                                {'label': 'Include metabolites', 'value': 1}
                            ]
                        )
                    ]
                )
            )
        ),
        html.Br(),
        dbc.Button('Make Network', id='make-network'),
        html.Br(),
        dcc.Loading(id='loading',
                    children=html.Div(id='make-network-message'),
                    type='dot'),
        html.Hr(),
        # dbc.Button('Explore Network', id='explore-button', href='/vis'),
        html.Br(),
        dbc.Button(
            html.A('Download Network(GraphML)',
                   id='download-link'
                   )
        ),
    ],
    fluid=True
)


def parse_gene_list(contents, filename):
    """Parses the uploaded gene list, returns a Bootstrap table and a Pandas DataFrame."""
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    genes_df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
    # Formatting DF
    # genes_df['pvalue'] = genes_df['pvalue'].map('{:.3g}'.format)  # Change to have 2 decimals
    # genes_df['padj'] = genes_df['padj'].map('{:.3g}'.format)
    # genes_df.loc[:, ['pvalue', 'padj']] = genes_df[['pvalue', 'padj']].astype(float)
    small_df = genes_df.head()  # smaller df to display on app
    cols = [col for col in small_df.columns if col not in ['pvalue', 'padj']]  # select all columns except pvalues
    small_df.loc[:, cols] = small_df[cols].round(2)
    table = dash_table.DataTable(
        data=small_df.to_dict('records'),
        columns=[{"name": i, "id": i} for i in small_df.columns],
        style_table={
            'maxHeight': '20vh',
        }
    )

    upload_contents = html.Div(
        [
            html.P('Your list was uploaded successfully!'),
            html.P(filename),
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


def make_network(network_type,
                 strain,
                 order,
                 detection_method,
                 metabolites,
                 rnaseq_filename,
                 rnaseq_contents,
                 tnseq_filename=None,
                 tnseq_contents=None):
    """Generates a BioNetwork. Serializes the result to JSON to use in the vis module."""
    start_time = datetime.now()
    upload_msg, genes_df = parse_gene_list(rnaseq_contents, rnaseq_filename)
    genes_df.rename(columns={genes_df.columns[0]: 'gene'}, inplace=True)
    gene_list = genes_df.gene.tolist()
    if network_type == 'basic':
        bio_network = ng.BioNetwork(gene_list=gene_list,
                                    strain=strain,
                                    order=order,
                                    detection_method=detection_method,
                                    metabolites=metabolites)
    elif network_type == 'DE':
        bio_network = ng.DENetwork(gene_list=gene_list,
                                   strain=strain,
                                   order=order,
                                   detection_method=detection_method,
                                   metabolites=metabolites,
                                   de_genes_df=genes_df)
    elif network_type == 'combined':
        upload_msg, tnseq_genes = parse_tnseq_list(tnseq_contents, tnseq_filename)
        bio_network = ng.CombinedNetwork(gene_list=gene_list,
                                         strain=strain,
                                         order=order,
                                         detection_method=detection_method,
                                         metabolites=metabolites,
                                         de_genes_df=genes_df,
                                         tnseq_gene_list=tnseq_genes)

    else:
        bio_network = None

    if metabolites:
        mapping_msg = html.Div('''{} genes were mapped to the network out of {} genes in your list.\n{} 
                                   metabolites were mapped to these genes.'''
                               .format(len(bio_network.mapped_genes),
                                       len(gene_list),
                                       len(bio_network.mapped_metabolites)
                                       )
                               )
    else:
        mapping_msg = html.Div('{} genes were mapped to the network out of {} genes in your list.'
                               .format(len(bio_network.mapped_genes),
                                       len(gene_list))
                               )
    end_time = datetime.now()

    print("order = {}, detection_method = {}, metabolites = {}".format(order, detection_method, str(metabolites)))
    print(end_time - start_time)

    return bio_network, mapping_msg


@app.callback(
    Output('tnseq-gene-list-upload', 'style'),
    [Input('network-type', 'value')]
)
def show_tnseq_upload(network_type):
    """Show TnSeq upload button when combined networks are selected."""
    if network_type == 'combined':
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    Output('data-upload-output', 'children'),
    [Input('gene-list-upload', 'contents')],
    [State('gene-list-upload', 'filename')]
)
def upload_message(contents, filename):
    """Returns a successful message after file was uploaded."""
    if contents:
        try:
            small_table, genes_df = parse_gene_list(contents, filename)
            return small_table
        except ValueError:
            return html.Div('There was a problem uploading your file. Check that it is the correct format.')


@app.callback(
    [Output('make-network-message', 'children'),
     Output('download-link', 'href'),
     Output('hidden-bionetwork', 'children'),
     Output('network-parameters', 'children')],
    [Input('make-network', 'n_clicks')],
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
def update_download_link(
        n_clicks,
        network_type,
        strain,
        order,
        detection_method,
        metabolites,
        rnaseq_contents,
        tnseq_contents,
        rnaseq_filename,
        tnseq_filename):
    """Generates a network every time the make network button is clicked."""
    if n_clicks is None:
        raise PreventUpdate
    bio_network, mapping_msg = make_network(network_type,
                                                strain,
                                                order,
                                                detection_method,
                                                metabolites,
                                                rnaseq_filename,
                                                rnaseq_contents,
                                                tnseq_filename,
                                                tnseq_contents)
    rel_filename = os.path.join('downloads', '{}_network.graphml'.format(rnaseq_filename[:-4]))
    abs_filename = os.path.join(os.getcwd(), rel_filename)
    bio_network.write_gml(abs_filename)

    json_network = json.dumps(json_graph.node_link_data(bio_network.network))
    network_params = {'strain': bio_network._strain}

    return mapping_msg, '/{}'.format(rel_filename), json_network, json.dumps(network_params)


# Create and send downloadable file
@app.server.route('/downloads/<path:path>')
def download_graphml(path):
    root_dir = os.getcwd()
    downloads_dir = os.path.join(root_dir, 'downloads')
    if not os.path.exists(downloads_dir):
        os.mkdir(downloads_dir)

    return flask.send_from_directory(
        os.path.join(downloads_dir),
        path,
        cache_timeout=-1  # Prevent browser from caching previous file
    )




