import base64
import io
import logging
import sys

import dash
from dash.dependencies import Output, Input, State
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import flask
import networkx as nx
import pandas as pd
import pandas.errors


import bio_networks.network_generator as ng

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    html.H1('PaIntDB'),
    dcc.Upload(
        id='gene-list-upload',
        children=html.Button(
            id='upload-button',
            children='Upload Gene List')
    ),
    html.Div(id='data-upload-output'),
    html.Hr(),
    html.Div('Select the strain:'),
    dcc.RadioItems(
        id='strain',
        options=[
            {'label': 'PAO1', 'value': 'PAO1'},
            {'label': 'PA14', 'value': 'PA14'},
        ],
        value='PAO1'
    ),
    html.Br(),
    html.Div('Select the network order:'),
    dcc.RadioItems(
        id='order',
        options=[
            {'label': 'Zero-order', 'value': 0},
            {'label': 'First-order', 'value': 1},
        ],
        value=0
    ),
    html.Br(),
    html.Div('Select the interaction detection method:'),
    dcc.RadioItems(
        id='detection-method',
        options=[
            {'label': 'All', 'value': 3},
            {'label': 'Experimental', 'value': 2},
            {'label': 'Mixed', 'value': 1},
            {'label': 'Computational', 'value': 0},
        ],
        value=3
    ),
    html.Br(),
    dcc.Checklist(
        id='metabolites',
        options=[
            {'label': 'Include metabolites', 'value': 1}
        ],
        value=[]
    ),
    html.Br(),
    dcc.Loading(id='loading', children=html.Div(id='make-network-output'), type='dot'),
    html.Hr(),
    html.A('Download Network(GraphML)',
           id='download-link'
           )
])


def parse_gene_list(contents, filename):
    """Parses the uploaded gene list, returns a Dash table and a Pandas DataFrame."""
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        gene_list = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
    except pd.errors.ParserError:
        return 'There was a problem uploading your file. Check that it is the correct format.'

    small_df = gene_list.head()
    children = html.Div([
        html.Div('Your list was uploaded successfully!'),
        html.Br(),
        html.Div(filename),
        dash_table.DataTable(
            data=small_df.to_dict('records'),
            columns=[{'name': i, 'id': i} for i in small_df.columns]
        )
    ])
    return children, gene_list


def make_network(strain, order, detection_method, metabolites, contents, filename):
    """Generates a network."""
    upload_msg, gene_list = parse_gene_list(contents, filename)
    gene_list = gene_list.rename(columns={gene_list.columns[0]: 'gene'})
    genes = list(gene_list.gene)
    metabolites = True if metabolites else False
    bio_network = ng.BioNetwork(gene_list=genes,
                                strain=strain,
                                order=order,
                                detection_method=detection_method,
                                metabolites=metabolites)
    if metabolites:
        mapping_msg = html.Div('''{} genes were mapped to the network out of {} genes in your list.\n{} 
                                   metabolites were mapped to these genes.'''
                               .format(len(ng.BioNetwork.get_mapped_genes(bio_network)),
                                       len(gene_list.index),
                                       len(ng.BioNetwork.get_mapped_metabolites(bio_network))
                                       )
                               )
        return bio_network, mapping_msg

    else:
        mapping_msg = html.Div('{} genes were mapped to the network out of {} genes in your list.'
                               .format(len(ng.BioNetwork.get_mapped_genes(bio_network)),
                                       len(gene_list.index))
                               )
        return bio_network, mapping_msg


@app.callback(
    Output('data-upload-output', 'children'),
    [Input('gene-list-upload', 'contents')],
    [State('gene-list-upload', 'filename')]
)
def upload_message(contents, filename):
    """Returns a successful message after file was uploaded."""
    if contents is not None:
        try:
            children, df = parse_gene_list(contents, filename)
            return children
        except ValueError:
            return html.Div('There was a problem uploading your file. Check that it is the correct format.')


@app.callback(
    Output('make-network-output', 'children'),
    [Input('strain', 'value'),
     Input('order', 'value'),
     Input('detection-method', 'value'),
     Input('metabolites', 'value'),
     Input('gene-list-upload', 'contents')],
    [State('gene-list-upload', 'filename')]
)
def network_message(strain, order, detection_method, metabolites, contents, filename):
    """Returns mapping message after file is uploaded or network parameters changed."""
    if contents is not None:
        network, mapping_msg = make_network(strain, order, detection_method, metabolites, contents, filename)
        return mapping_msg


@app.callback(
    Output('download-link', 'href'),
    [Input('strain', 'value'),
     Input('order', 'value'),
     Input('detection-method', 'value'),
     Input('metabolites', 'value'),
     Input('gene-list-upload', 'contents')],
    [State('gene-list-upload', 'filename')]
)
def update_link(strain, order, detection_method, metabolites, contents, filename):
    if contents is not None:
        return ('/dash/urlToDownload?strain={}&order={}&detection_method={}'
                '&metabolites={}&contents={}&filename={}').format(strain,
                                                                    order,
                                                                    detection_method,
                                                                    metabolites,
                                                                    contents,
                                                                    filename)


@app.server.route('/dash/urlToDownload')
def download_graphml():
    params = ['strain', 'order', 'detection_method', 'metabolites', 'contents', 'filename']
    args = [flask.request.args.get(param) for param in params]
    args[1] = int(args[1])
    args[2] = int(args[2])
    args[3] = True if args[3] else False
    bio_network, mapping_msg = make_network(*args)
    network = bio_network.get_network()
    str_io = io.StringIO()
    print('miculo')
    gml_string = ('\n'.join(nx.generate_graphml(network)))
    str_io.write(gml_string)
    mem = io.BytesIO()
    mem.write(str_io.getvalue().encode('utf-8'))
    mem.seek(0)
    return flask.send_file(mem,
                           mimetype='application/xml',
                           attachment_filename='network.graphml',
                           as_attachment=True)

# @app.route('/print')
# def printMsg():
#     app.logger.warning('testing warning log')
#     app.logger.error('testing error log')
#     app.logger.info('testing info log')
#     return "Check your console"

if __name__ == '__main__':
    app.run_server(debug=True)
