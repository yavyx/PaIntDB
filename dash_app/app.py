import base64
import io

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
    html.Div(['Select the network type:',
              dcc.RadioItems(
                  id='network-type',
                  options=[
                      {'label': 'No experimental info', 'value': 'basic'},
                      {'label': 'Differential Expression', 'value': 'DE'},
                  ],
                  value='basic'
              )]),
    html.Br(),
    html.Div([
        html.P('Your genes must be in the first column of a CSV file.'),
        dcc.Upload(
            id='gene-list-upload',
            children=html.Button(
                id='upload-button',
                children='Upload Gene List')
        )],
        style={'width': '29%', 'display': 'inline-block'}
    ),

    html.Div(id='data-upload-output',
             style={'width': '69%', 'display': 'inline-block'}),
    html.Hr(),
    html.Div(['Select the strain:',
              dcc.RadioItems(
                  id='strain',
                  options=[
                      {'label': 'PAO1', 'value': 'PAO1'},
                      {'label': 'PA14', 'value': 'PA14'},
                  ],
                  value='PAO1'
              )],
             style={'width': '20%', 'display': 'inline-block', 'vertical-align': 'top'}),
    html.Div(['Select the network order:',
              dcc.RadioItems(
                  id='order',
                  options=[
                      {'label': 'Zero-order', 'value': 0},
                      {'label': 'First-order', 'value': 1},
                  ],
                  value=0
              )],
             style={'width': '20%', 'display': 'inline-block', 'vertical-align': 'top'}),
    html.Div(['Select the interaction detection method:',
              dcc.RadioItems(
                  id='detection-method',
                  options=[
                      {'label': 'All', 'value': 3},
                      {'label': 'Experimental', 'value': 2},
                      {'label': 'Mixed', 'value': 1},
                      {'label': 'Computational', 'value': 0},
                  ],
                  value=3
              )],
             style={'width': '20%', 'display': 'inline-block', 'vertical-align': 'top'}),
    html.Br(),
    dcc.Checklist(
        id='metabolites',
        options=[
            {'label': 'Include metabolites', 'value': 1}
        ]
    ),
    html.Br(),
    dcc.Loading(id='loading', children=html.Div(id='make-network-message'), type='dot'),
    html.Hr(),
    html.Button(html.A('Download Network(GraphML)',
                       id='download-link'
                       )
                )
])


def parse_gene_list(contents, filename):
    """Parses the uploaded gene list, returns a Dash table and a Pandas DataFrame."""
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        genes_df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
    except pd.errors.ParserError:
        return 'There was a problem uploading your file. Check that it is the correct format.'

    small_df = genes_df.head()
    small_table = html.Div([
        html.Div('Your list was uploaded successfully!'),
        html.Br(),
        html.Div(filename),
        dash_table.DataTable(
            data=small_df.to_dict('records'),
            columns=[{'name': i, 'id': i} for i in small_df.columns]
        )
    ])
    return small_table, genes_df


def make_network(network_type, strain, order, detection_method, metabolites, contents, filename):
    """Generates a network."""
    upload_msg, genes_df = parse_gene_list(contents, filename)
    genes_df.rename(columns={genes_df.columns[0]: 'gene'}, inplace=True)
    gene_list = list(genes_df.gene)
    metabolites = True if metabolites else False
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
                                   genes_df=genes_df)
    else:
        bio_network = 'miguebo'

    if metabolites:
        mapping_msg = html.Div('''{} genes were mapped to the network out of {} genes in your list.\n{} 
                                   metabolites were mapped to these genes.'''
                               .format(len(ng.BioNetwork.get_mapped_genes(bio_network)),
                                       len(gene_list),
                                       len(ng.BioNetwork.get_mapped_metabolites(bio_network))
                                       )
                               )
    else:
        mapping_msg = html.Div('{} genes were mapped to the network out of {} genes in your list.'
                               .format(len(ng.BioNetwork.get_mapped_genes(bio_network)),
                                       len(gene_list))
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
            small_table, genes_df = parse_gene_list(contents, filename)
            return small_table
        except ValueError:
            return html.Div('There was a problem uploading your file. Check that it is the correct format.')


@app.callback(
    [Output('make-network-message', 'children'),
     Output('download-link', 'href')],
    [Input('network-type', 'value'),
     Input('strain', 'value'),
     Input('order', 'value'),
     Input('detection-method', 'value'),
     Input('metabolites', 'value'),
     Input('gene-list-upload', 'contents')],
    [State('gene-list-upload', 'filename')]
)
def update_link(network_type, strain, order, detection_method, metabolites, contents, filename):
    """Generates a network and updates the download link every time parameters are changed."""
    if contents is not None:
        bio_network, mapping_msg = make_network(network_type, strain, order, detection_method, metabolites, contents,
                                                filename)
        nx_network = bio_network.get_network()
        gml_string = ('\n'.join(nx.generate_graphml(nx_network)))
        download_url = '/dash/download?gml={}'.format(gml_string)
        #download_url = '/download/<{}>'.format(network)
    else:
        mapping_msg, download_url = [None, None]  # Workaround for now (Callback State not functioning properly)
    return mapping_msg, download_url


@app.server.route('/dash/download')
def download_graphml():
    # params = ['strain', 'order', 'detection_method', 'metabolites', 'contents', 'filename']
    # args = [flask.request.args.get(param) for param in params]
    # args[1] = int(args[1])  # Change order to integer (Flask returns string)
    # args[2] = int(args[2])  # Change detection method to integer (Flask returns string)
    # args[3] = True if args[3] is 1 else False  # Change metabolites to True/False
    # bio_network, mapping_msg = make_network(*args)
    gml_string = flask.request.args.get('gml')
    str_io = io.StringIO()
    str_io.write(gml_string)
    mem = io.BytesIO()
    mem.write(str_io.getvalue().encode('utf-8'))
    mem.seek(0)
    return flask.send_file(mem,
                           mimetype='application/xml',
                           attachment_filename='network.graphml',
                           as_attachment=True)


if __name__ == '__main__':
    app.run_server(debug=True)
