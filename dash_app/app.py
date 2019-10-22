import base64
import io
import urllib.parse
from datetime import datetime

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

radio_buttons_style = {'width': '20%',
                       'display': 'inline-block',
                       'verticalAlign': 'top'}

app.layout = html.Div([
    html.H1('PaIntDB'),
    html.Div(['Select your input data: ',
              dcc.RadioItems(
                  id='network-type',
                  options=[
                      {'label': 'Gene List', 'value': 'basic'},
                      {'label': 'Differential Expression Gene List', 'value': 'DE'},
                      {'label': 'Combined (DE genes and TnSeq genes)', 'value': 'combined'}
                  ],
                  value='basic'),
             html.Br(),
             html.P('Your genes must be in the first column of a CSV file.'),
             dcc.Upload(
                 id='gene-list-upload',
                 children=html.Button(
                     id='upload-button',
                     children='Upload Gene List'))
              ],
             style={'width': '25%', 'display': 'inline-block'},
             ),
    html.Div(id='data-upload-output',
             style={'height': '30vh',
                    'display': 'inline-block',
                    'verticalAlign': 'top'}),
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
             style=radio_buttons_style),
    html.Div(['Select the network order: ',
              html.Abbr("?", title=('Zero-order: Maps direct interactions between your queried genes. '
                                    'Recommended for long lists (>200 genes).\n\n'
                                    'First-order: Uses your queried genes as "seed" genes and finds any interaction '
                                    'between them and the other genes in the database. '
                                    'Recommended for short lists (<200 genes).')),
              dcc.RadioItems(
                  id='order',
                  options=[
                      {'label': 'Zero-order', 'value': 0},
                      {'label': 'First-order', 'value': 1},
                  ],
                  value=0
              )],
             style=radio_buttons_style),
    html.Div(['Select the interaction detection method: ',
              html.Abbr("?", title=('Choose which interactions you want to use to generate the network.\n\n'
                                    'Experimentally-verified interactions have the highest confidence, but '
                                    'result in smaller networks.\n\n '
                                    'If mixed, the detection method is ambiguous.')),
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
             style=radio_buttons_style),
    html.Br(),
    dcc.Checklist(
        id='metabolites',
        options=[
            {'label': 'Include metabolites', 'value': 1}
        ]
    ),
    html.Br(),
    dcc.Loading(id='loading',
                children=html.Div(id='make-network-message'),
                type='dot'),
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
        genes_df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
    except pd.errors.ParserError:
        return 'There was a problem uploading your file. Check that it is the correct format.'

    small_df = genes_df.head()
    children = [
        html.P('Your list was uploaded successfully!'),
        html.Br(),
        html.P(filename),
        dash_table.DataTable(
            data=small_df.to_dict('records'),
            columns=[{'name': i, 'id': i} for i in small_df.columns]
        )]
    return children, genes_df


def make_network(network_type, strain, order, detection_method, metabolites, contents, filename):
    """Generates a BioNetwork."""
    start_time = datetime.now()
    upload_msg, genes_df = parse_gene_list(contents, filename)
    genes_df.rename(columns={genes_df.columns[0]: 'gene'}, inplace=True)
    gene_list = list(genes_df.gene)
    metabolites = True if metabolites else False
    bio_network = None
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
    # elif network_type == 'combined':
    #     bio_network = ng.CombinedNetwork()
    else:
        bio_network = ''

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
    end_time = datetime.now()
    # with open('performance.txt', 'a') as f:
    #     f.write("order = {}, detection_method = {}, metabolites = {}\n".format(order, detection_method, str(metabolites)))
    #     f.write('{}\n\n'.format(end_time - start_time))
    print("order = {}, detection_method = {}, metabolites = {}".format(order, detection_method, str(metabolites)))
    print(end_time - start_time)
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
def update_download_link(network_type, strain, order, detection_method, metabolites, contents, filename):
    """Generates a network every time parameters are changed."""
    if contents is not None:
        bio_network, mapping_msg = make_network(network_type, strain, order, detection_method, metabolites, contents,
                                                filename)
        nx_network = bio_network.get_network()
        url = '/dash/download?gml={}'.format(('\n'.join(nx.generate_graphml(nx_network))))  # Create URL string
    else:
        mapping_msg, url = [None, None]  # Workaround for now (Callback State not functioning properly)
    return mapping_msg, url


@app.server.route('/dash/download')
def download_graphml():
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
