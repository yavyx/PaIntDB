import base64
import io

import dash
from dash.dependencies import Output, Input, State
import dash_core_components as dcc
import dash_html_components as html
import dash_table
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
            {'label': 'All', 'value': 'all'},
            {'label': 'Experimental', 'value': 'experimental'},
            {'label': 'Mixed', 'value': 'mixed'},
            {'label': 'Computational', 'value': 'computational'},
        ],
        value='all'
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

    html.Div(id='make-network-output'),

    html.Hr(),

    html.Button('Download network (.graphml)')
])


def parse_gene_list(contents, filename):
    """Parses the uploaded gene list, returns a Dash table and a Pandas DataFrame."""
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
        small_df = df.head()
    except pd.errors.ParserError:
        return 'There was a problem uploading your file. Check that it is the correct format.'

    children = html.Div([
        html.Div('Your list was uploaded successfully!'),
        html.Div(filename),
        dash_table.DataTable(
            data=small_df.to_dict('records'),
            columns=[{'name': i, 'id': i} for i in small_df.columns]
        )
    ])
    return children, df


@app.callback(
    Output('data-upload-output', 'children'),
    [Input('gene-list-upload', 'contents')],
    [State('gene-list-upload', 'filename')]
)
def update_upload(contents, filename):
    """Returns a successful message after file was uploaded."""
    if contents is not None:
        children, df = parse_gene_list(contents, filename)
        return children


@app.callback(
    Output('make-network-output', 'children'),
    [Input('strain', 'value'),
     Input('order', 'value'),
     Input('detection-method', 'value'),
     Input('metabolites', 'value'),
     Input('gene-list-upload', 'contents')],
    [State('gene-list-upload', 'filename')]
)
def make_network(strain, order, detection_method, metabolites, contents, filename):
    """Generates a network after a file is uploaded."""
    if contents is not None:
        children, df = parse_gene_list(contents, filename)
        df = df.rename(columns={df.columns[0]: 'gene'})
        genes = list(df.gene)
        metabolites = True if metabolites else False

        bio_network = ng.BioNetwork(gene_list=genes,
                                    strain=strain,
                                    order=order,
                                    detection_method=detection_method,
                                    metabolites=metabolites)
        network = bio_network.get_network()
        return html.Div('{} genes were mapped to the network out of {} genes in your list.'.format(len(network),
                                                                                                   len(df.index)))


if __name__ == '__main__':
    app.run_server(debug=True)