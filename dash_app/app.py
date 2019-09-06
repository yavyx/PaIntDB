import base64
import io

import dash
from dash.dependencies import Output, Input
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import pandas as pd

import bio_networks.network_generator as ng

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    html.H1('PaIntDB'),

    dcc.Upload(html.Button('Upload Gene List'), id='gene-list-upload'),

    html.Div(id='output-data-upload'),

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

    html.Div('Select the network order:'),
    dcc.RadioItems(
        id='order',
        options=[
            {'label': 'Zero-order', 'value': 0},
            {'label': 'First-order', 'value': 1},
        ],
        value=0
    ),

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

    html.Div('Include metabolites?'),
    dcc.RadioItems(
        id='metabolites',
        options=[
            {'label': 'Yes', 'value': 'True'},
            {'label': 'No', 'value': 'False'},
        ],
        value='False'
    ),

    html.Button('Make network'),

    html.Hr(),

    html.Button('Download network (.graphml)')
])


def parse_gene_list(contents):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
    return html.Div([
        dash_table.DataTable(
            data=df.to_dict('records'),
            columns=[{'name': i, 'id': i} for i in df.columns]
        )
    ])


@app.callback(
    Output('output-data-upload', 'children'),
    [Input('gene-list-upload', 'contents')]
)
def update_output(contents):
    children = parse_gene_list(contents)
    return children


if __name__ == '__main__':
    app.run_server(debug=True)