import dash
import dash_core_components as dcc
import dash_html_components as html


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    html.H1(children='PaIntDB'),

    html.Div(children='Select the strain:'),
    dcc.RadioItems(
        options=[
            {'label': 'PAO1', 'value': 'PAO1'},
            {'label': 'PA14', 'value': 'PA14'},
        ],
        value='PAO1'
    ),

    html.Div(children='Select the network order:'),
    dcc.RadioItems(
        options=[
            {'label': 'Zero-order', 'value': 0},
            {'label': 'First-order', 'value': 1},
        ],
        value=0
    ),

    html.Div(children='Select the interaction detection method:'),
    dcc.RadioItems(
        options=[
            {'label': 'All', 'value': 'all'},
            {'label': 'Experimental', 'value': 'experimental'},
            {'label': 'Mixed', 'value': 'mixed'},
            {'label': 'Computational', 'value': 'computational'},
        ],
        value='all'
    ),

    html.Div(children='Include metabolites?'),
    dcc.RadioItems(
        options=[
            {'label': 'Yes', 'value': 'True'},
            {'label': 'No', 'value': 'False'},
        ],
        value='False'
    ),
])

if __name__ == '__main__':
    app.run_server(debug=True)