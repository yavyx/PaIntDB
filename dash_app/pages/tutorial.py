import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

layout = html.Div(
    [
        dbc.Jumbotron(
            [
                html.H1('Tutorial'),
                dcc.Markdown(
                    'Generating Networks with PaIntDB',
                    className='lead',
                ),
                html.Hr(),
                html.H3('Step 1: Data Upload'),
                dcc.Markdown(
                    '',
                    style={'width': '60vw',
                           'font-size': '20px'}
                ),
                dbc.Button('Get Started', color='primary', id='start', href='/menu', size='lg')
            ],
            style={'margin': '10px',
                   'backgroundColor': '#d1d1d1'}
        ),
    ]
)