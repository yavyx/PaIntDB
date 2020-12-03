import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

from dash_app.app import app

layout = html.Div(
    [
        dbc.Jumbotron(
            [
                html.H1('Welcome to PaIntDB!'),
                dcc.Markdown(
                    '*Pseudomonas aeruginosa* Interactions Database',
                    className='lead',
                ),
                html.Hr(className='my-2'),
                dcc.Markdown(
                    'PaIntDB contains more than 157,000 protein-protein and protein-metabolite interactions in '
                    '*Pseudomonas aeruginosa* strains PAO1 and PA14.\n\n'
                    'It takes a list of significant genes identified through high-throughput experiments, maps the '
                    'interactions between them and returns a network that can be explored visually and filtered '
                    'to find putative biological pathways.\n\n'
                    '',
                    style={'width': '60vw',
                           'font-size': '20px'}
                ),
                dbc.Button('Get Started', color='primary', id='start', href='/menu', size='lg')
            ],
            style={'margin': '10px',
                   'backgroundColor': '#d1d1d1'}
        ),
        html.Div(
            [
                html.A(
                    html.Img(src=app.get_asset_url('hancock-lab-logo.svg'),
                             style={'width': '175px'}),
                    href='http://cmdr.ubc.ca/bobh/'
                ),
                dcc.Markdown('PaIntDB is being developed by the [Hancock Laboratory](http://cmdr.ubc.ca/bobh/) at the ' 
                             'University of British Columbia.\n\nFunding is currently provided by the Canadian '
                             'Institutes for Health Research FDN-154287.',
                             style={'margin-top': '10px'})
            ],
            style={'font-size': '14px',
                   'margin-left': '10px',
                   'padding': '7px'}
        )
    ],
    style={'background-color': '#ededed',
           'height': '95vh'}
)

