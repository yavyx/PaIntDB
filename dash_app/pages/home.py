import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

layout = dbc.Jumbotron(
    [
        html.H3('Welcome to PaIntDB!', className='display-3'),
        dcc.Markdown(
            'Pseudomonas aeruginosa Interactions Database',
            className='lead',
        ),
        html.Hr(className='my-2'),
        dcc.Markdown(
            'PaIntDB contains more than 157K protein-protein and protein-metabolite interactions in *Pseudomonas '
            'aeruginosa* strains PAO1 and PA14.'
            '\n\nIt takes high-throughput (RNASeq and/or TnSeq) experimental results as input, maps the interactions '
            'between the significant genes, and returns a network that can be explored visually and analyzed further.',
            style={'width': '60vw'}
        ),
        html.P(dbc.Button('Get Started', color='primary', id='start', href='/menu'), className='lead'),
    ]
)

