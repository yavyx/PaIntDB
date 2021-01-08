import json

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash.dash import no_update
from networkx.readwrite import json_graph
import pandas as pd

from dash_app.app import app, server
from dash_app.pages import home, menu, vis, tutorial

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    dbc.Navbar(
        id='top-bar',
        color='#1a3775',
        sticky='sticky',
        style={'height': '65px'},
        children=[
            dbc.Nav(
                dbc.Row(
                    [
                        html.A(
                            html.Img(src=app.get_asset_url('PaintDB-logo-small.png'),
                                     id='legend',
                                     height='60px',
                                     ),
                            href='/'),
                        dbc.NavLink('Build Network', href='/menu', active='exact'),
                        dbc.NavLink('Explore Network', href='/vis', id='explore', disabled=True),
                        dbc.DropdownMenu(
                            children=[
                                dbc.DropdownMenuItem('Tutorial', href='/tutorial'),
                                dbc.DropdownMenuItem('About', href='/about'),
                            ],
                            nav=True,
                            in_navbar=True,
                            label='More'
                        ),
                    ],
                    align='center')
            )
        ],
    ),
    # Different pages are shown here
    html.Div(id='page-content'),

    # Hidden divs to store and share JSON data across callbacks and pages
    html.Div(id='hidden-bionetwork', style={'display': 'none'}),
    html.Div(id='network-parameters', style={'display': 'none'}),
    html.Div(id='node-details-df', style={'display': 'none'}),
    html.Div(id='enrichment-results', style={'display': 'none'}),
    html.Div(id='cyto-network', style={'display': 'none'}),
    html.Div(id='genes-of-interest', style={'display': 'none'}),
])


@app.callback(
    Output('explore', 'disabled'),
    [Input('hidden-bionetwork', 'children'),
     Input('enrichment-results', 'children')]
)
def enable_explore_tab(bio_network, enrichment=None):
    """Disables Explore tab if there is no network to explore."""
    return False if bio_network and enrichment else True


@app.callback(
    [Output('page-content', 'children'),
     Output('cyto-network', 'children')],
    [Input('url', 'pathname')],
    [State('hidden-bionetwork', 'children'),
     State('node-details-df', 'children'),
     State('enrichment-results', 'children'),
     State('network-parameters', 'children')]
)
def display_page(pathname, bio_network, json_df, json_enrichment_results, network_params):
    """Navigates to the selected app page. Generates vis layout depending on BioNetwork type."""
    if pathname == '/':
        return home.layout, no_update
    elif pathname == '/menu':
        return menu.layout, no_update
    elif pathname == '/vis':
        if bio_network:
            # Load JSON data
            network = json_graph.node_link_graph(json.loads(bio_network))
            network_df = pd.read_json(json_df)
            enrichment_results = pd.read_json(json_enrichment_results)
            # Create cyto network for visualization
            cyto_network = vis.make_cyto_elements(network)
            # Generate layout using generated data
            vis_layout = vis.make_vis_layout(network_df, enrichment_results, cyto_network, network_params)
            return vis_layout, json.dumps(cyto_network)
        else:
            return dbc.Alert('Warning: You need to build a network first.',
                             color='warning',
                             style={'display': 'inline-block', 'margin': '10px'}), no_update
    elif pathname == '/tutorial':
        return tutorial.layout, no_update
    elif pathname == '/about':
        return dbc.Alert('Coming Soon.',
                         color='primary',
                         style={'display': 'inline-block', 'margin': '10px'}), no_update
    else:
        return '404'


if __name__ == '__main__':
    app.run_server(debug=True)
