import json

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash.dash import no_update
from networkx.readwrite import json_graph
import pandas as pd

from dash_app.app import app  # Loads app variable from app script
from dash_app.pages import home, menu, vis
from go_enrichment.go_enrichment import run_go_enrichment


app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    dbc.Navbar(
        id='top-bar',
        color='dark',
        dark=True,
        sticky='sticky',
        children=[
            dbc.Nav([
                dbc.NavbarBrand('PaintDB', href='/'),
                dbc.NavLink('Build Network', href='/menu'),
                dbc.NavLink('Explore Network', href='/vis', id='explore'),
                dbc.DropdownMenu(
                    children=[
                        dbc.DropdownMenuItem('Tutorial', href='#'),
                        dbc.DropdownMenuItem('About', href='#'),
                    ],
                    nav=True,
                    in_navbar=True,
                    label='More'
                )
            ])
        ],
    ),
    html.Div(id='page-content'),

    # Hidden divs to store and share data across callbacks and pages
    html.Div(id='hidden-bionetwork', style={'display': 'none'}),
    html.Div(id='network-parameters', style={'display': 'none'}),
    html.Div(id='node-details-df', style={'display': 'none'}),
    html.Div(id='enrichment-results', style={'display': 'none'}),
    html.Div(id='cyto-network', style={'display': 'none'}),
    html.Div(id='genes-of-interest', style={'display': 'none'}),
])


@app.callback(
    Output('explore', 'disabled'),
    [Input('hidden-bionetwork', 'children')]
)
def enable_explore_tab(bio_network):
    """Disables Explore tab if there is no network to explore."""
    return True if bio_network is None else False


def load_network(network_params, bio_network, genes_of_interest):
    """Loads the Bionetwork for use in the vis module."""
    network = json_graph.node_link_graph(json.loads(bio_network))

    cyto_network = dict()
    cyto_network['elements'], cyto_network['nodes'], cyto_network['edges'], network = \
        vis.make_cyto_elements(network, 5, 1000)

    network_params = json.loads(network_params)
    strain = network_params['strain']
    print('{} genes of interest'.format(len(json.loads(genes_of_interest))))
    enrichment_results, goea_results = run_go_enrichment(strain, json.loads(genes_of_interest))
    # Keep only overrepresented terms (remove underrepresented)
    enrichment_results = enrichment_results.loc[enrichment_results['enrichment'] == 'e', :]

    return json.dumps(cyto_network), enrichment_results.to_json()


@app.callback(
    [Output('page-content', 'children'),
     Output('cyto-network', 'children'),
     Output('enrichment-results', 'children')],
    [Input('url', 'pathname')],
    [State('hidden-bionetwork', 'children'),
     State('node-details-df', 'children'),
     State('network-parameters', 'children'),
     State('genes-of-interest', 'children')]
)
def display_page(pathname, bio_network, json_df, network_params, genes_of_interest):
    """Navigates to the selected app page. Generates vis layout depending on BioNetwork attributes."""
    if pathname == '/':
        return home.layout, no_update, no_update
    elif pathname == '/menu':
        return menu.layout, no_update, no_update
    elif pathname == '/vis':
        if bio_network:
            # Load JSON data
            network_df = pd.read_json(json_df)
            json_cyto_network, json_enrichment_results = load_network(network_params, bio_network, genes_of_interest)
            # Deserialize JSON
            cyto_network, enrichment_results = json.loads(json_cyto_network), pd.read_json(json_enrichment_results)
            # Generate layout using generated data
            vis_layout = vis.make_vis_layout(network_df, enrichment_results, cyto_network, network_params)
            return vis_layout, json_cyto_network, json_enrichment_results
        else:
            return html.Div('You need to create a network first.'), no_update, no_update
    else:
        return '404'


if __name__ == '__main__':
    app.run_server(debug=True)
