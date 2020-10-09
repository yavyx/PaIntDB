import json
import time

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash.dash import no_update
import networkx as nx
from networkx.algorithms import approximation
from networkx.readwrite import json_graph
import pandas as pd

from dash_app.app import app  # Loads app variable from app script
from dash_app.apps import app_home, app_menu, app_vis
from go_enrichment.go_enrichment import run_go_enrichment


app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    dbc.Navbar(
        id='top-bar',
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
        color='dark',
        dark=True,
        sticky='fixed'
    ),
    html.Div(id='page-content'),

    # Hidden divs to store and share data across callbacks
    html.Div(id='hidden-bionetwork', style={'display': 'none'}),
    html.Div(id='metric-closure', style={'display': 'none'}),
    html.Div(id='network-parameters', style={'display': 'none'}),
    html.Div(id='node-details-df', style={'display': 'none'}),
    html.Div(id='enrichment-results', style={'display': 'none'}),
    html.Div(id='cyto-network', style={'display': 'none'}),
    html.Div(id='genes-of-interest', style={'display': 'none'}),
    html.Div(id='hidden-div', style={'display': 'none'})
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
    all_nodes = len(network.nodes())

    # Extract and use just main component
    nodes = [component for component in nx.connected_components(network)][0]
    main_comp_network = network.subgraph(nodes)
    new_nodes = len(main_comp_network.nodes())
    print('Lost {} nodes after keeping main component.'.format(all_nodes - new_nodes))

    print('Calculating metric closure')
    t_start = time.time()
    metric_closure = approximation.metric_closure(main_comp_network)
    json_metric_closure = json.dumps(nx.node_link_data(metric_closure))
    print(time.time() - t_start)

    cyto_network = dict()
    cyto_network['elements'], cyto_network['nodes'], cyto_network['edges'], network = \
        app_vis.make_cyto_elements(main_comp_network, 2, 1000)

    network_params = json.loads(network_params)
    strain = network_params['strain']
    print('{} genes of interest'.format(len(json.loads(genes_of_interest))))
    enrichment_results, goea_results = run_go_enrichment(strain, json.loads(genes_of_interest))

    return json.dumps(cyto_network), enrichment_results.to_json(), json_metric_closure


@app.callback(
    [Output('page-content', 'children'),
     Output('cyto-network', 'children'),
     Output('metric-closure', 'children'),
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
        return app_home.layout, no_update, no_update, no_update
    elif pathname == '/menu':
        return app_menu.layout, no_update, no_update, no_update
    elif pathname == '/vis':
        if bio_network:
            # Load JSON data
            network_df = pd.read_json(json_df)
            json_cyto_network, json_enrichment_results, json_metric_closure = \
                load_network(network_params, bio_network, genes_of_interest)
            # Deserialize JSON
            cyto_network, enrichment_results = json.loads(json_cyto_network), pd.read_json(json_enrichment_results)
            # Generate layout using data
            vis_layout = app_vis.make_vis_layout(network_df, enrichment_results, cyto_network, network_params)
            return vis_layout, json_cyto_network, json_metric_closure, json_enrichment_results
        else:
            return html.Div('You need to create a network first.'), no_update, no_update, no_update
    else:
        return '404'


if __name__ == '__main__':
    app.run_server(debug=True)
