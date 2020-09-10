import json

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import networkx as nx
from networkx.readwrite import json_graph
import pandas as pd

from dash_app.app import app  # Loads app variable from app script
from dash_app.apps import app_home, app_menu, app_vis
import go_enrichment.go_enrichment as goe


app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    dbc.Navbar(
        id='top-bar',
        children=[
            dbc.Nav([
                dbc.NavbarBrand('PaintDB', href='/home'),
                dbc.NavLink("Build Network", href='/menu'),
                dbc.NavLink("Explore Network", href='/vis', id='explore'),
                dbc.DropdownMenu(
                    children=[
                        dbc.DropdownMenuItem("Tutorial", href='#'),
                        dbc.DropdownMenuItem("About", href="#"),
                    ],
                    nav=True,
                    in_navbar=True,
                    label="More"
                )
            ])
        ],
        color="dark",
        dark=True,
        sticky='fixed'
    ),
    html.Div(id='page-content'),

    # Hidden divs to store and share data across callbacks
    html.Div(id='hidden-bionetwork', style={'display': 'none'}),
    html.Div(id='network-parameters', style={'display': 'none'}),
    html.Div(id='node-details-df', style={'display': 'none'}),
    html.Div(id='enrichment-results', style={'display': 'none'}),
    html.Div(id='hidden-div', style={'display': 'none'})
])


@app.callback(
    Output('explore', 'disabled'),
    [Input('hidden-bionetwork', 'children')]
)
def enable_explore_tab(bio_network):
    return True if bio_network is None else False


def load_network(network_params, bio_network):
    "Loads the Bionetwork for use in the vis module."
    network = json_graph.node_link_graph(json.loads(bio_network))
    all_nodes = len(network.nodes())

    # Extract and use just main component
    nodes = [component for component in nx.connected_components(network)][0]
    main_comp_network = network.subgraph(nodes)
    new_nodes = len(main_comp_network.nodes())
    print('Lost {} nodes after keeping main component.'.format(all_nodes - new_nodes))

    cyto_elements, cyto_nodes, cyto_edges, network_main_comp = app_vis.make_cyto_elements(main_comp_network)
    network_df = app_vis.make_network_df(main_comp_network)

    network_params = json.loads(network_params)
    strain = network_params['strain']
    enrichment_results, goea_results = goe.run_go_enrichment(strain, list(network_df.index))

    return [cyto_elements, network_df.to_json(), enrichment_results.to_json()]


@app.callback(
    Output('page-content', 'children'),
    [Input('url', 'pathname')],
    [State('hidden-bionetwork', 'children'),
     State('node-details-df', 'children'),
     State('enrichment-results', 'children'),
     State('network-parameters', 'children')]
)
def display_page(pathname, bio_network, network_df, enrichment_results, network_params):
    """Navigates to the selected app page. Generates vis layout depending on BioNetwork attributes."""
    if pathname == '/home':
        return app_home.layout
    elif pathname == '/menu':
        return app_menu.layout
    elif pathname == '/vis':
        if bio_network is not None:
            elements, json_df, json_enrichment_results = load_network(network_params, bio_network)
            network_df, enrichment_results = pd.read_json(json_df), pd.read_json(json_enrichment_results)
            vis_layout = app_vis.make_vis_layout(network_df, enrichment_results, elements)
            return vis_layout
        else:
            return html.Div('You need to create a network first.')
    else:
        return '404'


if __name__ == '__main__':
    app.run_server(debug=True)
