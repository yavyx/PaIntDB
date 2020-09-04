import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash
from dash.dependencies import Input, Output, State

import json

from dash_app.app import app  # Loads app variable from app script
from dash_app.apps import app_home, app_menu, app_vis


app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    dbc.Navbar(
        id='top-bar',
        children=[
            dbc.Nav([
                dbc.NavbarBrand('PaintDB', href='/home'),
                dbc.NavLink("Build Network", href="/menu"),
                dbc.NavLink("Explore Network", href="/vis", id='explore'),
                dbc.DropdownMenu(
                    children=[
                        dbc.DropdownMenuItem("Help", header=True),
                        dbc.DropdownMenuItem("Page 2", href="#"),
                        dbc.DropdownMenuItem("Page 3", href="#")
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
    html.Div(id='hidden-bionetwork', style={'display': 'none'}),
    html.Div(id='network-parameters', style={'display': 'none'}),
    html.Div(id='page-content')
])


@app.callback(
    Output('explore', 'disabled'),
    [Input('hidden-bionetwork', 'children')]
)
def enable_explore_tab(bio_network):
    if bio_network is None:
        return True


@app.callback(
    Output('page-content', 'children'),
    [Input('url', 'pathname')],
    [State('hidden-bionetwork', 'children')]
)
def display_page(pathname, bio_network):
    """Navigates to the selected app page."""
    if pathname == '/home':
        return app_home.layout
    elif pathname == '/menu':
        return app_menu.layout
    elif pathname == '/vis':
        if bio_network is not None:
            return app_vis.layout
        else:
            return html.Div('You need to create a network first.')
    else:
        return '404'


# def get_layout(bio_network):
#     """Adds hidden divs to generate a vis layout depending on the BioNetwork."""
#     return [
#         ,
#         app_vis.layout
#     ]

if __name__ == '__main__':
    app.run_server(debug=True)
