import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from dash_app.app import app  # Loads app variable from app script
from dash_app.apps import app_menu, app_vis


app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    dbc.Navbar(
        id='top-bar',
        children=[
            dbc.Nav([
                dbc.NavbarBrand('PaintDB', href='/menu'),
                dbc.NavLink("Build Network", href="/menu"),
                dbc.NavLink("Explore Network", href="/vis"),
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
        dark=True
    ),
    html.Div(id='page-content')
])


@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/home':
        return app_menu.layout
    elif pathname == '/menu':
        return app_menu.layout
    elif pathname == '/vis':
        return app_vis.layout
    else:
        return '404'


if __name__ == '__main__':
    app.run_server(debug=True)
