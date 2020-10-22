from dash.dependencies import Output, Input, State
import dash_bootstrap_components as dbc
import dash_html_components as html

layout = dbc.Container(
    [
        dbc.Row(
            dbc.Col(
                html.H2('Welcome to PaIntDB!')
            )
        )
    ]
)
