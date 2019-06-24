
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from app import app
from apps import D2420, D2460,BplusD0,PhaseSpace
import flask

import dash

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
])



@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/apps/D2460':
         return D2460.layout
    elif pathname == '/apps/D2420':
         return D2420.layout
    elif pathname == '/apps/Dst':
         return PhaseSpace.layout
    elif pathname == '/apps/D0':
         return BplusD0.layout

    else:
        return '404'


"""
if __name__ == '__main__':
    app.run_server(debug=True)


app = dash.Dash(
    __name__, 
    external_stylesheets=['https://codepen.io/chriddyp/pen/bWLwgP.css']
)

url_bar_and_content_div = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
])


layout_index = html.Div([
    dcc.Link('D2420', href='/apps/D2420'),
    html.Br(),
    dcc.Link('D2460', href='/apps/D2460'),
])



layout_page_1=D2420.layout
layout_page_2=D2460.layout


def serve_layout():
    if flask.has_request_context():
        return url_bar_and_content_div
    return html.Div([
        url_bar_and_content_div,
        layout_index,
        layout_page_1,
        layout_page_2,
    ])

app.layout = serve_layout



# Index callbacks
@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == "/apps/D2420":
        return layout_page_1
    elif pathname == "/apps/D2460":
        return layout_page_2
    else:
        return layout_index

"""

if __name__ == '__main__':
    app.run_server(debug=True)

