import dash

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__)
app.css.append_css({'external_url': 'https://codepen.io/plotly/pen/EQZeaW.css'})
server = app.server




app.config.suppress_callback_exceptions = True
