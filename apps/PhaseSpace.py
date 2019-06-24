



"""
Created on Wed Feb 20 18:14:41 2019

@author: ke
"""


from app import app
from textwrap import dedent as d
import json
import dash
from dash.dependencies import Input, Output,State
import dash_html_components as html
import dash_core_components as dcc
import pandas as pd
import os
import root_pandas
import pandas as pd
import hepvector
from hepvector.numpyvector import Vector3D,LorentzVector

from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from numpy import cos,sin,tan,sqrt,absolute,real,conjugate,imag,abs,max,min

import hepvector
from hepvector.numpyvector import Vector3D,LorentzVector

import plotly
import plotly.graph_objs as go
import plotly.plotly as py 
import plotly.tools as tls
from plotly.graph_objs import Data, Layout, Figure
from plotly.graph_objs import Scatter



"""
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__)#, external_stylesheets=external_stylesheets)
app.css.append_css({'external_url': 'https://codepen.io/plotly/pen/EQZeaW.css'})
server = app.server
"""


dft=root_pandas.read_root('/home/ke/tests/apps/model_tree.root',key='DecayTree')




styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}
dft['W_PX_TRUE']=dft['B_PX_TRUE']-dft['Dst_PX_TRUE']
dft['W_PY_TRUE']=dft['B_PY_TRUE']-dft['Dst_PY_TRUE']
dft['W_PZ_TRUE']=dft['B_PZ_TRUE']-dft['Dst_PZ_TRUE']
dft['W_E_TRUE']=dft['B_E_TRUE']-dft['Dst_E_TRUE']
df=dft.head(50000)

B=LorentzVector(df['B_PX_TRUE'],df['B_PY_TRUE'],df['B_PZ_TRUE'],df['B_E_TRUE'])

W=LorentzVector(df['W_PX_TRUE'],df['W_PY_TRUE'],df['W_PZ_TRUE'],df['W_E_TRUE'])

Dst=LorentzVector(df['Dst_PX_TRUE'],df['Dst_PY_TRUE'],df['Dst_PZ_TRUE'],df['Dst_E_TRUE'])
tau=LorentzVector(df['Tau_PX_TRUE'],df['Tau_PY_TRUE'],df['Tau_PZ_TRUE'],df['Tau_E_TRUE'])
D0=LorentzVector(df['D0_PX_TRUE'],df['D0_PY_TRUE'],df['D0_PZ_TRUE'],df['D0_E_TRUE'])
nuB=LorentzVector(df['B_nu_PX_TRUE'],df['B_nu_PY_TRUE'],df['B_nu_PZ_TRUE'],df['B_nu_E_TRUE'])
K=LorentzVector(df['D0_K_PX_TRUE'],df['D0_K_PY_TRUE'],df['D0_K_PZ_TRUE'],df['D0_K_E_TRUE'])
piDst=LorentzVector(df['Dst_Pi_PX_TRUE'],df['Dst_Pi_PY_TRUE'],df['Dst_Pi_PZ_TRUE'],df['Dst_Pi_E_TRUE'])
piK=LorentzVector(df['D0_Pi_PX_TRUE'],df['D0_Pi_PY_TRUE'],df['D0_Pi_PZ_TRUE'],df['D0_Pi_E_TRUE'])
pitau1=LorentzVector(df['Tau_Pi1_PX_TRUE'],df['Tau_Pi1_PY_TRUE'],df['Tau_Pi1_PZ_TRUE'],df['Tau_Pi1_E_TRUE'])
pitau2=LorentzVector(df['Tau_Pi2_PX_TRUE'],df['Tau_Pi2_PY_TRUE'],df['Tau_Pi2_PZ_TRUE'],df['Tau_Pi2_E_TRUE'])
pitau3=LorentzVector(df['Tau_Pi3_PX_TRUE'],df['Tau_Pi3_PY_TRUE'],df['Tau_Pi3_PZ_TRUE'],df['Tau_Pi3_E_TRUE'])
nutau=LorentzVector(df['Tau_nu_PX_TRUE'],df['Tau_nu_PY_TRUE'],df['Tau_nu_PZ_TRUE'],df['Tau_nu_E_TRUE'])




def change_frame(COM):
	newB=B.boost(-COM.boostp3)
	newDst=Dst.boost(-COM.boostp3)
	newtau=tau.boost(-COM.boostp3)
	newD0=D0.boost(-COM.boostp3)
	newnuB=nuB.boost(-COM.boostp3)
	newK=K.boost(-COM.boostp3)
	newpiDst=piDst.boost(-COM.boostp3)
	newpitau1=pitau1.boost(-COM.boostp3)
	newpitau2=pitau2.boost(-COM.boostp3)
	newpitau3=pitau3.boost(-COM.boostp3)
	newnutau=nutau.boost(-COM.boostp3)
	newpiK=piK.boost(-COM.boostp3)
	res=[newB,newDst,newtau,newD0,newnuB,newK,newpiDst,newpitau1,newpitau2,newpitau3,newnutau,newpiK]
	return res




nouvtau=tau.boost(-(tau+nuB).boostp3)
nouvnu=nuB.boost(-(tau+nuB).boostp3)
nouvpi=piDst.boost(-(piDst+D0).boostp3)
nouvD0=D0.boost(-(piDst+D0).boostp3)
nouvDst=D0.boost(-B.boostp3)
unittau=(nouvtau.p3).unit
unitnu=(nouvnu.p3).unit
unitDst=(nouvDst.p3).unit
unitD0=(nouvD0.p3).unit
nnewtau=tau.boost(-B.boostp3)
nnewD0=D0.boost(-B.boostp3)
unitau=nnewtau.unit
uniD0=nnewD0.unit
nnormal1=unitDst.cross(uniD0)
normal1=nnormal1.unit
nnormal2=unitDst.cross(unitau)
normal2=nnormal2.unit
pparallel=normal1.cross(unitDst)
parallel=pparallel.unit

co = normal1.dot(normal2)
si = parallel.dot(normal2)




chi = np.arctan2(si,co)



costhetast=unitD0.dot(unitDst)
costhetal=unitDst.dot(unittau)



thetast=np.arccos(costhetast)
thetal=np.arccos(costhetal)





nouvtau=tau.boost(-(tau+nuB).boostp3)
nouvnu=nuB.boost(-(tau+nuB).boostp3)
nouvpi=piDst.boost(-(piDst+D0).boostp3)
nouvD0=D0.boost(-(piDst+D0).boostp3)
nouvDst=D0.boost(-B.boostp3)
unittau=(nouvtau.p3).unit
unitnu=(nouvnu.p3).unit
unitDst=(nouvDst.p3).unit
unitD0=(nouvD0.p3).unit
nnewtau=tau.boost(-B.boostp3)
nnewD0=D0.boost(-B.boostp3)
unitau=nnewtau.unit
uniD0=nnewD0.unit
nnormal1=unitDst.cross(uniD0)
normal1=nnormal1.unit
nnormal2=unitDst.cross(unitau)
normal2=nnormal2.unit
pparallel=normal1.cross(unitDst)
parallel=pparallel.unit

co = normal1.dot(normal2)
si = parallel.dot(normal2)




chi = np.arctan2(si,co)



costhetast=unitD0.dot(unitDst)
costhetal=unitDst.dot(unittau)
trace_phase=go.Scatter3d(
	    x=costhetast,
	    y=costhetal,
	    z=chi,  
	    mode='markers',
	    marker=dict(
		size=5,
		color=W.mag2,
#		color=chi,                
		colorscale='Electric',  
		opacity=0.8,
		colorbar=dict(thickness=20)
	    )
)	



############################################################################         
layout = html.Div(children=[
	html.Div([
		html.Div([
			html.Div([
				html.H2('B to D* decay visualisation',
                		style={
                    		'position': 'relative',
                    		'top': '0px',
                    		'left': '10px',
                    		'font-family': 'Dosis',
                    		'display': 'inline',
				'font-size':'5.0rem',
                    		'color': '#4D637F'
               			 }),
               		 ]),
			html.Div([
                        	html.Br(),
                        	html.P('Choose ranges :',
                                style={
                                'font-family': 'Dosis',
                                'display': 'inline',     #block?
                                'font-size': '3rem',
                                'color': '#4D637F'
                                 })
				],style={'marginLeft':0}),
			html.Br(),            
        	    	html.P('costheta* : ',
				style={
					'display':'inline-block',
					'verticalAlign': 'top',
					'marginRight': '10px'
				}
           		 ),
            		html.Div([
				dcc.RangeSlider(
					id='choose-thetast_Dst',
					min=-1,
					max=1,
			        	value=[-1,1], step=0.1
				),
    	            		html.Div(id='output-range-thetast_Dst')
                 	 ],style={'width':300,'display':'inline-block','marginBottom':10,'marginLeft':8}),
      	       		html.Br(),
                        html.P('costheta_l :',
					style={
						'display':'inline-block',
						'verticalAlign': 'top',
						'marginRight': '10px'
					}
           	 	),
            		html.Div([
				dcc.RangeSlider(
					id='choose-thetal_Dst',
					min=-1, max=1, value=[-1,1], step=0.1,
				),
                		html.Div(id='output-range-thetal_Dst')
            		], style={'width':300, 'display':'inline-block', 'marginBottom':11}),

           		 html.Br(),
           		 html.P('chi : ',
					style={
						'display':'inline-block',
						'verticalAlign': 'top',
						'marginRight': '10px'
					}
           		 ),
           		 html.Div([
				dcc.RangeSlider(
					id='choose-chi_Dst',
					min=-3.14, max=3.14, value=[-3.14,3.14], step=0.1,
				),
              			html.Div(id='output-range-chi_Dst')
            		], style={'width':300, 'display':'inline-block', 'marginBottom':10,'marginLeft':46}),
                
        	], style={'margin':0} ),
        
       		html.P('Phase space of selected ranges',
			style={      
                    'font-family': 'Dosis',
                    'display': 'inline',
                    'font-size': '3.0rem',
                    'color': '#4D637F'
                }
        	),
        
        	dcc.Graph(id = 'phase-space_Dst'),	
                html.Div([
            dcc.Markdown(d("""
                **Click Data**

                Click on points in the graph.
            """)),
            html.Pre(id='click-result_Dst', style=styles['pre']),
        ], className='six columns'),    
        html.Br(),
	html.Br(),
	html.Br(),
	html.Br(),
	html.Br(),
	html.Br(),
	html.Br(),
	html.Br(),
	html.Br(),
	html.Br(),
	html.Br(),
	html.Br(),
	html.Br(),
	html.Br(),
	html.Br(),

        html.P('Show angles',
                 style={
                        'font-family': 'Dosis',
                        'display': 'inline',
                        'font-size': '3rem',
                        'color': '#4D637F'}),

	html.Br(),
	html.P('All flight distances are fixed to 1.',
                         style={
                        'font-family': 'Dosis',
                        'display': 'inline',
                        'font-size': '1.5rem',
                        'color': '#4D637F'}),



	dcc.Dropdown(
                    options=[{'label': 'Show theta* (D* frame)', 'value': 'Dst'},
                                        {'label': 'Show theta_l (W frame)', 'value': 'W'},
                                        {'label': 'Show chi (B frame)', 'value': 'B'},
                                        {'label': 'Show theta* and theta_l', 'value': 'two'}],

                        value='Dst',
                        id='dropdown-angle_Dst'
                ),


  	dcc.Graph(id='show-angle_Dst'),

	], className='six columns', style={'margin':0}), 
       
	
	html.Div([
		dcc.RadioItems(
	            options=[
			{'label':'3D one event','value':'3D'},
			{'label':'XY projection','value':'XY'},
			{'label':'YZ projection','value':'YZ'},
			{'label':'ZX projection','value':'ZX'}
		    ],
		    value='3D',
		    id='which-D_Dst'
		),	
		html.Br(),
		html.P('Select frame:',
		 style={
                    'font-family': 'Dosis',
                    'display': 'inline',
                    'font-size': '3rem',
                    'color': '#4D637F'
		}),

		dcc.Dropdown(
		    options=[{'label': 'Lab', 'value': 'lab'},
					{'label': 'COM of B', 'value': 'B'},
					{'label': 'COM of D*', 'value': 'Dst'},
			     		{'label': 'COM of W', 'value': 'W'}],
				
			value='lab',
			id='dropdown-frame_Dst'
		),

		dcc.Graph(id='frame-graph_Dst',style={'height':500,'width':500}),
		html.Br(),
		html.P('The line color in the 3D plot stands for the magnitude of the momentum in the LAB frame, darker color means larger momentum.',
			 style={
                    	'font-family': 'Dosis',
                    	'display': 'inline',
                    	'font-size': '1.5rem',
                    	'color': '#4D637F'}),
                html.Br(),
                html.P('Histograms',
                         style={
                        'font-family': 'Dosis',
                        'display': 'inline',
                        'font-size': '3rem',
                        'color': '#4D637F'}) 
		], className='six columns', style={'margin':0}),
	html.Br(),
	html.Div([
		html.Div([

			dcc.Dropdown(
                               options=[{'label': 'chi', 'value': 'chid'},
                                        {'label': 'costheta*', 'value': 'std'},
					{'label': 'costhetal', 'value': 'ld'},
                                        {'label': 'q2', 'value': 'q2'}],


                        		value='chid',	
                        		id='dropdown-dist_Dst'
                			),
			dcc.Graph(id='histog_Dst'),
		], className='six columns', style={'margin':0})]),
    html.Br(),
    dcc.Link('Navigate to "D2*(2460)0"', href='/apps/D2460'),
    html.Br(),
    dcc.Link('Navigate to "D1(2420)0" ', href='/apps/D2420'),
    html.Br(),
    dcc.Link('Navigate to "D0"', href='/apps/D0'),


])
		
#--------------------------------------------------------------
#Callback for the showangle graph, clickData and Dropdown Menu |
#--------------------------------------------------------------

@app.callback(
Output('show-angle_Dst', 'figure'),
[Input('phase-space_Dst', 'clickData'),
Input('dropdown-angle_Dst','value')])


def drawangle_Dst(selection,choice):
	if selection is None:
		return {}
	else:
		i=selection['points'][0]['pointNumber']

		if choice=='two':
			angleBx,angleBy,angleBz=(0,0,0)
			angleDstx,angleDsty,angleDstz=(nouvDst.z[i]/nouvDst.p[i],nouvDst.x[i]/nouvDst.p[i],nouvDst.y[i]/nouvDst.p[i])
			angleDx,angleDy,angleDz=(nouvD0.z[i]/nouvD0.p[i]+angleDstx,nouvD0.x[i]/nouvD0.p[i]+angleDsty,nouvD0.y[i]/nouvD0.p[i]+angleDstz)
			anglepix,anglepiy,anglepiz=(-nouvD0.z[i]/nouvD0.p[i]+angleDstx,-nouvD0.x[i]/nouvD0.p[i]+angleDsty,-nouvD0.y[i]/nouvD0.p[i]+angleDstz)
			angleQx,angleQy,angleQz=(-angleDstx,-angleDsty,-angleDstz)
			angletaux,angletauy,angletauz=(angleQx+nouvtau.z[i]/nouvtau.p[i],angleQy+nouvtau.x[i]/nouvtau.p[i],angleQz+nouvtau.y[i]/nouvtau.p[i])             
			anglenux,anglenuy,anglenuz=(angleQx-nouvtau.z[i]/nouvtau.p[i],angleQy-nouvtau.x[i]/nouvtau.p[i],angleQz-nouvtau.y[i]/nouvtau.p[i])             
			d1x,d1y,d1z=(2*angleQx,2*angleQy,2*angleQz)
			d2x,d2y,d2z=(2*angleDstx,2*angleDsty,2*angleDstz)
			dash1=go.Scatter3d(x=[angleQx,d1x],y=[angleQy,d1y],z=[angleQz,d1z],mode='lines',line = dict(color = ('rgb(20, 20, 255)'),width = 3,dash='dot'))
			dash2=go.Scatter3d(x=[angleDstx,d2x],y=[angleDsty,d2y],z=[angleDstz,d2z],mode='lines',line = dict(color = ('rgb(20, 20, 255)'),width = 3,dash='dot'))


			mesh1 = go.Mesh3d(x=[d1x,angletaux,angleQx],y=[d1y,angletauy,angleQy],z=[d1z,angletauz,angleQz],                   
                   		opacity=0.4,
                   		color='#3E3A3A')

			mesh2 = go.Mesh3d(x=[d2x,angleDx,angleDstx],y=[d2y,angleDy,angleDsty],z=[d2z,angleDz,angleDstz],
                                opacity=0.4,
                                color='#3E3A3A')
			traceQ=go.Scatter3d(x=[angleBx,angleQx],y=[angleBy,angleQy],z=[angleBz,angleQz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['B', 'W'],textposition='top right',line = dict(color = ('rgb(0, 0, 255)'),width = 3))

			traceDst=go.Scatter3d(x=[angleBx,angleDstx],y=[angleBy,angleDsty],z=[angleBz,angleDstz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'D*'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))

			tracenu=go.Scatter3d(x=[angleQx,anglenux],y=[angleQy,anglenuy],z=[angleQz,anglenuz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'nu'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))

			tracetau=go.Scatter3d(x=[angleQx,angletaux],y=[angleQy,angletauy],z=[angleQz,angletauz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'tau'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))

			traceD0=go.Scatter3d(x=[angleDstx,angleDx],y=[angleDsty,angleDy],z=[angleDstz,angleDz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'D0'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))

			tracepi=go.Scatter3d(x=[angleDstx,anglepix],y=[angleDsty,anglepiy],z=[angleDstz,anglepiz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'pi'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))

			dataangle=[traceQ,traceDst,tracenu,tracetau,traceD0,tracepi,dash1,dash2,mesh1,mesh2]


			
		else: 
			COM=LorentzVector(df[choice+'_PX_TRUE'],df[choice+'_PY_TRUE'],df[choice+'_PZ_TRUE'],df[choice+'_E_TRUE'])
			liste_part=change_frame(COM)
			[newB,newDst,newtau,newD0,newnuB,newK,newpiDst,newpitau1,newpitau2,newpitau3,newnutau,newpiK]=liste_part
			newW=W.boost(-COM.boostp3)

			if choice=='Dst':			
				angleDstx,angleDsty,angleDstz=(0,0,0)
				angleBx,angleBy,angleBz=(newB.z[i]/newB.p[i],newB.x[i]/newB.p[i],newB.y[i]/newB.p[i])
	
	
				angleQx,angleQy,angleQz=(newW.z[i]/newW.p[i]+angleBx,newW.x[i]/newW.p[i]+angleBy,newW.y[i]/newW.p[i]+angleBz)

				d1x,d1y,d1z=(-angleBx,-angleBy,-angleBz)
				angleDx,angleDy,angleDz=(newD0.z[i]/newD0.p[i]+angleDstx,newD0.x[i]/newD0.p[i]+angleDsty,newD0.y[i]/newD0.p[i]+angleDstz)
				anglepix,anglepiy,anglepiz=(newpiDst.z[i]/newpiDst.p[i]+angleDstx,newpiDst.x[i]/newpiDst.p[i]+angleDsty,newpiDst.y[i]/newpiDst.p[i]+angleDstz)
				angletaux,angletauy,angletauz=(angleQx+newtau.z[i]/newtau.p[i],angleQy+newtau.x[i]/newtau.p[i],angleQz+newtau.y[i]/newtau.p[i])
				anglenux,anglenuy,anglenuz=(angleQx+newnuB.z[i]/newnuB.p[i],angleQy+newnuB.x[i]/newnuB.p[i],angleQz+newnuB.y[i]/newnuB.p[i])


				traceQ=go.Scatter3d(x=[angleBx,angleQx],y=[angleBy,angleQy],z=[angleBz,angleQz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['B', 'W'],textposition='top right',line = dict(color = ('rgb(0, 0, 255)'),width = 3))

				traceDst=go.Scatter3d(x=[angleBx,angleDstx],y=[angleBy,angleDsty],z=[angleBz,angleDstz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'D*'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))
				tracenu=go.Scatter3d(x=[angleQx,anglenux],y=[angleQy,anglenuy],z=[angleQz,anglenuz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'nu'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))

				tracetau=go.Scatter3d(x=[angleQx,angletaux],y=[angleQy,angletauy],z=[angleQz,angletauz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'tau'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))
				traceD0=go.Scatter3d(x=[angleDstx,angleDx],y=[angleDsty,angleDy],z=[angleDstz,angleDz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'D0'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))
				tracepi=go.Scatter3d(x=[angleDstx,anglepix],y=[angleDsty,anglepiy],z=[angleDstz,anglepiz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'pi'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))


				dash1=go.Scatter3d(x=[angleDstx,d1x],y=[angleDsty,d1y],z=[angleDstz,d1z],mode='lines',line = dict(color = ('rgb(20, 20, 255)'),width = 3,dash='dot'))
				mesh1 = go.Mesh3d(x=[d1x,0,angleDx],y=[d1y,0,angleDy],z=[d1z,0,angleDz],
                                opacity=0.4,
                                color='#3E3A3A')
				dataangle=[traceQ,traceDst,tracenu,tracetau,traceD0,tracepi,dash1,mesh1]


			if choice=='W':
				angleQx,angleQy,angleQz=(0,0,0)
				angleBx,angleBy,angleBz=(newB.z[i]/newB.p[i],newB.x[i]/newB.p[i],newB.y[i]/newB.p[i])
				angleDstx,angleDsty,angleDstz=(newDst.z[i]/newDst.p[i]+angleBx,newDst.x[i]/newDst.p[i]+angleBy,newDst.y[i]/newDst.p[i]+angleBz)
				angleDx,angleDy,angleDz=(newD0.z[i]/newD0.p[i]+angleDstx,newD0.x[i]/newD0.p[i]+angleDsty,newD0.y[i]/newD0.p[i]+angleDstz)
				anglepix,anglepiy,anglepiz=(newpiDst.z[i]/newpiDst.p[i]+angleDstx,newpiDst.x[i]/newpiDst.p[i]+angleDsty,newpiDst.y[i]/newpiDst.p[i]+angleDstz)
				angletaux,angletauy,angletauz=(angleQx+newtau.z[i]/newtau.p[i],angleQy+newtau.x[i]/newtau.p[i],angleQz+newtau.y[i]/newtau.p[i])
				anglenux,anglenuy,anglenuz=(angleQx+newnuB.z[i]/newnuB.p[i],angleQy+newnuB.x[i]/newnuB.p[i],angleQz+newnuB.y[i]/newnuB.p[i])
				d1x,d1y,d1z=(-angleBx,-angleBy,-angleBz)
				dash1=go.Scatter3d(x=[angleQx,d1x],y=[angleQy,d1y],z=[angleQz,d1z],mode='lines',line = dict(color = ('rgb(20, 20, 255)'),width = 3,dash='dot'))
				traceQ=go.Scatter3d(x=[angleBx,angleQx],y=[angleBy,angleQy],z=[angleBz,angleQz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['B', 'W'],textposition='top right',line = dict(color = ('rgb(0, 0, 255)'),width = 3))

				traceDst=go.Scatter3d(x=[angleBx,angleDstx],y=[angleBy,angleDsty],z=[angleBz,angleDstz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'D*'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))

				tracenu=go.Scatter3d(x=[angleQx,anglenux],y=[angleQy,anglenuy],z=[angleQz,anglenuz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'nu'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))
				tracetau=go.Scatter3d(x=[angleQx,angletaux],y=[angleQy,angletauy],z=[angleQz,angletauz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'tau'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))
				traceD0=go.Scatter3d(x=[angleDstx,angleDx],y=[angleDsty,angleDy],z=[angleDstz,angleDz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'D0'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))
				tracepi=go.Scatter3d(x=[angleDstx,anglepix],y=[angleDsty,anglepiy],z=[angleDstz,anglepiz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'pi'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))



				mesh1 = go.Mesh3d(x=[d1x,0,angletaux],y=[d1y,0,angletauy],z=[d1z,0,angletauz],
				opacity=0.4,
				color='#3E3A3A')
				dataangle=[traceQ,traceDst,tracenu,tracetau,traceD0,tracepi,dash1,mesh1]
        
			if choice=='B':
				angleBx,angleBy,angleBz=(0,0,0)
				angleDstx,angleDsty,angleDstz=(newDst.x[i]/newDst.p[i]+angleBx,newDst.y[i]/newDst.p[i]+angleBy,newDst.z[i]/newDst.p[i]+angleBz)

				angleQx,angleQy,angleQz=(newW.x[i]/newW.p[i]+angleBx,newW.y[i]/newW.p[i]+angleBy,newW.z[i]/newW.p[i]+angleBz)

				angleDx,angleDy,angleDz=(newD0.x[i]/newD0.p[i]+angleDstx,newD0.y[i]/newD0.p[i]+angleDsty,newD0.z[i]/newD0.p[i]+angleDstz)
				anglepix,anglepiy,anglepiz=(newpiDst.x[i]/newpiDst.p[i]+angleDstx,newpiDst.y[i]/newpiDst.p[i]+angleDsty,newpiDst.z[i]/newpiDst.p[i]+angleDstz)
				angletaux,angletauy,angletauz=(angleQx+newtau.x[i]/newtau.p[i],angleQy+newtau.y[i]/newtau.p[i],angleQz+newtau.z[i]/newtau.p[i])
				anglenux,anglenuy,anglenuz=(angleQx+newnuB.x[i]/newnuB.p[i],angleQy+newnuB.y[i]/newnuB.p[i],angleQz+newnuB.z[i]/newnuB.p[i])

				vec1=Vector3D(newW.x[i]/newW.p[i],newW.y[i]/newW.p[i],newW.z[i]/newW.p[i])      #unit vector W
				vec2=Vector3D(newtau.x[i]/newtau.p[i],newtau.y[i]/newtau.p[i],newtau.z[i]/newtau.p[i])	#unit vector tau
				vec3=Vector3D(newD0.x[i]/newD0.p[i],newD0.y[i]/newD0.p[i],newD0.z[i]/newD0.p[i])	#unit vector D0
				vec4=Vector3D(newDst.x[i]/newDst.p[i],newDst.y[i]/newDst.p[i],newDst.z[i]/newDst.p[i])	#unit vector D*
			
				cross1=vec1.cross((vec2.cross(vec1)).unit)
				cross1=cross1.unit
				cross2=vec4.cross((vec3.cross(vec4)).unit)
				cross2=cross2.unit
				d1x,d1y,d1z=(cross1.x,cross1.y,cross1.z)
				d2x,d2y,d2z=(cross2.x,cross2.y,cross2.z)

				traceQ=go.Scatter3d(x=[angleBx,angleQx],y=[angleBy,angleQy],z=[angleBz,angleQz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['B', 'W'],textposition='top right',line = dict(color = ('rgb(0, 0, 255)'),width = 3))
				traceDst=go.Scatter3d(x=[angleBx,angleDstx],y=[angleBy,angleDsty],z=[angleBz,angleDstz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'D*'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))
				tracenu=go.Scatter3d(x=[angleQx,anglenux],y=[angleQy,anglenuy],z=[angleQz,anglenuz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'nu'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))
				tracetau=go.Scatter3d(x=[angleQx,angletaux],y=[angleQy,angletauy],z=[angleQz,angletauz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'tau'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))
				traceD0=go.Scatter3d(x=[angleDstx,angleDx],y=[angleDsty,angleDy],z=[angleDstz,angleDz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'D0'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))
				tracepi=go.Scatter3d(x=[angleDstx,anglepix],y=[angleDsty,anglepiy],z=[angleDstz,anglepiz],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'pi'],textposition='top left',line = dict(color = ('rgb(0, 0, 255)'),width = 3))



		

				mesh1 = go.Mesh3d(x=[angleQx,angleDstx,newtau.x[i]/newtau.p[i]+newDst.x[i],newtau.x[i]/newtau.p[i]-newDst.x[i]],y=[angleQy,angleDsty,newtau.y[i]/newtau.p[i]+newDst.y[i],newtau.y[i]/newtau.p[i]-newDst.y[i]],z=[angleQz,angleDstz,newtau.z[i]/newtau.p[i]+newDst.z[i],newtau.z[i]/newtau.p[i]-newDst.z[i]],
					opacity=0.4,
					color='#00cc51')
				mesh2 = go.Mesh3d(x=[angleQx,angleDstx,newD0.x[i]/newD0.p[i]+newDst.x[i],newD0.x[i]/newD0.p[i]-newDst.x[i]],y=[angleQy,angleDsty,newD0.y[i]/newD0.p[i]+newDst.y[i],newtau.y[i]/newtau.p[i]-newDst.y[i]],z=[angleQz,angleDstz,newD0.z[i]/newD0.p[i]+newDst.z[i],newD0.z[i]/newD0.p[i]-newDst.z[i]],
                                        opacity=0.4,
                                        color='#00cc51')
				mesh3=go.Mesh3d(x=[0,d1x,d2x],y=[0,d1y,d2y],z=[0,d1z,d2z],
                                        opacity=0.4,
                                        color='#3E3A3A')

				dataangle=[traceQ,traceDst,tracenu,tracetau,traceD0,tracepi,mesh1,mesh2]

		layoutangle = go.Layout(
                                showlegend=False,
                                paper_bgcolor = '#F4F4F8',
                                plot_bgcolor = '#F4F4F8',
                                autosize=True,
                                margin=dict(t=10, b=10, l=20, r=10),
                                scene=dict(
                                camera=dict(eye=dict(x=1.5,y=1.5,z=1.5)),
                                xaxis=dict(
                              
                                title='Z',
                                titlefont=dict(
                                family='Arial, sans-serif',
                                size=18,
                                color='black'
                                ),
                                gridcolor='#bdbdbd',
                                gridwidth=2,
                                zerolinecolor='#969696',
                                zerolinewidth=4,
                                linecolor='#636363',
                                linewidth=4,
                                showbackground=True,
                                backgroundcolor='rgb(230, 230,230)'
                                ),
                                yaxis=dict(
                               
                                title='X',
                                titlefont=dict(
                                family='Arial, sans-serif',
                                size=18,
                                color='black'
                                ),
                                gridcolor='#bdbdbd',
                                gridwidth=2,
                                zerolinecolor='#969696',
                                zerolinewidth=4,
                                linecolor='#636363',
                                linewidth=4,
                                showbackground=True,
                                backgroundcolor='rgb(230, 230, 230)'
                                ),
                                zaxis=dict(
                             
                                title='Y',
                                titlefont=dict(
                                family='Arial, sans-serif',
                                size=18,
                                color='black'
                                ),
                        gridcolor='#bdbdbd',
                        gridwidth=2,
                        zerolinecolor='#969696',
                        zerolinewidth=4,
                        linecolor='#636363',
                        linewidth=4,

                        showbackground=True,
                        backgroundcolor='rgb(230, 230,230)'
                        ),
                        aspectratio = dict(x=1, y=1, z=0.7),
                        aspectmode = 'manual'
                                )
                        )



		return {'data':dataangle,'layout':layoutangle}


#--------------------------------
#                                |
#Callback for the range sliders  |
#                                |
#--------------------------------

@app.callback(
    dash.dependencies.Output('output-range-thetast_Dst', 'children'),
    [dash.dependencies.Input('choose-thetast_Dst', 'value')])

def update_output1_Dst(value):
	return 'You have selected "{}"'.format(value)
    

            
@app.callback(
    dash.dependencies.Output('output-range-thetal_Dst', 'children'),
    [dash.dependencies.Input('choose-thetal_Dst', 'value')])
def update_output2_Dst(value):
	return 'You have selected "{}"'.format(value)

            
@app.callback(
    dash.dependencies.Output('output-range-chi_Dst', 'children'),
    [dash.dependencies.Input('choose-chi_Dst', 'value')])
def update_output3_Dst(value):  
	return 'You have selected "{}"'.format(value)
                
@app.callback(
	dash.dependencies.Output('phase-space_Dst', 'figure'),
	[Input('choose-thetast_Dst', 'value'),
	Input('choose-thetal_Dst', 'value'),
	Input('choose-chi_Dst', 'value')])
def plot_phase_space_Dst(rangest,rangel,rangechi):	
	
	layout1=go.Layout(
		hovermode = 'closest',
		paper_bgcolor = '#F4F4F8',
		plot_bgcolor = '#F4F4F8',
		margin=dict(t=50, b=50, l=50, r=40),
		width=450,
		height=600,
		clickmode='event+select',
                scene=dict(camera=dict(eye=dict(x=2,y=2,z=2)),
			   xaxis=dict(title='costhetast',range=rangest),
			   yaxis=dict(title='costhetal',range=rangel),
			   zaxis=dict(title='chi',range=rangechi)
                                        ))

	return {'data': [trace_phase], 'layout':layout1}
#---------------------------
#                           |
#Callback for dropdown dist |
#                           |
#---------------------------

@app.callback(
Output('histog_Dst', 'figure'),
[Input('dropdown-dist_Dst','value')])
def distribution_Dst(choice):
	if choice=='chid':
		return {'data':[go.Histogram(x=chi)],'layout':{'title':'chi','paper_bgcolor' : '#F4F4F8', 'plot_bgcolor' : '#F4F4F8'}}
	if choice=='std':
		return {'data':[go.Histogram(x=costhetast)],'layout':{'title':'costheta*','paper_bgcolor' : '#F4F4F8', 'plot_bgcolor' : '#F4F4F8'}}
	if choice=='ld':
		return {'data':[go.Histogram(x=costhetal)],'layout':{'title':'costhetal','paper_bgcolor' : '#F4F4F8', 'plot_bgcolor' : '#F4F4F8'}}
	if choice=='q2':
		q=B-Dst
		q2=q.mag2
		return {'data':[go.Histogram(x=q2)],'layout':{'title':'q2','paper_bgcolor' : '#F4F4F8', 'plot_bgcolor' : '#F4F4F8'}}



#------------------------
#			 |
#Callback for clickData  |
#                        |
#------------------------



@app.callback(
    Output('click-result_Dst', 'children'),
    [Input('phase-space_Dst', 'clickData')])
def display_click_data_Dst(clickData):
    return json.dumps(clickData, indent=2)


@app.callback(					
Output('frame-graph_Dst', 'figure'),
[Input('phase-space_Dst', 'clickData'),
Input('which-D_Dst','value'),
Input('dropdown-frame_Dst','value')])
def drawevent_Dst(selection,radio,frame):
	if selection is None:
		return {}
	else: 
		i=selection['points'][0]['pointNumber']
		dis=(df['B_FD_TRUE'][i]+df['Dst_FD_TRUE'][i]+df['Tau_FD_TRUE'][i]+df['D0_FD_TRUE'][i])/4		#the average flight distance
		cmax=B.p[i]
		def colorval(vec):              #determine the color
			val=vec.p[i]
			if val==cmax:
				return 'rgb(0,0,0)'
			if val<cmax and val>cmax/4:
				return 'rgb(0,0,83)'
			if val<=cmax/4 and val>cmax/16 :
				return 'rgb(0,0,128)'
			if val<=cmax/16 and val>cmax/64:
				return 'rgb(0,0,196)'
			else :
				return 'rgb(0,0,255)'



		if frame=='lab':
					
			PV_X,PV_Y,PV_Z=(df['B_Ori_z_TRUE'][i],df['B_Ori_x_TRUE'][i],df['B_Ori_y_TRUE'][i])
			B_X,B_Y,B_Z=(df['B_End_z_TRUE'][i],df['B_End_x_TRUE'][i],df['B_End_y_TRUE'][i])    
			Dst_X,Dst_Y,Dst_Z=(df['Dst_End_z_TRUE'][i],df['Dst_End_x_TRUE'][i],df['Dst_End_y_TRUE'][i])
			tau_X,tau_Y,tau_Z=[df['Tau_End_z_TRUE'][i],df['Tau_End_x_TRUE'][i],df['Tau_End_y_TRUE'][i]]
			D0_X,D0_Y,D0_Z=[df['D0_End_z_TRUE'][i],df['D0_End_x_TRUE'][i],df['D0_End_y_TRUE'][i]]
			nu_X=df['B_nu_PZ_TRUE'][i]*dis/df['B_nu_P_TRUE'][i]+B_X
			nu_Y=df['B_nu_PX_TRUE'][i]*dis/df['B_nu_P_TRUE'][i]+B_Y
			nu_Z=df['B_nu_PY_TRUE'][i]*dis/df['B_nu_P_TRUE'][i]+B_Z
			K_X=df['D0_K_PZ_TRUE'][i]*dis/df['D0_K_P_TRUE'][i]+D0_X
			K_Y=df['D0_K_PX_TRUE'][i]*dis/df['D0_K_P_TRUE'][i]+D0_Y
			K_Z=df['D0_K_PY_TRUE'][i]*dis/df['D0_K_P_TRUE'][i]+D0_Z
			piK_X=df['D0_Pi_PZ_TRUE'][i]*dis/df['D0_Pi_P_TRUE'][i]+D0_X
			piK_Y=df['D0_Pi_PX_TRUE'][i]*dis/df['D0_Pi_P_TRUE'][i]+D0_Y
			piK_Z=df['D0_Pi_PY_TRUE'][i]*dis/df['D0_Pi_P_TRUE'][i]+D0_Z
			piDst_X=df['Dst_Pi_PZ_TRUE'][i]*dis/df['Dst_Pi_P_TRUE'][i]+Dst_X
			piDst_Y=df['Dst_Pi_PX_TRUE'][i]*dis/df['Dst_Pi_P_TRUE'][i]+Dst_Y
			piDst_Z=df['Dst_Pi_PY_TRUE'][i]*dis/df['Dst_Pi_P_TRUE'][i]+Dst_Z
		
			pitau1_X=df['Tau_Pi1_PZ_TRUE'][i]*dis/df['Tau_Pi1_P_TRUE'][i]+tau_X
			pitau1_Y=df['Tau_Pi1_PX_TRUE'][i]*dis/df['Tau_Pi1_P_TRUE'][i]+tau_Y
			pitau1_Z=df['Tau_Pi1_PY_TRUE'][i]*dis/df['Tau_Pi1_P_TRUE'][i]+tau_Z
		
			pitau2_X=df['Tau_Pi2_PZ_TRUE'][i]*dis/df['Tau_Pi2_P_TRUE'][i]+tau_X
			pitau2_Y=df['Tau_Pi2_PX_TRUE'][i]*dis/df['Tau_Pi2_P_TRUE'][i]+tau_Y
			pitau2_Z=df['Tau_Pi2_PY_TRUE'][i]*dis/df['Tau_Pi2_P_TRUE'][i]+tau_Z
			pitau3_X=df['Tau_Pi3_PZ_TRUE'][i]*dis/df['Tau_Pi3_P_TRUE'][i]+tau_X
			pitau3_Y=df['Tau_Pi3_PX_TRUE'][i]*dis/df['Tau_Pi3_P_TRUE'][i]+tau_Y
			pitau3_Z=df['Tau_Pi3_PY_TRUE'][i]*dis/df['Tau_Pi3_P_TRUE'][i]+tau_Z
			nutau_X=df['Tau_nu_PZ_TRUE'][i]*dis/df['Tau_nu_P_TRUE'][i]+tau_X
			nutau_Y=df['Tau_nu_PX_TRUE'][i]*dis/df['Tau_nu_P_TRUE'][i]+tau_Y
			nutau_Z=df['Tau_nu_PY_TRUE'][i]*dis/df['Tau_nu_P_TRUE'][i]+tau_Z





		else:


			COM=LorentzVector(df[frame+'_PX_TRUE'],df[frame+'_PY_TRUE'],df[frame+'_PZ_TRUE'],df[frame+'_E_TRUE'])
			liste_part=change_frame(COM)

			[newB,newDst,newtau,newD0,newnuB,newK,newpiDst,newpitau1,newpitau2,newpitau3,newnutau,newpiK]=liste_part
			PV_X,PV_Y,PV_Z=(0,0,0)
			if frame=='B':
				B_X,B_Y,B_Z=(0,0,0)
			else:
				B_X,B_Y,B_Z=[newB.z[i]*df['B_FD_TRUE'][i]/newB.p[i],newB.x[i]*df['B_FD_TRUE'][i]/newB.p[i],newB.y[i]*df['B_FD_TRUE'][i]/newB.p[i]]
	
			Dst_X,Dst_Y,Dst_Z=[newDst.z[i]*df['Dst_FD_TRUE'][i]/newDst.p[i]+B_X,newDst.x[i]*df['Dst_FD_TRUE'][i]/newDst.p[i]+B_Y,newDst.y[i]*df['Dst_FD_TRUE'][i]/newDst.p[i]+B_Z]
			tau_X,tau_Y,tau_Z=[newtau.z[i]*df['Tau_FD_TRUE'][i]/newtau.p[i]+B_X,newtau.x[i]*df['Tau_FD_TRUE'][i]/newtau.p[i]+B_Y,newtau.y[i]*df['Tau_FD_TRUE'][i]/newtau.p[i]+B_Z]

			D0_X,D0_Y,D0_Z=[newD0.z[i]*df['D0_FD_TRUE'][i]/newD0.p[i]+Dst_X,newD0.x[i]*df['D0_FD_TRUE'][i]/newD0.p[i]+Dst_Y,newD0.y[i]*df['D0_FD_TRUE'][i]/newD0.p[i]+Dst_Z]



			nu_X=newnuB.z[i]*dis/newnuB.p[i]+B_X
			nu_Y=newnuB.x[i]*dis/newnuB.p[i]+B_Y
			nu_Z=newnuB.y[i]*dis/newnuB.p[i]+B_Z
			K_X=newK.z[i]*dis/newK.p[i]+D0_X
			K_Y=newK.x[i]*dis/newK.p[i]+D0_Y
			K_Z=newK.y[i]*dis/newK.p[i]+D0_Z
			piK_X=newpiK.z[i]*dis/newpiK.p[i]+D0_X
			piK_Y=newpiK.x[i]*dis/newpiK.p[i]+D0_Y
			piK_Z=newpiK.y[i]*dis/newpiK.p[i]+D0_Z
			piDst_X=newpiDst.z[i]*dis/newpiDst.p[i]+Dst_X
			piDst_Y=newpiDst.x[i]*dis/newpiDst.p[i]+Dst_Y
			piDst_Z=newpiDst.y[i]*dis/newpiDst.p[i]+Dst_Z

			pitau1_X=newpitau1.z[i]*dis/newpitau1.p[i]+tau_X
			pitau1_Y=newpitau1.x[i]*dis/newpitau1.p[i]+tau_Y
			pitau1_Z=newpitau1.y[i]*dis/newpitau1.p[i]+tau_Z

			pitau2_X=newpitau2.z[i]*dis/newpitau2.p[i]+tau_X
			pitau2_Y=newpitau2.x[i]*dis/newpitau2.p[i]+tau_Y
			pitau2_Z=newpitau2.y[i]*dis/newpitau2.p[i]+tau_Z

			pitau3_X=newpitau3.z[i]*dis/newpitau3.p[i]+tau_X
			pitau3_Y=newpitau3.x[i]*dis/newpitau3.p[i]+tau_Y
			pitau3_Z=newpitau3.y[i]*dis/newpitau3.p[i]+tau_Z
                      
			nutau_X=newnutau.z[i]*dis/newnutau.p[i]+tau_X 
			nutau_Y=newnutau.x[i]*dis/newnutau.p[i]+tau_Y            
			nutau_Z=newnutau.y[i]*dis/newnutau.p[i]+tau_Z		


		xrange=max(abs([PV_X,B_X,Dst_X,D0_X,tau_X,K_X,nutau_X,nu_X,pitau1_X,
                                pitau2_X,pitau3_X,piDst_X,piK_X]))
		yrange=max(abs([PV_Y,B_Y,Dst_Y,D0_Y,tau_Y,K_Y,nutau_Y,nu_Y,pitau1_Y,
                         pitau2_Y,pitau3_Y,piDst_Y,piK_Y]))
		zrange=max(abs([PV_Z,B_Z,Dst_Z,D0_Z,tau_Z,K_Z,nutau_Z,nu_Z,pitau1_Z,
                                pitau2_Z,pitau3_Z,piDst_Z,piK_Z]))

		if radio=='3D':
			[traceB,traceDst,tracetau,traceD0,tracenuB,traceK]=[go.Scatter3d(x=[PV_X,B_X],y=[PV_Y,B_Y],z=[PV_Z,B_Z],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['PV', ''],textposition='top left',line = dict(color =colorval(B) ,width = 3)),go.Scatter3d(x=[B_X,Dst_X],y=[B_Y,Dst_Y],z=[B_Z,Dst_Z],mode='lines+markers+text',marker=dict(size=5,color="rgb(5,200,5)",opacity=0.8),text=['', 'D*'],textposition='top left',line = dict(width=3,color=colorval(Dst))),go.Scatter3d(x=[B_X,tau_X],y=[B_Y,tau_Y],z=[B_Z,tau_Z],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'tau'],textposition='top left',line = dict(color = colorval(tau),width = 3)),go.Scatter3d(x=[Dst_X,D0_X],y=[Dst_Y,D0_Y],z=[Dst_Z,D0_Z],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'D0'],textposition='top left',line = dict(color = colorval(D0),width = 3)),go.Scatter3d(x=[B_X,nu_X],y=[B_Y,nu_Y],z=[B_Z,nu_Z],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)",opacity=0.8),text=['', 'nuB'],textposition='top left',line = dict(color = colorval(nuB),width = 3)),go.Scatter3d(x=[D0_X,K_X],y=[D0_Y,K_Y],z=[D0_Z,K_Z],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'K'],textposition='top left',line = dict(color =colorval(K),width = 3))]
			tracepiK=go.Scatter3d(x=[D0_X,piK_X],y=[D0_Y,piK_Y],z=[D0_Z,piK_Z],mode='lines+markers+text',marker=dict(size=5,color="rgb(5,200,5)", opacity=0.8),text=['', 'pi'],textposition='top left',line = dict(color=colorval(piK),width = 3))
    		
			tracepiDst=go.Scatter3d(x=[Dst_X,piDst_X],y=[Dst_Y,piDst_Y],z=[Dst_Z,piDst_Z],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'pi'],textposition='top left',line = dict(color = colorval(piDst),width = 3))

			tracepitau1=go.Scatter3d(x=[tau_X,pitau1_X],y=[tau_Y,pitau1_Y],z=[tau_Z,pitau1_Z],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'pi'],textposition='top left',line = dict(color = colorval(pitau1),width = 3))

			tracepitau2=go.Scatter3d(x=[tau_X,pitau2_X],y=[tau_Y,pitau2_Y],z=[tau_Z,pitau2_Z],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'pi'],textposition='top left',line = dict(color = colorval(pitau2),width = 3))


			tracepitau3=go.Scatter3d(x=[tau_X,pitau3_X],y=[tau_Y,pitau3_Y],z=[tau_Z,pitau3_Z],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'pi'],textposition='top left',line = dict(color = colorval(pitau3),width = 3))


			tracenutau=go.Scatter3d(x=[tau_X,nutau_X],y=[tau_Y,nutau_Y],z=[tau_Z,nutau_Z],mode='lines+markers+text',marker=dict(size=5,color= "rgb(5,200,5)", opacity=0.8),text=['', 'nu'],textposition='top left',line = dict(color = colorval(nutau),width = 3))	


 
			layout_event = go.Layout(
				showlegend=False,	
				#width=400,
				#height=400,
				paper_bgcolor = '#F4F4F8',
				plot_bgcolor = '#F4F4F8',
				autosize=True,
				margin=dict(t=10, b=10, l=20, r=10),	
				scene=dict(
				camera=dict(eye=dict(x=1.5,y=1.5,z=1.5)),
				xaxis=dict(
				range=[-abs(xrange), abs(xrange)],
				title='Z direction [nm]',
				titlefont=dict(
				family='Arial, sans-serif',
				size=18,
				color='black'
				),
				gridcolor='#bdbdbd',
				gridwidth=2,
				zerolinecolor='#969696',
				zerolinewidth=4,
				linecolor='#636363',
				linewidth=4,
				showbackground=True,
				backgroundcolor='rgb(230, 230,230)'
				),
				yaxis=dict(
				range=[-abs(yrange), abs(yrange)],
				title='X direction [nm]',
				titlefont=dict(
				family='Arial, sans-serif',
				size=18,
				color='black'
				),
				gridcolor='#bdbdbd',
				gridwidth=2,
				zerolinecolor='#969696',
				zerolinewidth=4,
				linecolor='#636363',
				linewidth=4,
				showbackground=True,
				backgroundcolor='rgb(230, 230, 230)'
				),
				zaxis=dict(
				range=[-abs(zrange), abs(zrange)],
				title='Y direction [nm]',
				titlefont=dict(
				family='Arial, sans-serif',
				size=18,
				color='black'
				),
			gridcolor='#bdbdbd',
                        gridwidth=2,
			zerolinecolor='#969696',
			zerolinewidth=4,
			linecolor='#636363',
			linewidth=4,

			showbackground=True,
			backgroundcolor='rgb(230, 230,230)'
			),
			aspectratio = dict(x=1, y=1, z=0.7),
			aspectmode = 'manual'
				)
			)
                	
		elif radio=='ZX':

			traceB=go.Scatter(x=[PV_X,B_X],y=[PV_Y,B_Y],mode='lines+markers+text',text=['PV','B'],textposition='top center',line=dict(color='darkblue',width=2))


			traceDst=go.Scatter(x=[B_X,Dst_X],y=[B_Y,Dst_Y],mode='lines+markers+text',text=['','D*'],textposition='top center',line=dict(color='darkblue',width=2))


			tracetau=go.Scatter(x=[B_X,tau_X],y=[B_Y,tau_Y],mode='lines+markers+text',text=['','tau'],textposition='top center',line=dict(color='darkblue',width=2))

			tracenuB=go.Scatter(x=[B_X,nu_X],y=[B_Y,nu_Y],mode='lines+markers+text',text=['','nu'],textposition='top center',line=dict(color='darkblue',width=2))

			traceD0=go.Scatter(x=[Dst_X,D0_X],y=[Dst_Y,D0_Y],mode='lines+markers+text',text=['','D0'],textposition='top center',line=dict(color='darkblue',width=2))

			traceK=go.Scatter(x=[D0_X,K_X],y=[D0_Y,K_Y],mode='lines+markers+text',text=['','K'],textposition='top center',line=dict(color='darkblue',width=2))

			tracepiK=go.Scatter(x=[D0_X,piK_X],y=[D0_Y,piK_Y],mode='lines+markers+text',text=['', 'pi'],textposition='top center',line=dict(color='darkblue',width=2))

			tracepiDst=go.Scatter(x=[Dst_X,piDst_X],y=[Dst_Y,piDst_Y],mode='lines+markers+text',text=['', 'pi'],textposition='top center',line=dict(color='darkblue',width=2))

			tracepitau1=go.Scatter(x=[tau_X,pitau1_X],y=[tau_Y,pitau1_Y],mode='lines+markers+text',text=['', 'pi'],textposition='top center',line=dict(color='darkblue',width=2))

			tracepitau2=go.Scatter(x=[tau_X,pitau2_X],y=[tau_Y,pitau2_Y],mode='lines+markers+text',text=['', 'pi'],textposition='top center',line=dict(color='darkblue',width=2))

			tracepitau3=go.Scatter(x=[tau_X,pitau3_X],y=[tau_Y,pitau3_Y],mode='lines+markers+text',text=['', 'pi'],textposition='top center',line=dict(color='darkblue',width=2))
			tracenutau=go.Scatter(x=[tau_X,nutau_X],y=[tau_Y,nutau_Y],mode='lines+markers+text',text=['', 'nu'],textposition='top center',line=dict(color='darkblue',width=2))



			layout_event = go.Layout(
    			showlegend=False,
			paper_bgcolor = '#F4F4F8',
			plot_bgcolor = '#F4F4F8',
			xaxis=dict(
				title='Z direction [nm]',
        			showgrid=True,
        			zeroline=True,
        			showline=True,
        			mirror='ticks',
        			gridcolor='#bdbdbd',
        			gridwidth=2,
        			zerolinecolor='#969696',
        			zerolinewidth=4,
        			linecolor='#636363',
        			linewidth=6
    			),
    			yaxis=dict(
				title='X direction [nm]',
				showgrid=True,
        			zeroline=True,
        			showline=True,
        			mirror='ticks',
        			gridcolor='#bdbdbd',
        			gridwidth=2,
        			zerolinecolor='#969696',
        			zerolinewidth=4,
        			linecolor='#636363',
        			linewidth=6
    				)
			)


		elif radio=='XY':

			traceB=go.Scatter(x=[PV_Y,B_Y],y=[PV_Z,B_Z],mode='lines+markers+text',text=['PV','B'],textposition='top center',line=dict(color='darkblue',width=2))

			traceDst=go.Scatter(x=[B_Y,Dst_Y],y=[B_Z,Dst_Z],mode='lines+markers+text',text=['','D*'],textposition='top center',line=dict(color='darkblue',width=2))

			tracetau=go.Scatter(x=[B_Y,tau_Y],y=[B_Z,tau_Z],mode='lines+markers+text',text=['','tau'],textposition='top center',line=dict(color='darkblue',width=2))

			tracenuB=go.Scatter(x=[B_Y,nu_Y],y=[B_Z,nu_Z],mode='lines+markers+text',text=['','nu'],textposition='top center',line=dict(color='darkblue',width=2))

			traceD0=go.Scatter(x=[Dst_Y,D0_Y],y=[Dst_Z,D0_Z],mode='lines+markers+text',text=['','D0'],textposition='top center',line=dict(color='darkblue',width=2))

			traceK=go.Scatter(x=[D0_Y,K_Y],y=[D0_Z,K_Z],mode='lines+markers+text',text=['','K'],textposition='top center',line=dict(color='darkblue',width=2))

			tracepiK=go.Scatter(x=[D0_Y,piK_Y],y=[D0_Z,piK_Z],mode='lines+markers+text',text=['', 'pi'],textposition='top center',line=dict(color='darkblue',width=2))

			tracepiDst=go.Scatter(x=[Dst_Y,piDst_Y],y=[Dst_Z,piDst_Z],mode='lines+markers+text',text=['', 'pi'],textposition='top center',line=dict(color='darkblue',width=2))

			tracepitau1=go.Scatter(x=[tau_Y,pitau1_Y],y=[tau_Z,pitau1_Z],mode='lines+markers+text',text=['', 'pi'],textposition='top center',line=dict(color='darkblue',width=2))

			tracepitau2=go.Scatter(x=[tau_Y,pitau2_Y],y=[tau_Z,pitau2_Z],mode='lines+markers+text',text=['', 'pi'],textposition='top center',line=dict(color='darkblue',width=2))

			tracepitau3=go.Scatter(x=[tau_Y,pitau3_Y],y=[tau_Z,pitau3_Z],mode='lines+markers+text',text=['', 'pi'],textposition='top center',line=dict(color='darkblue',width=2))

			tracenutau=go.Scatter(x=[tau_Y,nutau_Y],y=[tau_Z,nutau_Z],mode='lines+markers+text',text=['', 'nu'],textposition='top center',line=dict(color='darkblue',width=2))


			layout_event = go.Layout(
				showlegend=False,
				paper_bgcolor = '#F4F4F8',
				plot_bgcolor = '#F4F4F8',
				xaxis=dict(
				title='X direction [nm]',
        			showgrid=True,
        			zeroline=True,
        			showline=True,
        			mirror='ticks',
        			gridcolor='#bdbdbd',
        			gridwidth=2,
        			zerolinecolor='#969696',
        			zerolinewidth=4,
        			linecolor='#636363',
        			linewidth=6
    				),
    				yaxis=dict(
				title='Y direction [nm]',
        			showgrid=True,
        			zeroline=True,
        			showline=True,
        			mirror='ticks',
        			gridcolor='#bdbdbd',
        			gridwidth=2,
        			zerolinecolor='#969696',
        			zerolinewidth=4,
        			linecolor='#636363',
        			linewidth=6
    				)
			)			

		elif radio=='YZ':

			traceB=go.Scatter(x=[PV_Z,B_Z],y=[PV_X,B_X],mode='lines+markers+text',text=['PV','B'],textposition='top center',line=dict(color='darkblue',width=2))
			traceDst=go.Scatter(x=[B_Z,Dst_Z],y=[B_X,Dst_X],mode='lines+markers+text',text=['','D*'],textposition='top center',line=dict(color='darkblue',width=2))

			tracetau=go.Scatter(x=[B_Z,tau_Z],y=[B_X,tau_X],mode='lines+markers+text',text=['','tau'],textposition='top center',line=dict(color='darkblue',width=2))

			tracenuB=go.Scatter(x=[B_Z,nu_Z],y=[B_X,nu_X],mode='lines+markers+text',text=['','nu'],textposition='top center',line=dict(color='darkblue',width=2))


			traceD0=go.Scatter(x=[Dst_Z,D0_Z],y=[Dst_X,D0_X],mode='lines+markers+text',text=['','D0'],textposition='top center',line=dict(color='darkblue',width=2))

			traceK=go.Scatter(x=[D0_Z,K_Z],y=[D0_X,K_X],mode='lines+markers+text',text=['','K'],textposition='top center',line=dict(color='darkblue',width=2))

			tracepiK=go.Scatter(x=[D0_Z,piK_Z],y=[D0_X,piK_X],mode='lines+markers+text',text=['', 'pi'],textposition='top center',line=dict(color='darkblue',width=2))
			tracepiDst=go.Scatter(x=[Dst_Z,piDst_Z],y=[Dst_X,piDst_X],mode='lines+markers+text',text=['', 'pi'],textposition='top center',line=dict(color='darkblue',width=2))

			tracepitau1=go.Scatter(x=[tau_Z,pitau1_Z],y=[tau_X,pitau1_X],mode='lines+markers+text',text=['', 'pi'],textposition='top center',line=dict(color='darkblue',width=2))
			tracepitau2=go.Scatter(x=[tau_Z,pitau2_Z],y=[tau_X,pitau2_X],mode='lines+markers+text',text=['', 'pi'],textposition='top center',line=dict(color='darkblue',width=2))

			tracepitau3=go.Scatter(x=[tau_Z,pitau3_Z],y=[tau_X,pitau3_X],mode='lines+markers+text',text=['', 'pi'],textposition='top center',line=dict(color='darkblue',width=2))

			tracenutau=go.Scatter(x=[tau_Z,nutau_Z],y=[tau_X,nutau_X],mode='lines+markers+text',text=['', 'nu'],textposition='top center',line=dict(color='darkblue',width=2))


			layout_event = go.Layout(
				showlegend=False,
				paper_bgcolor = '#F4F4F8',
				plot_bgcolor = '#F4F4F8',
    				xaxis=dict(
				title='Y direction [nm]',
        			showgrid=True,
        			zeroline=True,
        			showline=True,
        			mirror='ticks',
        			gridcolor='#bdbdbd',
       	 			gridwidth=2,
        			zerolinecolor='#969696',
        			zerolinewidth=4,
        			linecolor='#636363',
        			linewidth=6
   				),
    				yaxis=dict(
				title='Z direction [nm]',
        			showgrid=True,
        			zeroline=True,
        			showline=True,
        			mirror='ticks',
        			gridcolor='#bdbdbd',
        			gridwidth=2,
        			zerolinecolor='#969696',
        			zerolinewidth=4,
        			linecolor='#636363',
        			linewidth=6
    				)
			)



		data_event=[traceB,traceDst,tracetau,traceD0,tracenuB,traceK,tracepiK,tracepiDst,tracepitau1,tracepitau2,tracepitau3,tracenutau]
		return {'data': data_event, 'layout':layout_event    

        			}		




   
""" 
if __name__ == '__main__':
    app.run_server(debug=True)
"""
