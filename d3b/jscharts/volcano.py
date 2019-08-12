#!/usr/bin/python

# Turn on debug mode.
import sys
import cgi
import cgitb
cgitb.enable()
import csv
import numpy as np
from sklearn import linear_model
#from sklearn import decomposition
import scipy.stats as stats
#from scipy.spatial import distance
import json
import os.path
import collections
import math
#import matplotlib
#matplotlib.use('Agg')
import skbio.diversity as skdiv
import skbio.stats as skstats
import pandas as pd
#import statsmodels.api as sm
#from statsmodels.formula.api import ols

import d3bf


form = cgi.FieldStorage()
id = "emap"
d3bf.chdir( form.getvalue( "datapath" ) )
dgroup = form.getvalue( "dgroup", "none" )
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
level = form.getvalue( "level" )
ptype = form.getvalue( "ptype", "volcano" )
mmethod = form.getvalue( "mmethod", "anova" )
samples = form.getlist( "samples" )
resolution = form.getvalue( "resolution", "low" )

ilevel = int( level ) 
( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, ilevel )
( findex, mtags ) = d3bf.processtags_m( volumes, tags, dfilter )
( edata, site_ids, species_ids ) = d3bf.load_edata_m( data, ilevel, mn, ml, kdict, volumes, findex, kdnames )
aedata = np.array( edata, dtype=float )
aenorm = np.sum( aedata, axis=1 )
aedata /= aenorm.reshape( len(edata), 1 )

if len( set( mtags[ dgroup ] ) ) == 2:
	samples = list( set( mtags[ dgroup ] ) )

if len( samples ) != 2:
	print("2 groups should be selected")
	sys.exit( 0 )

logfc = []
scores = []
selspecies = []
meanlevels = []

numkeys = len( kdict )
gtag = mtags[ dgroup ]
numsamples = sum( [ ( 1 if gtag[ i ] in samples else 0 ) for i in range( len( aedata ) ) ] )

for sc in range( numkeys ):
	tnlists = [ [], [] ]
	for i in range( len( aedata ) ):
		itn = gtag[ i ]
		if not itn in samples:
			continue
		tnlists[ samples.index( itn ) ].append( aedata[i][ sc ] )
	fflag = False
	try:
		c_logfc = math.log( sum( tnlists[ 0 ] ) ) - math.log( sum( tnlists[ 1 ] ) )
	except:
		fflag = True
		
	try:
		c_level = math.log( ( sum( tnlists[0] ) + sum( tnlists[1] ) ) / numsamples )
	except:
		fflag = True
		
	if fflag:
		continue
	
	maxpv = 50.
	if mmethod == "anova":
		try:
			( fstat, pvalue ) = stats.f_oneway( *tnlists )
			c_score = min( -math.log( pvalue ), maxpv ) if pvalue > 0 else maxpv
		except:
			fflag = True
	elif mmethod == "ttest":
		try:
			( tstat, pvalue ) = stats.ttest_ind( tnlists[ 0 ], tnlists[ 1 ] )
			c_score = min( -math.log( pvalue ), maxpv ) if pvalue > 0 else maxpv
		except:
			fflag = True
	elif shist == "wilcoxon":
		try:
			( stat, pvalue ) = stats.mannwhitneyu( tnlists[ 0 ], y=tnlists[ 1 ], use_continuity=False, alternative='two-sided' )
			c_score = min( -math.log( pvalue ), maxpv ) if pvalue > 0 else maxpv
		except:
			fflag = True
	
	if not fflag:
		logfc.append( c_logfc )
		scores.append( c_score )
		selspecies.append( species_ids[ sc ] )
		meanlevels.append( c_level )
		
outlen = len( logfc )
if ptype == "volcano":
	outorder = sorted( list(range( outlen)), key = lambda k: scores[k] )
else:
	outorder = sorted( list(range( outlen)), key = lambda k: meanlevels[k] )

if resolution == "high":
	print("<svg width=\"1800\" height=\"1800\" id=\"normal\"></svg>")
else:
	print("<svg width=\"800\" height=\"800\" id=\"normal\"></svg>")
print("<script type=\"text/javascript\">")
			   


if ptype == "volcano":
	print("var ydata = %s;" % json.dumps( [ scores[ outorder[ k ] ] for k in range( outlen ) ] ))
else:
	zerolevel = min( meanlevels )
	print("var ydata = %s;" % json.dumps( [ meanlevels[ outorder[ k ] ] - zerolevel for k in range( outlen ) ] ))
print("var xdata = %s;" % json.dumps( [ logfc[ outorder[ k ] ] for k in range( outlen ) ] ))
print("var species = %s;" % json.dumps( [ selspecies[ outorder[ k ] ] for k in range( outlen ) ] ))

print("""
var margin = {top: 100, right: 100, bottom: 100, left: 100};
var svg0 = d3.select( "#normal" );
var diameter = +svg0.attr("width") - margin.top - margin.bottom;
var chart = svg0.append("g")
	.attr("transform", "translate(" + margin.left + "," + margin.top + ")");
var fs = diameter / 80;
var width = +svg0.attr( "width" ) - margin.top - margin.bottom;;
var height = +svg0.attr( "height" ) - margin.left - margin.right;

//var xaxtitles = { "linear" : "Rank", "logarithmic" : "Log( Rank )", "sqr-log" : "Sqr( Log( Rank ) )" };
//var yaxtitles = { "linear" : "Count", "logarithmic" : "Log( Count )" };
    
var xscale = d3.scaleLinear().range( [ 0, width - fs * 8 ] );
var yscale = d3.scaleLinear().range( [ height - fs * 8, 0 ] );
var xmin = d3.min( xdata );
var xmax = d3.max( xdata );
var ymin = d3.min( ydata );
var ymax = d3.max( ydata );
var xmargin = 0.1 * ( xmax - xmin );
var ymargin = 0.1 * ( ymax - ymin );

xscale.domain( [ xmin, xmax ] );	
yscale.domain( [ 0, ymax ] );	

var color = d3.scaleOrdinal(d3.schemeCategory10);

var symbol = d3.symbol();
var csymbol = symbol.type( d3.symbolCircle ).size( 0.5 * fs * fs );

chart.selectAll("dot")
	.data( xdata )
	.enter().append("circle")
	.attr("r", 2.)
	//.attr("transform", function(d,i) { return "translate(" + xscale( d ) + "," + yscale( ydata[i] ) + ")"; })
	.attr("cx", function( d ) { return xscale( d ); } )
	.attr("cy", function( d, i ) { return yscale( ydata[i] ); } )
	.style("fill", "black" ); //function( d, i ) { return ( ydata[i] < 2 ) ? "black" : ( ( d > 0 ) ? color[ "left" ] : color[ "right" ] ); } )
	;

var xax = chart.append("g").attr("transform", "translate(0," + ( height - fs * 8 ) + ")").call(d3.axisBottom(xscale));
xax.selectAll("text").style("text-anchor", "middle").style("font", 1.5 * fs + "px sans-serif" );
xax.selectAll("line").style("stroke-width", "3px");
xax.selectAll("path").style("stroke-width", "3px");

var yax = chart.append("g").call(d3.axisLeft(yscale));
yax.selectAll("text").style("text-anchor", "left").style("font", 1.5 * fs + "px sans-serif" );
yax.selectAll("line").style("stroke-width", "3px");
yax.selectAll("path").style("stroke-width", "3px");

var tfont = fs * 1.4 + "px sans-serif";

chart.append("text")             
	.attr("transform",
		"translate(" + (width/2) + " ," + 
						(height - 3 * fs ) + ")")
	.style("text-anchor", "middle")
	.style("font", 2 * fs + "px sans-serif" )
	.text( "logFC" );

chart.append("text")
	.attr("transform", "rotate(-90)")
	.attr("y", 0 - margin.left + fs * 5 )
	.attr("x",0 - (height / 2))
	.attr("dy", "1em")
	.style("text-anchor", "middle")
	.style("font", 2 * fs + "px sans-serif" )
	.text( "p-value" );      

""")


print("</script>")

