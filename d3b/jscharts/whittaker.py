#!/usr/bin/python

# Turn on debug mode.
import sys
tfp = '/home/sferanchuk/.local/lib/python2.7/site-packages'
if tfp in sys.path:
	sys.path.remove( tfp )
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
import matplotlib
matplotlib.use('Agg')
import skbio.diversity as skdiv
import pandas as pd
#import statsmodels.api as sm
#from statsmodels.formula.api import ols

import d3bf


form = cgi.FieldStorage()
id = "emap"
dgroup = form.getvalue( "dgroup", "none" )
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
level = form.getvalue( "level" )
xscale = form.getvalue( "xscale", "linear" )
yscale = form.getvalue( "yscale", "linear" )
curvetype = form.getvalue( "curvetype", "direct" )
dmarks = form.getvalue( "dmarks", "no" )
regression = form.getvalue( "regression", "no" )
resolution = form.getvalue( "resolution", "low" )

( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
( findex, gtags ) = d3bf.processtags( volumes, tags, dfilter, dgroup )
ilevel = int( level ) 
( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, ilevel )
edata = d3bf.load_edata( data, ilevel, ml, kdict, findex, gtags )

dlengths = []
gscnt = 0
for gkey in gtags:
	sedata = sorted( edata[ gscnt ], reverse=True )
	if 0 in sedata:
		dlengths.append( sedata.index( 0 ) )
	else:
		dlengths.append( len( sedata ) )
	gscnt += 1
		
		
if xscale == "logarithmic":
	dmaxx = [ math.log( x ) for x in dlengths ]
elif xscale == "sqr-log":
	dmaxx = [ math.log( x ) * math.log( x ) for x in dlengths ]
else:
	dmaxx = dlengths
	
maxx = max( dmaxx )
	
print maxx

def calc_regression( cedata ):
	sitot = sorted( cedata, reverse=True )
	if 0 in sitot:
		zind = sitot.index( 0 )
		distr = sitot[ 0 : zind ]
	else:
		distr = sitot
	datax = []
	datay = []
	for i in range( len( distr ) ):
		vy = distr[i]
		if yscale == "logarithmic":
			vy = math.log( distr[i] )
		vx = i + 1
		if xscale == "logarithmic":
			vx = math.log( i + 1 )
		elif xscale == "sqr-log":
			vx = math.log( i + 1 ) * math.log( i + 1 )
		datax.append( vx )
		datay.append( vy )
	nmin = 3
	nmax = len( distr ) / 3
	pvbest = 1
	nbest = nmax
	abest = 0
	bbest = 0
	for cnsegm in range( nmin, nmax ):
		(a_s,b_s,r,tt,stderr) = stats.linregress( datax[:cnsegm], datay[:cnsegm] )
		if tt < pvbest:
			pvbest = tt
			nbest = cnsegm
			abest = a_s
			bbest = b_s
#	print "%d %g %s" % ( nbest, abest, bbest )
	ry = []
	for x in range( nbest ):
		ry.append( datax[ x ] * abest + bbest )
	return ( [ datax[ -1 ] - datax[ 0 ], datax[ -1 ] - datax[ nbest - 1 ] ], [ ry[ 0 ], ry[ -1 ] ] )

def calc_jakovenko( cedata, cmaxx ):
	global dbest
	sitot = sorted( cedata, reverse=True )
	if 0 in sitot:
		zind = sitot.index( 0 )
		distr = sitot[ 0 : zind ]
	else:
		distr = sitot
	datax = []
	datay = []
	xn = float( len( distr ) )
	yn = float( distr[0] )
	for i in range( len( distr ) ):
		vy = math.log( distr[i] / yn )
		vx = math.log( ( i + 1 ) / xn )
		datax.append( vx )
		datay.append( vy )
	nmin = 3
	nmax = len( distr ) / 3
	pvbest = 1
	nbest = nmax
	for cnsegm in range( nmin, nmax ):
		(a_s,b_s,r,tt,stderr) = stats.linregress( datax[:cnsegm], datay[:cnsegm] )
		if tt < pvbest:
			pvbest = tt
			nbest = cnsegm
	npmin = nbest
	npmax = len( distr )
	npbest = npmin
	npcmin = 0
	npcmax = 2 * len( distr )
	npcbest = 0
	npvbest = 1
	abest = 0
	bbest = 0
	npc = 0
	for cpbest in range( nbest + 1, npmax - 2, ( npmax - nbest ) / 8 ):
		for npc in range( npcmin, npcmax, ( npcmax - npcmin ) / 8 ):
			cdx = []
			for cx in range( cpbest, npmax ):
				cdx.append( math.log( npmax + npc ) - math.log( npc + npmax - cx ) )
			cdy = distr[ cpbest :  ]
			(a_s,b_s,r,tt,stderr) = stats.linregress( cdx, cdy )
			#model = sm.OLS( cdy, cdx )
			#rr = model.fit()
			#tt = rr.f_pvalue
			#a_s = rr.params[0]
			#b_s = 0
			if tt < npvbest:
				npvbest = tt
				npbest = cpbest
				npcbest = npc
				abest = a_s
				bbest = b_s
				dbest = [ cdx, cdy, abest, bbest, tt ]
	xres = []
	yres = []
	mxres = []
	myres = []
	for x in range( npbest, len( distr ) ):
		mx = math.log( npmax + npcbest ) - math.log( npcbest + npmax - x ) 
		my = max( 1e-10, abest * mx + bbest )
		vx = x + 1
		if xscale == "logarithmic":
			vx = math.log( vx )
		elif xscale == "sqr-log":
			lx = math.log( vx )
			vx = lx * lx
		vy = my
		if yscale == "logarithmic":
			vy = math.log( my )
		xres.append( cmaxx - vx )
		yres.append( max( vy, 0 ) )
		vdy = distr[x]
		if yscale == "logarithmic":
			vdy = math.log( distr[x] )
		mxres.append( mx )
		myres.append( [ my, vy, vdy ] )
	dbest += [ mxres, myres ]
	return ( xres, yres, npbest, npcbest, npvbest )

colflag = 1 if ( dmarks == "color" or dmarks == "both" ) else 0
shapeflag = 1 if ( ( dmarks == "shape" or dmarks == "both" ) and len( edata ) < 5 ) else 0

(a,b) = calc_regression( edata[0] )		
#curvetupe = "direct"
if resolution == "high":
	print "<svg width=\"1800\" height=\"1800\" id=\"normal\"></svg>"
else:
	print "<svg width=\"800\" height=\"800\" id=\"normal\"></svg>"
print "<script type=\"text/javascript\">"
print "var data =  [ "
gsize = len( gtags )
gbeg = 1
gscnt = 0
for gkey in gtags:
	sedata = sorted( edata[ gscnt ], reverse=True )
	for i in range( len( sedata ) ):
		v = sedata[i]
		if v > 0:
			vy = v
			if yscale == "logarithmic":
				vy = math.log( v )
			vx = i + 1
			if xscale == "logarithmic":
				vx = math.log( i + 1 )
			elif xscale == "sqr-log":
				vx = math.log( i + 1 ) * math.log( i + 1 )
			vb = ","
			if gbeg == 1:
				vb = " "
				gbeg = 0
			ckey = gkey if colflag == 1 else ""
			skey = gkey if shapeflag == 1 else ""
			print "%s[ %g, %g, \"%s\", \"%s\" ]" % ( vb, dmaxx[ gscnt ] - vx, vy, ckey, skey )
	gscnt = gscnt + 1

print "];"

npres = {}

print "var rdata = [ "
gbeg = 1
gscnt = 0
for gkey in gtags:
	rx = []
	ry = []
	r2x = []
	r2y = []
	if regression != "no":
		( rx, ry ) = calc_regression( edata[ gscnt ] )
	if regression == "2 parts":
		( r2x, r2y, nbest, npbest, npvbest  ) = calc_jakovenko( edata[ gscnt ], dmaxx[ gscnt ] )
		npres[ gkey ] = [ nbest, npbest, npvbest ]
	vb = ","
	if gbeg == 1:
		vb = " "
		gbeg = 0
	print "%s { rx: %s, ry: %s, r2x: %s, r2y: %s }" % ( vb, json.dumps( rx ), json.dumps( ry ), json.dumps( r2x ), json.dumps( r2y ) )
	gscnt = gscnt + 1

print "];"
print "var xaxtype = \"%s\";" % xscale
print "var yaxtype = \"%s\";" % yscale

print """
var margin = {top: 100, right: 100, bottom: 100, left: 100};
var svg0 = d3.select( "#normal" );
var diameter = +svg0.attr("width") - margin.top - margin.bottom;
var svg = svg0.append("g")
	.attr("transform", "translate(" + margin.left + "," + margin.top + ")");
var fs = diameter / 80;
var width = +svg0.attr( "width" ) - margin.top - margin.bottom;;
var height = +svg0.attr( "height" ) - margin.left - margin.right;

var xaxtitles = { "linear" : "Rank", "logarithmic" : "Log( Rank )", "sqr-log" : "Sqr( Log( Rank ) )" };
var yaxtitles = { "linear" : "Count", "logarithmic" : "Log( Count )" };
var xaxtitle = xaxtitles[ xaxtype ];
var yaxtitle = yaxtitles[ yaxtype ];
    
var xscale = d3.scaleLinear().range([0,width - fs * 8 ]);
var yscale = d3.scaleLinear().range([0,height - fs * 8 ]);
var xmin = d3.min(data, function(d) { return d[0]; } );
var xmax = d3.max(data, function(d) { return d[0]; } );
var ymin = d3.min(data, function(d) { return d[1]; } );
var ymax = d3.max(data, function(d) { return d[1]; } );
var xmargin = 0.1 * ( xmax - xmin );
var ymargin = 0.1 * ( ymax - ymin );

xscale.domain( [ xmax + xmargin, xmin - xmargin] );	
yscale.domain( [ ymax + ymargin, ymin - ymargin ] );	
//yscale.domain([d3.min(data, function(d) { return d[1]; } )-0.2, d3.max(data, function(d) { return d[1]; } )+0.2]);	
//yscale.domain([d3.max(data, function(d) { return d[1]; } ) * 1.1, d3.min(data, function(d) { return d[1]; } ) - 1000 ]);	
	
var cValue = function(d) { return d[2]; };
var color = d3.scaleOrdinal(d3.schemeCategory10);
	
var slabels = d3.set(data.map(function(d) { return d[3];})).values();	
var shapeScale = d3.scaleOrdinal()
			.domain(slabels)
            .range([ d3.symbolCircle, d3.symbolCross, d3.symbolDiamond, d3.symbolSquare, d3.symbolTrianle ]);

           
var valueline = d3.line()
	.x(function(d) { return xscale(d[0]); })
	.y(function(d) { return yscale(d[1]); });
		
for ( i = 0; i < rdata.length; i++ )
{
	var rline = "M" + xscale( rdata[i].rx[0] ) + " " + yscale( rdata[i].ry[0] ) + " L" + xscale( rdata[i].rx[1] ) + " " + yscale( rdata[i].ry[1] );	
	svg.append( "path" ).attr("d", rline).style("stroke-width", "2px").style("stroke","black" );
	if ( rdata[i].r2x.length > 1 )
	{
		//var r2line = "M " + xscale( rdata[i].r2x[0] ) + " " + yscale( rdata[i].r2y[0] );
		var k;
		var stepvalue = ( Math.max( 1, rdata[i].r2x.length / 30 ) ).toFixed( 0 ) * 2;
		var sv2 = stepvalue / 2;
		//console.log( stepvalue )
		//console.log( sv2 )
		for ( k = sv2; k < rdata[i].r2x.length; k+= stepvalue  )
		{
			var r2step = "M " + xscale( rdata[i].r2x[ k - sv2 ] ) + " " + yscale( rdata[i].r2y[ k - sv2 ] ) + " L " + xscale( rdata[i].r2x[k] ) + " " + yscale( rdata[i].r2y[k] );
			//r2line += " L " + xscale( rdata[i].r2x[k] ) + " " + yscale( rdata[i].r2y[k] );
			svg.append( "path" ).attr("d", r2step).style("stroke-width", fs * 0.2 + "px" ).style("stroke","black" ).style( "fill", "none" );
			//console.log( r2step );
		}
		//console.log( r2line )
		//svg.append( "path" ).attr("d", r2line).style("stroke-width", fs * 0.2 + "px" ).style("stroke","black" ).style( "fill", "none" );
	}
}
	
var symbol = d3.symbol();

svg.selectAll("dot")
	.data(data)
	.enter().append("path")
	.attr("r", 5.)
	.attr("transform", function(d) { return "translate(" + xscale( d[0] ) + "," + yscale( d[1] ) + ")"; })
	.attr("d", symbol.type( function(d){ 
	return shapeScale( d[3] );
		} ).size( 0.5 * fs * fs ) ) 
	.attr("cx", function( d ) { return xscale( d[0] ); } )
	.attr("cy", function( d ) { return yscale( d[1] ); } )
	.style("fill", function(d) { return color( cValue( d ) ); } )
	;
	
k = 0;
var dline = "";

for ( i = 0; i <= data.length; i++ )
{
	if ( i == 0 || i == data.length || data[i][0] > data[ i - 1 ][0] )
	{
		if ( i != 0 ) svg.append( "path" ).attr("d", dline ).style("stroke-width", fs * 0.1 + "px" ).style("stroke", color( data[ i - 1 ][2] ) ).style( "fill", "none" );
		if ( i != data.length ) dline = "M " + xscale( data[i][0] ) + " " + yscale( data[i][1] );
	}
	else
	{
		dline += " L " + xscale( data[i][0] ) + " " + yscale( data[i][1] );
	}
}

var xax = svg.append("g").attr("transform", "translate(0," + ( height - fs * 8 ) + ")").call(d3.axisBottom(xscale));
xax.selectAll("text").style("text-anchor", "middle").style("font", 1.5 * fs + "px sans-serif" );
xax.selectAll("line").style("stroke-width", "3px");
xax.selectAll("path").style("stroke-width", "3px");

var yax = svg.append("g").call(d3.axisLeft(yscale));
yax.selectAll("text").style("text-anchor", "left").style("font", 1.5 * fs + "px sans-serif" );
yax.selectAll("line").style("stroke-width", "3px");
yax.selectAll("path").style("stroke-width", "3px");

svg.append("text")             
	.attr("transform",
		"translate(" + (width/2) + " ," + 
						(height - 3 * fs ) + ")")
	.style("text-anchor", "middle")
	.style("font", 2 * fs + "px sans-serif" )
	.text( xaxtitle );

svg.append("text")
	.attr("transform", "rotate(-90)")
	.attr("y", 0 - margin.left + fs * 5 )
	.attr("x",0 - (height / 2))
	.attr("dy", "1em")
	.style("text-anchor", "middle")
	.style("font", 2 * fs + "px sans-serif" )
	.text( yaxtitle );      



"""


print "</script>"

if len( npres ) > 0:
	print "<br><br>"
	for gkey in npres:
		print "%s %d %d %g<br>" % ( gkey, npres[gkey][0], npres[gkey][1], npres[gkey][2] )
	print "<br>"
#	print dbest
	
