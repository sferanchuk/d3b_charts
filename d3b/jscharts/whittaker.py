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
#xscale = form.getvalue( "xscale", "linear" )
#yscale = form.getvalue( "yscale", "linear" )
curvetype = form.getvalue( "curvetype", "direct" )
dmarks = form.getvalue( "dmarks", "no" )
regression = form.getvalue( "regression", "no" )
resolution = form.getvalue( "resolution", "low" )
color = form.getvalue( "color", "none" )
shape = form.getvalue( "shape", "none" )
ptype = form.getvalue( "ptype", "none" )

( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
ilevel = int( level ) 
( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, ilevel )
( findex, gtags ) = d3bf.processtags( volumes, tags, dfilter, dgroup )
if dgroup != "none":
	edata = d3bf.load_edata( data, ilevel, ml, kdict, findex, gtags )
else:
	( findex, mtags ) = d3bf.processtags_m( volumes, tags, dfilter )
	( edata, site_ids, species_ids ) = d3bf.load_edata_m( data, ilevel, mn, ml, kdict, volumes, findex, kdnames )


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

def calc_rarefaction( si ): 
	n_indiv = sum( si )
	n_otu = len( si )
	#def rcount( sn, n, x, i ):
		#return sn -  np.sum( comb( n-x, i, exact = False ) ) / comb( n, i, exact = False ) 
		#print >>sys.stderr, ( x, n, i, sn )
	#	return sn -  sum( [ comb( n-xk, i, exact = False ) for xk in x ] )  / comb( n, i, exact = False ) 
	def subsample( si, i ):
		ssi = skstats.subsample_counts( si, i )
		return np.count_nonzero( ssi )

	#def errfn(p, n, y):
	#	return ( ( ( p[0] * n / (p[1] + n ) ) - y ) ** 2 ).sum()
	#	#return ( ( p[0] * ( 1. - np.exp( n / p[1] ) ) - y ) ** 2 ).sum()
	
	i_step = max( n_indiv / 200, 1 )
	num_repeats = max( 2000 / i_step, 1 ) 
	print >>sys.stderr, ( i_step, num_repeats )
	S_max_guess = n_otu
	B_guess = int( round( n_otu / 2 ) )
	params_guess = ( S_max_guess, B_guess )
	xvals = np.arange( 1, n_indiv, i_step )
	ymtx = np.empty( ( num_repeats, len( xvals ) ), dtype=int )
	for i in range( num_repeats ):
		ymtx[i] = np.asarray( [ subsample( si, n ) for n in xvals ], dtype=int )
	yvals = ymtx.mean(0)
	return ( xvals.tolist(), yvals.tolist() )

colflag = 1 if ( dmarks == "color" or dmarks == "both" ) else 0
shapeflag = 1 if ( ( dmarks == "shape" or dmarks == "both" ) and len( edata ) <= 5 ) else 0

#(a,b) = calc_regression( edata[0] )		
#curvetupe = "direct"
if resolution == "high":
	print "<svg width=\"1800\" height=\"1800\" id=\"normal\"></svg>"
else:
	print "<svg width=\"800\" height=\"800\" id=\"normal\"></svg>"
print "<script type=\"text/javascript\">"

ldata = []
rdata = []

xscale = "---"
yscale = "---"

gtvalues = sorted( gtags.values() )
for gscnt in range( len( gtvalues ) ):
	gnum = gtvalues[ gscnt ]
	gtag = gtags.keys()[ gtags.values().index( gnum ) ]
	setot = sorted( edata[ gscnt ], reverse=True )
	if 0 in setot:
		zind = setot.index( 0 )
		sedata = setot[ 0 : zind ]
	else:
		sedata = setot
		
	if dgroup != "none":
		ckey = gtag if colflag == 1 else ""
		skey = gtag if shapeflag == 1 else ""
	else:
		gind = mtags[ "name" ].index( gtag )
		ckey = mtags[ color ][ gind ] if color != "none" else ""
		skey = mtags[ shape ][ gind ] if shape != "none" else ""
	
	if ptype == "rarefaction":
		( xdata, ydata ) = calc_rarefaction( sedata )
		xscale = "count of reads"
		yscale = "count of phylotypes"
	elif ptype == "lorentz":
		xdata = ( np.arange( len( sedata ), dtype = 'float' ) / len( sedata ) ).tolist()
		ydata = []
		ysum = 0
		ynorm = float( sum( sedata ) )
		for k in range( len( sedata ) ):
			ysum += sedata[ -( k + 1 ) ]
			ydata.append( ysum / ynorm ) 
		xscale = "relative rank"
		yscale = "relative cumulative abundance"
	else:
		ydata = np.log( sedata ).tolist()
		xrdata = range( 1, len( ydata ) + 1 )
		if ptype == "log-log":
			xdata = np.log( xrdata ).tolist()
			xscale = "log(rank)"
		elif xscale == "log-sqrlog":
			xdata = np.sqr( np.log( xrdata ) ).tolist()
			xscale = "sqr(log(rank))"
		else: 
			xdata = xrdata
			xscake = "rank"
		( rx, ry ) = ( [], [] )
		( r2x, r2y ) = ( [], [] )
		if regression != "no":
			( rx, ry ) = calc_regression( edata[ gscnt ] )
			rdata.append( [ rx, ry, r2x, r2y ] )
		yscale = "log(count)"

	ldata.append( [ xdata[:], ydata[:], ckey, skey ]  )
	
			   

print "var ldata = %s;" % json.dumps( ldata )
print "var rdata = %s;" % json.dumps( rdata )

print "var xaxtitle = \"%s\";" % xscale
print "var yaxtitle = \"%s\";" % yscale

print """
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
var xmin = d3.min( ldata, function(d) { return d3.min( d[0] ); } );
var xmax = d3.max( ldata, function(d) { return d3.max( d[0] ); } );
var ymin = d3.min( ldata, function(d) { return d3.min( d[1] ); } );
var ymax = d3.max( ldata, function(d) { return d3.max( d[1] ); } );
var xmargin = 0.1 * ( xmax - xmin );
var ymargin = 0.1 * ( ymax - ymin );

xscale.domain( [ xmin - xmargin, xmax + xmargin] );	
yscale.domain( [ ymin - ymargin, ymax + ymargin ] );	
//yscale.domain([d3.min(ldata, function(d) { return d[1]; } )-0.2, d3.max(data, function(d) { return d[1]; } )+0.2]);	
//yscale.domain([d3.max(ldata, function(d) { return d[1]; } ) * 1.1, d3.min(data, function(d) { return d[1]; } ) - 1000 ]);	
	
//var cValue = function(d) { return d[2]; };
var color = d3.scaleOrdinal(d3.schemeCategory10);
var clabels = d3.set( ldata.map(function(d) { return d[2];})).values();
var slabels = d3.set( ldata.map(function(d) { return d[3];})).values();	
var shapeScale = d3.scaleOrdinal()
			.domain(slabels)
            .range([ d3.symbolCircle, d3.symbolCross, d3.symbolDiamond, d3.symbolSquare, d3.symbolTriangle ]);
//var sValue = function(d) { return d[3];};

           
var valueline = d3.line()
	.x(function(d) { return xscale(d[0]); })
	.y(function(d) { return yscale(d[1]); });
		
for ( i = 0; i < rdata.length; i++ )
{
	var rline = "M" + xscale( rdata[i][0][0] ) + " " + yscale( rdata[i][1][0] ) + " L" + xscale( rdata[i][0][1] ) + " " + yscale( rdata[i][1][1] );	
	chart.append( "path" ).attr("d", rline).style("stroke-width", "2px").style("stroke","black" );
	if ( rdata[i][2].length > 1 )
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
			chart.append( "path" ).attr("d", r2step).style("stroke-width", fs * 0.2 + "px" ).style("stroke","black" ).style( "fill", "none" );
			//console.log( r2step );
		}
		//console.log( r2line )
		//chart.append( "path" ).attr("d", r2line).style("stroke-width", fs * 0.2 + "px" ).style("stroke","black" ).style( "fill", "none" );
	}
}
	
var symbol = d3.symbol();

	

for ( k = 0; k < ldata.length; k++ )
{
	chart.selectAll("dot")
		.data( ldata[k][0] )
		.enter().append("path")
		.attr("r", 5.)
		.attr("transform", function(d,i) { return "translate(" + xscale( d ) + "," + yscale( ldata[k][1][i] ) + ")"; })
		.attr("d", symbol.type( shapeScale( ldata[k][3] ) ).size( 0.5 * fs * fs ) ) 
		.attr("cx", function( d ) { return xscale( d ); } )
		.attr("cy", function( d, i ) { return yscale( ldata[k][1][i] ); } )
		.style("fill", color( ldata[k][2] ) )
		;

	var dline = "M " + xscale( ldata[k][0][0] ) + " " + yscale( ldata[k][1][0] );
	for ( i = 0; i < ldata[k][0].length; i++ )
	{
		dline += " L " + xscale( ldata[k][0][i] ) + " " + yscale( ldata[k][1][i] );
	}
	chart.append( "path" ).attr("d", dline ).style( "stroke-width", fs * 0.1 + "px" ).style( "stroke", color( ldata[ k ][2] ) ).style( "fill", "none" );
}

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
	.text( xaxtitle );

chart.append("text")
	.attr("transform", "rotate(-90)")
	.attr("y", 0 - margin.left + fs * 5 )
	.attr("x",0 - (height / 2))
	.attr("dy", "1em")
	.style("text-anchor", "middle")
	.style("font", 2 * fs + "px sans-serif" )
	.text( yaxtitle );      

var l2t = 0;
var lheight = fs * 2
if ( clabels.length > 1 )
{
	var legend = chart.selectAll(".legend")
		.data( color.domain() )
		.enter().append("g")
			.attr("class", "legend")
			.attr("transform", function(d, i) { return "translate(0," + i * lheight + ")"; });

	legend.append("rect")
		.attr("x", width - 0.9 * lheight )
		.attr("width", 0.9 * lheight )
		.attr("height", 0.9 * lheight )
		.style("fill", color); 

	legend.append("text")
		.attr("x", width - lheight * 1.2 )
		.attr("y", lheight * 0.6 )
		.style("text-anchor", "end")
		.style("font", tfont )
		.text(function(d) { return d;});
	
	l2t = color.domain().length;
}

if ( slabels.length > 1 )
{
	var slegend = chart.selectAll(".slegend")
		.data( shapeScale.domain() )
		.enter().append("g")
			.attr("class", "slegend")
			.attr("transform", function(d, i) { return "translate(0," + ( i + l2t ) * lheight + ")"; });

	slegend.append("path")
		.attr('stroke', 'black')
		.attr('stroke-width', 1)
		.attr('transform', 'translate(' + ( width - lheight * 0.4 ) + ',' + lheight * 0.4  + ')')
		.attr("d", symbol.type( function(d,i){ 	return shapeScale( d ); } ) );
	
	slegend.append("text")
		.attr("x", width - lheight * 1.2 )
		.attr("y", lheight * 0.6 )
		.style("text-anchor", "end")
		.style("font", tfont )
		.text(function(d) { return d;});
}


"""


print "</script>"

