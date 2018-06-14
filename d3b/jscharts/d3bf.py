
import sys
import cgi
import cgitb
cgitb.enable()
import csv
import numpy as np
import json
import os.path
import collections
from operator import itemgetter
import scipy.stats as stats
from scipy.spatial import distance
import math
import re
import skbio
import skbio.diversity as skdiv

def sorted_alnum( l ):
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)


def loaddata( fname ):
	reader = csv.reader(open( fname, "r"), delimiter='\t')
	data = list(reader)
	volumes = []
	for i in range( 0, len( data[0] ) ):
		if len( data[0][i] ) > 0:
			mn = i
			c0 = data[0][i][0]
			if c0 >= '0' and c0 <= '9':
				ml = i
			else:
				volumes.append( data[0][i] )
	return ( data, volumes, mn, ml )


def loadtags( tfname, volumes ):
	tags = {}
	if os.path.isfile( tfname ):
		tfreader = csv.reader(open( tfname, "r"), delimiter='\t')
		tfdata = list( tfreader )
		vkindex = []
		for j in range( 1, len( tfdata ) ):
			vkindex.append( tfdata[ j ][ 0 ] )
		for i in range( len( tfdata[0] ) - 1 ):
			tval = [ "" ] * len( volumes )
			for j in range( 1, len( tfdata ) ):
				vind = volumes.index( vkindex[ j - 1 ] )
				tval[ vind ] = tfdata[ j ][ i ] 
			tags[ tfdata[0][i] ] = tval
		tags[ "none" ] = [ "" ] * len( volumes )
	else:
		tags[ "name" ] = volumes
		tags[ "none" ] = [ "" ] * len( volumes )
	tkeys = sorted_alnum( tags.keys() )
	return ( tags, tkeys )

def loadfilters( sfname, spfilter ):
	if os.path.isfile( sfname ):
		sfreader = csv.reader(open( sfname, "r"), delimiter='\t')
		sfdata = list( sfreader )
		for j in range( 0, len( sfdata ) ):
			if sfdata[j][0] == spfilter and len( sfdata[j] ) == 3:
				reverse = int( sfdata[j][1] )
				splist = sfdata[j][2].split( ";" )
				splist.remove( "" )
				return ( splist, reverse )
	return ( [], 0 )

def processtags( volumes, tags, dfilter, dgroup ):
	findex = {}
	gtags = {}
	gcnt = 0
	for j in range( len( volumes ) ):
		if dfilter == "none" or len( tags[ dfilter ][ j ] ):
			cgtag = volumes[j]		
			if dgroup != "none":
				cgtag = tags[ dgroup ][ j ]
			if not cgtag in gtags:
				gtags[cgtag] = gcnt
				findex[j] = cgtag
				gcnt = gcnt + 1
			else:
				findex[j] = cgtag
	return ( findex, gtags )

def processtags_m( volumes, tags, dfilter ):
	if dfilter == "none":
		findex = range( 0, len( volumes ) )
		mtags = tags
	else:
		findex = []
		for j in range( len( volumes ) ):
			if len( tags[ dfilter ][ j ] ):
				findex.append( j )
		mtags = {}
		for tkey in tags:
			mtval = []
			for j in range( len( volumes ) ):
				if j in findex: 				
					mtval.append( tags[tkey][ j ] ) 
			mtags[ tkey ] = mtval
	return ( findex, mtags )

def loadtaxonomy( data, ml, spfilter, ilevel ):
	kdict = {}
	kdnames = {}
	kgnames = {}
	cm = 0
	knorder = []
	kdata = {}
	if isinstance( spfilter, (list,)):
		excludelist = []
		reverse = 0
	else:
		excludelist = spfilter[0]
		reverse = spfilter[1]
	dgnames = []
	minglevel = 0 if ml < 2 else 1
	maxglevel = max( 4, ilevel - 2 )
	nkgnames = [ {} for x in xrange( minglevel, maxglevel ) ]
	ndgnames = [ [] for x in xrange( minglevel, maxglevel ) ]
	
	for d in data[1:]:
		if len( d ) < len( data[0] ):
			continue
		ckey = "".join( d[0:ilevel] )
		ckname = ""
		cgname = d[0]
		exflag = reverse
		if len( excludelist ) > 0:
			for k in range( ml ):
				if len( d[k] ) > 0 and d[k] in excludelist:
					exflag = 1 - reverse
					break
		if exflag == 1:
			continue
		for k in range( ilevel - 1, -1, -1 ):
			if len( d[k] ) > 0 and d[k] != "*" and d[k] != "Unknown":
				ckname = d[k]
				break
		if ml > 0 and len( d[1] ) > 0:
			cgname = d[1]
		if not ckey in kdict:
			kdict[ ckey ] = cm
			kdnames[ ckey ] = ckname
			kgnames[ ckey ] = cgname
			for k in range( maxglevel - minglevel ):
				dkey = d[ k + minglevel ]
				nkgnames[ k ][ ckey ] = dkey
				if not dkey in ndgnames[k]:
					ndgnames[k].append( dkey )
			cm = cm + 1
			knorder.append( ckey )
			kdata[ ckey ] = d[0:ilevel] 
	for k in range( 0, maxglevel - minglevel ):
		if k + minglevel < len( ndgnames ) and len( ndgnames[ k + minglevel ] ) > 2:
			break
	return ( kdict, kdnames, nkgnames[k], knorder, kdata )

def select_toptax( edata, kdict, num_best ):
	if num_best > 0 and num_best < len( edata[0] ):
		aedata = np.array( edata, dtype=float )
		aenorm = np.sum( aedata, axis=1 )
		aedata /= aenorm.reshape( len(edata), 1 )
		ssum = np.sum( aedata, axis=0 )
		ssorted = sorted( ssum.tolist(), reverse=True )
		smax = ssorted[ num_best ]
		nkdict = {}
		nedata = []
		tcnt = 0
		for key in kdict:
			if ssum[ kdict[key] ] > smax and tcnt < num_best:
				for k in range( len( edata ) ):
					if tcnt == 0:
						nedata.append( [] )
					nedata[k].append( edata[k][ kdict[key] ] )
				nkdict[ key ] = tcnt
				tcnt += 1
		return ( nedata, nkdict )
	else:
		return ( edata, kdict )
		
def load_edata( data, ilevel, ml, kdict, findex, gtags ):
	edata = []
	for tagnum in sorted( gtags.values() ):
		cgtag = gtags.keys()[ gtags.values().index( tagnum ) ]
		crow = [ 0 ] * len( kdict )
		for i in sorted( findex ):
			if findex[i] == cgtag:
				for d in data[1:]:
					if len( d ) < max( findex.keys() ) + ml:
						continue
					ckey = "".join( d[0:ilevel] )
					if ckey in kdict:
						cind = kdict[ ckey ]
						v = int( float( d[ i + ml + 1 ] ) )
						crow[ cind ] += v
		edata.append( crow )
	return edata

def load_edata_m( data, ilevel, mn, ml, kdict, volumes, findex, kdnames ):
	edata = []
	site_ids = []
	species_ids = [ "" ] * len( kdict )
	dim2_ids = species_ids
	for i in range( ml + 1, mn + 1 ):
		if ( i - ml - 1 ) in findex:
			crow = [ 0 ] * len( kdict )
			site_ids.append( volumes[ i - ml - 1 ] )
			for d in data[1:]:
				ckey = "".join( d[0:ilevel] )
				if ckey in kdict:
					cind = kdict[ ckey ]
					crow[ cind ] += int( float( d[i] ) )
					if len( site_ids ) == 1:
						species_ids[ cind ] = kdnames[ ckey ] 
			edata.append( crow )
	return ( edata, site_ids, species_ids )


	
def morisitaHorn(data1, data2):
    X = sum(data1)
    Y = sum(data2)

    if len(data1) != len(data2):
        raise ValueError("Error while calculating MorisitaHorn similarity index. The two input data lists must have the same length!\n")
    
    if X == 0 or Y == 0:
        return -1

    sumXiYi = 0
    sumXiSqr = 0
    sumYiSqr = 0

    for i, x in enumerate(data1):
        y = data2[i]
        sumXiYi += x*y
        sumXiSqr += x*x
        sumYiSqr += y*y

    numerator = 2*sumXiYi
    denominator = ( float(sumXiSqr)/(X*X) + float(sumYiSqr)/(Y*Y) )*X*Y

    return float(numerator)/denominator

def calc_distances( edata, aedata, dmethod, kdata, knorder, rev ):	
	if dmethod ==  "Jaccard" or dmethod.find( "Unifrac" ) != -1:
		iedata = np.array( edata, dtype=int )
	if dmethod.find( "Unifrac" ) != -1:
		mkdata = {}
		for key in kdata:
			mkdata[ key ] = [ "Root" ] + kdata[ key ]
		tree = skbio.tree.TreeNode.from_taxonomy( mkdata.items() )
		for node in tree.preorder():
			node.length = 1.
	cdata = []
	irev = 0
	imax = 0
	eps = 1e-7
	for i in range( len( edata ) ):
		cdrow = [ 0 ] * len( edata )
		for j in range( i ):
			if dmethod in [ "Pearson", "Kendall", "Spearman" ]:
				if dmethod == "Pearson":
					corr, pvalue = stats.pearsonr( edata[i], edata[j] )
				elif dmethod == "Kendall":
					corr, pvalue = stats.kendalltau( edata[i], edata[j] )
				elif dmethod == "Spearman":
					corr, pvalue = stats.pearsonr( edata[i], edata[j] )
				v = 1 - corr
				imax = 1.
			elif dmethod == "Euclidean":
				v = distance.euclidean( aedata[i], aedata[j] )
			elif dmethod == "Bray-Curtis":
				v = distance.braycurtis( edata[i], edata[j] )
			elif dmethod == "Morisita-Horn":
				v = 1 - morisitaHorn( edata[i], edata[j] )
			elif dmethod == "Jaccard":
				bd1 = ( iedata[i] > 0 )
				bd2 = ( iedata[j] > 0 )
				v = distance.jaccard( bd1, bd2 )
			elif dmethod == "Unifrac-weighted":
				v = skdiv.beta.weighted_unifrac( iedata[i], iedata[j], knorder, tree, validate=True )
			elif dmethod == "Unifrac-unweighted":
				v = skdiv.beta.unweighted_unifrac( iedata[i], iedata[j], knorder, tree, validate=True )
			if math.isnan( v ):
				v = 1
			v = max( v, 0 )
			imax = max( imax, v )
			cdrow[ j ] = v 
		cdrow[ i ] = 0
		cdata.append( cdrow )
	#print rev
	for i in range( len( cdata ) ):
		for j in range( i ):
			if rev == 1:
				cdata[i][j] = ( imax - cdata[i][j] ) / imax
			cdata[j][i] = cdata[i][j]
#	print cdata
	return cdata

def print_popupstyle():
	print """
<style>
.axis line, .axis path {
    shape-rendering: crispEdges;
    stroke: black;
    fill: none;
}

circle {
    fill: steelblue;
}

div.tooltip {	
    position: absolute;			
    text-align: center;			
    width: 100px;					
    height: 16px;					
    padding: 2px;				
    font: 12px sans-serif;		
    background: #f4f2ec;	
    border: 1px solid #7E713D;		
    border-radius: 8px;			
    pointer-events: none;			
}
</style>
"""

	
def print_drawscatter():
	print """
	function drawscatter( val_id, data )
	{
		var margin = {top: 20, right: 15, bottom: 60, left: 60};
		var svg = d3.select( val_id )
			.attr('background-color', 'white' );
		var width = +svg.attr( "width" ) - margin.right - margin.left;
		var height = +svg.attr( "height" ) - margin.top - margin.bottom;
		var fs = Math.max( width / 100, 8 );
		var chart = svg.append("g")
				.attr('background-color', 'white' )
				.attr("transform", "translate(" + margin.left + "," + margin.top + ")");
    
		var xscale = d3.scale.linear().range([ 0, width ]);
		var yscale = d3.scale.linear().range([ height, 0 ]);
		var xAxis = d3.svg.axis().scale(xscale).orient('bottom').ticks(2);
		var yAxis = d3.svg.axis().scale(yscale).orient('left').ticks(2);
		var xValue = function(d) { return d[1];}; 
		var yValue = function(d) { return d[0];};
		var cValue = function(d) { return d[2];};
		var color = d3.scale.category10();
		var vmax = Math.max( -d3.min(data, xValue), d3.max(data, xValue), -d3.min(data, yValue), d3.max(data, yValue) ) * 1.6;
		xscale.domain([-vmax, vmax]);	
		yscale.domain([-vmax, vmax]);	
		xMap = function(d) { return xscale(xValue(d)); };
		yMap = function(d) { return yscale(yValue(d)); };
		var sValue = function(d) { return d[4];};
		var slabels = d3.set(data.map(function(d) { return d[4];})).values();	
		var clabels = d3.set(data.map(function(d) { return d[2];})).values();
		var shapeScale = d3.scale.ordinal()
            .domain(slabels)
            .range([ "circle", "cross", "diamond", "square", "triangle-down", "triangle-up" ]);
		var tfont = fs * 1.4 + "px sans-serif";
		var ssize = fs * fs * 0.7;
            
   	    

		chart.append('g')
			.attr('transform', 'translate(0,' + height / 2 + ')')
			.attr('class', 'x axis')
			.attr('shape-rendering', 'crispEdges')
			.attr('stroke', 'black')
			.attr("stroke-width", fs * 0.15 +"px")
			.attr('fill', 'none')	
			.call(xAxis)
			.append("text")
				.attr("class", "label")
				.attr("x", width)
				.attr("y", - fs * 0.6 )
				.style("text-anchor", "end")
				.style("font", tfont )
				.text( atitle + " 1 " + alabels[0] );	

		chart.append('g')
			.attr('transform', 'translate(' + width / 2 + ',0)')
			.attr('class', 'y axis')
			.attr('shape-rendering', 'crispEdges')
			.attr('stroke', 'black')
			.attr("stroke-width", fs * 0.15 +"px")
			.attr('fill', 'none')	
			.call(yAxis)
			.append("text")
				.attr("class", "label")
				.attr("transform", "rotate(-90)")
				.attr("y", fs * 1.4 )
				.style("text-anchor", "end")
				.style("font", tfont )
				.text( atitle + " 2 " + alabels[1] );

		
		chart.selectAll(".dodo")
		.data(data)
		.enter().append("text")
			.attr("class", "dodo")
			.attr("x", xMap)
			.attr("y", yMap)
			.style("font", tfont )
			.text(function(d) { return d[3];});
			
		if ( data[0].length > 5 )
		{
            var div = d3.select("#tt");
            
                    
			chart.selectAll("dot")
			.data(data)
			.enter().append("circle")
				.attr("cx", function(d) { return xMap(d); } )	
				.attr("cy", function(d) { return yMap(d); } )	
				.attr("r", 5 )
				.style("fill", function(d) { return color( cValue( d ) ); } )
				.on("mouseover", function(d) {		
					div.transition()		
						.duration(200)		
						.style("opacity", .9);		
					div	.html( d[5] )	
						.style("left", (d3.event.pageX) + "px")		
						.style("top", (d3.event.pageY - 28) + "px");	
					})					
				.on("mouseout", function(d) {		
					div.transition()		
						.duration(500)		
						.style("opacity", 0);	
				});
		}
		else
		{
			chart.selectAll(".point")
			.data(data)
			.enter().append("path")
				.attr("class", "point")
				.attr("transform", function(d) { return "translate(" + xMap(d) + "," + yMap( d ) + ")"; })
				.attr("d", function(d,i){ 
					return d3.svg.symbol().type( shapeScale( d[4] ) ).size( ssize )();
					} )
				.style("fill", function(d) { return color( cValue( d ) ); } )
				;
		}
			

    
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
				.attr("d", function(d,i){ 
					return d3.svg.symbol().type( shapeScale( d ) ).size( ssize )();
					} );
			
			slegend.append("text")
				.attr("x", width - lheight * 1.2 )
				.attr("y", lheight * 0.6 )
				.style("text-anchor", "end")
				.style("font", tfont )
				.text(function(d) { return d;});
		}
	}
"""






