#!/usr/bin/python

# Turn on debug mode.
import sys
import cgi
import cgitb
cgitb.enable()
import csv
import numpy as np
import json
import os.path
import collections
import scipy.stats as stats

import d3bf

form = cgi.FieldStorage()
id = "emap"
d3bf.chdir( form.getvalue( "datapath" ) )
dgroup = form.getvalue( "dgroup", "none" )
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
level = form.getvalue( "level" )
numbest = form.getvalue( "numbest", "inf" )
resolution = form.getvalue( "resolution", "low" )
cmethod = form.getvalue( "cmethod", "Pearson" )
cthreshold = form.getvalue( "cthreshold", "0.5" )
psizes = form.getvalue( "psizes", "equal" )

( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
( findex, gtags ) = d3bf.processtags( volumes, tags, dfilter, dgroup )
ilevel = int( level ) 

( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, ilevel )

edata = d3bf.load_edata( data, ilevel, ml, kdict, findex, gtags )

aedata = np.array( edata, dtype=float )
aenorm = np.sum( aedata, axis=1 )
aedata /= aenorm.reshape( len(edata), 1 )
	
inumbest = 0
if numbest.isdigit():
	inumbest = int( numbest )
if ( inumbest > 0 and inumbest < len( kdict ) ):
	( nedata, aedata, nkdict ) = d3bf.select_toptax( edata, kdict, num_best=inumbest )
	edata = nedata
	nknorder = []
	for cknum in range( len( knorder ) ):
		ckey = knorder[ cknum ]
		if ckey in nkdict:
			nknorder.append( ckey )
	kdict = nkdict
	knorder = nknorder
	
edges = []
nodes = []
ordlist = [] 

for ind1 in range( len( kdict )):
	ckey = knorder[ ind1 ]
	row1 = []
	for cgtag in list( gtags.keys() ):
		row1.append( aedata[ gtags[ cgtag ] ][ ind1 ] )
	dratio = sum( row1 ) / len( row1  )
	nodes.append( { "id" : ind1, "label" : kdnames[ ckey ], "module" : kgnames[ ckey ], "abundance": dratio } )
	if not kgnames[ ckey ] in ordlist and dratio > 0.01:
		ordlist.append( kgnames[ ckey ] )
	for ind2 in range( ind1 ):
		row2 = []
		for cgtag in list( gtags.keys() ):
			row2.append( aedata[ gtags[ cgtag ] ][ ind2 ] )
		if cmethod == "Pearson":
			corr, pvalue = stats.pearsonr( row1, row2 )
		elif cmethod == "Spearman":
			corr, pvalue = stats.pearsonr( row1, row2 )
		elif cmethod == "Kendall":
			corr, pvalue = stats.kendalltau( row1, row2 )
		elif cmethod == "Binary":
			corr, pvalue = stats.pearsonr( np.sign( np.array( row1 ) ), np.sign( np.array( row2 ) ) )
		if corr > float( cthreshold ):
			edges.append( { "id" : str( ind1 ) + "-" + str( ind2 ), "source" : ind1, "target" : ind2, "weight" : 1 } )
			
d3bf.print_popupstyle()

if resolution == "high":
	print("<svg width=\"2400\" height=\"2400\" id=\"normal\"></svg>")
else:
	print("<div id=\"tt\" class=\"tooltip\" style=\"opacity:0;\"></div><svg width=\"960\" height=\"960\" id=\"normal\"></svg>")
print("<script>")
print( "var nodes =  %s;" % json.dumps( nodes ) )
print( "var edges =  %s;" % json.dumps( edges ) )
print( "var ordlist = %s;" % json.dumps( ordlist ) )
print( "var psizes = %s;" % json.dumps( psizes ) )

print("""

var svg0 = d3.select( "#normal" );
var margin = {top: 100, right: 15, bottom: 100, left: 60};
var diameter = +svg0.attr("width") - margin.top - margin.bottom;
var width = diameter;
var svg = svg0.append("g")
var ssize = width * 0.05;
var color = d3.scale.category20();
var tooltip = d3.select("#tt");

var force = d3.layout.force().nodes(nodes).links(edges)
	.size( [ width, width ] )
	.charge( -0.3 * width )
	.chargeDistance( 0.3 * width )
	.linkStrength( 0.1 * ssize )
	.linkDistance( 3 * ssize );
	//.gravity(0.07);
		
var edgeEnter = svg.selectAll("g.edge")
	.data(edges, function (d) {return d.id})
	.enter()
	.append("g")
	.attr("class", "edge");

edgeEnter.append("line")
	.style( "stroke-width", 0.02 * ssize + "px" )
	.style( "stroke", "gray")
	.style( "pointer-events", "none" );

var nodeEnter = svg.selectAll("g.node")
	.data(nodes, function (d) {return d.id})
	.enter();

nodeEnter.append("circle")
	.attr( "r", function (d) { return ( ( psizes == "equal" ) ? ssize * 0.3 : ssize * Math.sqrt( d.abundance ) ) + "px"; } )
	.style( "fill", function (d) { return color( d.module ); } )
	.on("mouseover", function(d) {		
		tooltip.transition()		
			.duration(200)		
			.style("opacity", .9);		
		tooltip.html( d.label )	
			.style("left", (d3.event.pageX) + "px")		
			.style("top", (d3.event.pageY - 28) + "px");	
		})					
	.on("mouseout", function(d) {		
		tooltip.transition()		
			.duration(500)		
			.style("opacity", 0);	
	});

	
force.start();
for( var i = 0; i < 200; i++ ) force.tick();
force.stop();
updateNetwork();

console.log( "force:", force )

function updateNetwork(e) {

	svg.selectAll("line")
		.attr("x1", function (d) {return d.source.x})
		.attr("y1", function (d) {return d.source.y})
		.attr("x2", function (d) {return d.target.x})
		.attr("y2", function (d) {return d.target.y});

	svg.selectAll("circle")
		.attr("cx", function (d) {return d.x; })
		.attr("cy", function (d) {return d.y; });
		
}

var dlegend = "yes";
var fs = diameter / 100;
if ( dlegend == "yes" )
{
	var cellwidth = fs * 11;
	var cellheight = fs * 1.3;
	var height = diameter + cellheight * 4;
	svg.selectAll(".legendtext")
		.data( ordlist )
		.enter().append("text")
		.attr("class", "legendtext")
		.attr("x", function(d,i) { return i * cellwidth; } )
		.attr("y", function(d,i) { return height + 5 * cellheight + ( i % 3 ) * fs * 1.8; } )
		.style("font", 1.6 * fs + "px sans-serif" )
		.style("text-anchor", "start" ) 
		.text(function(d) { return d;}
		);      

	svg.selectAll( ".legendrect")
		.data( ordlist )
		.enter().append("rect")
		.attr("class", ".legendrect")
		.attr("x", function(d,i) { return i * cellwidth; })
		.attr("y", height + 2 * cellheight )
		.attr("width", cellwidth )
		.attr("height", cellheight )
		.style("stroke",  "black")
		.style("stroke-width", "0.5px" )
		.style("fill", function(d) { return color( d ); } 
		);
}

</script>
""")


