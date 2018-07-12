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
#from sklearn import manifold
#from sklearn import decomposition
#import scipy.stats as stats
#from scipy.spatial import distance
import json
import os.path
import collections
import math

import d3bf

form = cgi.FieldStorage()
id = "emap"
d3bf.chdir( form.getvalue( "datapath" ) )
dgroup = form.getvalue( "dgroup", "none" )
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
resolution = form.getvalue( "resolution", "low" )
numbest = int( form.getvalue( "numbest", "30" ) )
samples = form.getlist( "samples" )
ilevel = int( form.getvalue( "level" ) )

( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
( findex, gtags ) = d3bf.processtags( volumes, tags, dfilter, dgroup )
( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, ilevel )
edata = d3bf.load_edata( data, ilevel, ml, kdict, findex, gtags )
aedata = np.array( edata, dtype=float )
aenorm = np.sum( aedata, axis=1 )
aedata /= aenorm.reshape( len(edata), 1 )

if len( samples ) != 3:
	print "Exactly 3 samples should be selected"
	sys.exit( 1 )
s = [ [], [], [] ]
s[0] = aedata[ gtags[ samples[0] ] ]
s[1] = aedata[ gtags[ samples[1] ] ]
s[2] = aedata[ gtags[ samples[2] ] ]

ssum = s[0] + s[1] + s[2]
sssum = sorted( ssum, reverse=True )
minv = 0
if len( sssum ) > numbest:
	minv = sssum[ numbest ]
	
sval = [ [], [], [] ]
snames = []
sgnames = []
for k in range( len( ssum ) ):
	if ssum[k] < minv:
		continue
	for sk in range( 3 ):
		sval[sk].append( s[sk][k] )
	snames.append( kdnames[ knorder[ k ] ] )
	sgnames.append( kgnames[ knorder[ k ] ] )

d3bf.print_popupstyle()
	
if resolution == "high":
	print "<svg width=\"2400\" height=\"2400\" id=\"normal\"></svg>"
else:
	print "<div id=\"tt\" class=\"tooltip\" style=\"opacity:0;\"></div><svg width=\"960\" height=\"960\" id=\"normal\"></svg>"
print "<script>"
print "var data1 =  %s;" % json.dumps( sval[0] )
print "var data2 =  %s;" % json.dumps( sval[1] )
print "var data3 =  %s;" % json.dumps( sval[2] )
print "var alabels = %s;" % json.dumps( snames ) 
print "var agroups = %s;" % json.dumps( [ gn if len( gn ) > 0 else "Unassigned" for gn in sgnames ] ) 
print "var tlabels = %s;" % json.dumps( samples ) 

print """
	
	
(function() {
var svg0 = d3.select( "#normal" );
var margin = {top: 100, right: 15, bottom: 100, left: 60};
var diameter = +svg0.attr("width") - margin.top - margin.bottom;
var width = diameter;
var svg = svg0.append("g")
var fs = width / 110.; 

    var w = diameter * 0.6;
    var h = diameter * 0.52;
    var m = diameter * 0.1;

	var corners = [[m,h+m], [w+m,h+m], [(w/2)+m,m]]
	var tfont = fs * 2 + "px sans-serif";
	var ltfont = "bold " + fs * 2.6 + "px sans-serif";

	corners.forEach(function(corner, idx) { 
		var c1 = idx, c2 = idx + 1; if(c2 >= corners.length) { c2 = 0;}
		svg.append("line")
			.attr("x1", corners[c1][0])
			.attr("y1", corners[c1][1])
			.attr("x2", corners[c2][0])
			.attr("y2", corners[c2][1])
			.attr('stroke', 'black')
			.attr("stroke-width", fs * 0.15 +"px")
			.attr('fill', 'none');
		svg.append("text")
			.attr( "x", corners[c1][0] + ( idx == 0 ? 0 : ( idx == 1 ? 0 : 0 ) ) )
			.attr( "y", corners[c1][1] + ( idx == 0 ? 5 * fs : ( idx == 1 ? 5 * fs : -3 * fs ) ) )
			.text( tlabels[ idx ] )
			.attr("text-anchor", idx == 0 ? "end" : "start" )
			.style("font", ltfont );
			
	})

	var ticks = [0,20,40,60,80,100], n = ticks.length;
	ticks.forEach(function(v) {
		
		var coord1 = coord(v, 0, 100-v);
		var coord2 = coord(v, 100-v, 0);
		var coord3 = coord(0, 100-v, v);
		var coord4 = coord(100-v, 0, v);

		if(v !== 0 && v !== 100) {

			svg.append("line")
				.attr("x1", coord1[0])
				.attr("y1", coord1[1])
				.attr("x2", coord2[0])
				.attr("y2", coord2[1])
				.attr('stroke', 'black')
				.attr("stroke-width", fs * 0.1 +"px")
				.attr('fill', 'none');

			svg.append("line")
				.attr("x1", coord2[0])
				.attr("y1", coord2[1])
				.attr("x2", coord3[0])
				.attr("y2", coord3[1])
				.attr('stroke', 'black')
				.attr("stroke-width", fs * 0.1 +"px")
				.attr('fill', 'none');

			svg.append("line")
				.attr("x1", coord3[0])
				.attr("y1", coord3[1])
				.attr("x2", coord4[0])
				.attr("y2", coord4[1])
				.attr('stroke', 'black')
				.attr("stroke-width", fs * 0.1 +"px")
				.attr('fill', 'none');
		}
		svg.append("text")
			.attr("x", coord1[0] - 3 * fs)
			.attr("y", coord1[1]  )
			.text( function (d) { return v; })
			.style("font", tfont );
  
		svg.append("text")
			.attr("x", coord2[0] - fs)
			.attr("y", coord2[1] + 2 * fs )
			.text( function (d) { return (100- v); })
			.style("font", tfont );

		svg.append("text")
			.attr("x", coord3[0] + fs )
			.attr("y", coord3[1] )
			.text( function (d) { return v; })
			.style("font", tfont );

	});
	
	var color = d3.scale.category10();
	var cdata = [];
	for ( var i = 0; i < alabels.length; i++ )
	{
		cdata.push( point( data1[i], data2[i], data3[i], alabels[i], agroups[i] ) );
	}
	var clabels = d3.set( cdata.map(function(d) { return d[3];})).values();

	var circles = svg.selectAll("circle").data( cdata );
	var div = d3.select("#tt");
	var ssize = fs * 0.8;

	circles.enter().append("circle")
		.attr("cx", function (d) { return d[0]; })
		.attr("cy", function (d) { return d[1]; })
		.style("fill", function(d) { return color( d[3] ); } )
		.attr("r", ssize)
		.on("mouseover", function(d) {		
			div.transition()		
				.duration(200)		
				.style("opacity", .9);		
			div	.html( d[2] )	
				.style("left", (d3.event.pageX) + "px")		
				.style("top", (d3.event.pageY - 28) + "px");	
			})					
		.on("mouseout", function(d) {		
			div.transition()		
				.duration(500)		
				.style("opacity", 0);	
		});
	
	/*
	circles.enter().append("text")
		.attr("x", function (d) { return d[0]; })
		.attr("y", function (d) { return d[1]; })
		.text( function(d) { return d[2]; })
		.classed('label-text', true);
	*/

	function coord( a, b, c ){
		var sum, pos = [0,0];
	    sum = a + b + c;
	    if(sum !== 0) {
		    a /= sum;
		    b /= sum;
		    c /= sum;

			pos[0] =  corners[0][0]  * a + corners[1][0]  * b + corners[2][0]  * c;
			pos[1] =  corners[0][1]  * a + corners[1][1]  * b + corners[2][1]  * c;
		}
	    return pos;
	}

	function point(a, b, c, label, clabel ){
		var sum, pos = [0,0];
	    sum = a + b + c;
	    if(sum !== 0) {
		    a /= sum;
		    b /= sum;
		    c /= sum;

			pos[0] =  corners[0][0]  * a + corners[1][0]  * b + corners[2][0]  * c;
			pos[1] =  corners[0][1]  * a + corners[1][1]  * b + corners[2][1]  * c;
			pos[2] = label;
			pos[3] = clabel;
		}
	    return pos;
	}
	function scale(/* point */ p, factor) {
	    return [p[0] * factor, p[1] * factor];
	}
	
	var lheight = fs * 2;
	
	var legend = svg.selectAll(".legend")
		.data( color.domain() )
		.enter().append("g")
			.attr("class", "legend")
			.attr("transform", function(d, i) { return "translate(0," + ( i + 4 ) * lheight + ")"; });

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
	


})()

 	
 </script>
"""	