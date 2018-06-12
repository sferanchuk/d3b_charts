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
from sklearn import manifold
from sklearn import decomposition
from sklearn import preprocessing
import scipy.stats as stats
from scipy.spatial import distance
import json
import os.path
import sqlite3
import time
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import skbio.stats as skbiostats
from skbio import DistanceMatrix
import math
import skbio.stats.ordination as ordination
import collections

import d3bf

form = cgi.FieldStorage()
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
level = form.getvalue( "level" )
fmethod = form.getvalue( "fmethod", "PCA" )
color = form.getvalue( "color", "none" )
shape = form.getvalue( "shape", "none" )
labels = form.getvalue( "labels", "none" )
resolution = form.getvalue( "resolution", "low" )
pc1 = form.getvalue( "pc1", "1" )
pc2 = form.getvalue( "pc2", "2" )
ptype = "species"

ilevel = int( level ) 
( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, ilevel )
( findex, mtags ) = d3bf.processtags_m( volumes, tags, dfilter )
( edata, site_ids, species_ids ) = d3bf.load_edata_m( data, ilevel, mn, ml, kdict, volumes, findex, kdnames )


aedata = np.array( edata, dtype=float )
aenorm = np.sum( aedata, axis=1 )
aedata /= aenorm.reshape( len(edata), 1 )
ncomp = min( 5, len( volumes ) )
ptypes = "abundancies"

if False:
	import sys
	tfp = '/home/sferanchuk/.local/lib/python2.7/site-packages'
	if tfp in sys.path:
		sys.path.remove( tfp )
	import cgi
	import cgitb
	cgitb.enable()
	import csv
	import numpy as np
	from sklearn import manifold
	from sklearn import decomposition
	from sklearn import preprocessing
	import scipy.stats as stats
	from scipy.spatial import distance
	import json
	import os.path
	import collections
	import sqlite3
	import time
	import base64

	def npencode( ndarray ):
		return json.dumps([str(ndarray.dtype),base64.b64encode(ndarray),ndarray.shape])

	def npdecode( strval ):
		loaded = json.loads( strval )
		dtype = np.dtype( loaded[0] )
		arr = np.frombuffer( base64.decodestring( loaded[1]), dtype )
		if len(loaded) > 2:
			return arr.reshape( loaded[2] )
		return arr


	form = cgi.FieldStorage()
	id = form.getvalue( "id" )
	color = form.getvalue( "color", "none" )
	shape = form.getvalue( "shape", "none" )
	labels = form.getvalue( "labels", "none" )
	dfilter = form.getvalue( "filter", "none" )
	pc1 = form.getvalue( "pc1", "1" )
	pc2 = form.getvalue( "pc2", "2" )
	#fmethod = form.getvalue( "fmethod", "MDS" )
	#dmethod = form.getvalue( "dmethod", "Pearson" )
	level = form.getvalue( "level" )
	excludestr = form.getvalue( "exclude", "" )
	if form.getvalue( "noexclude" ):
		excludestr=""
	excludelist = excludestr.split( "--" )
	excludelist.remove( "" )
	revex = form.getvalue( "revex", "none" )
	ptypes = form.getvalue( "ptypes", "samples" )

	if id:
		reader = csv.reader(open( id + ".txt", "r"), delimiter='\t')
		data = list(reader)
		ml = 0
		mn = 0
		volumes = []
		for i in range( 0, len( data[0] ) ):
			if len( data[0][i] ) > 0:
				mn = i
				c0 = data[0][i][0]
				if c0 >= '0' and c0 <= '9':
					ml = i
				else:
					volumes.append( data[0][i] )
		if not level:
			level = str( min( ml + 1, 4 ) )	
		ilevel = int( level )
		cm = 0
		kdict = {}
		kdnames = {}
		kgnames = {}
		for d in data[1:]:
			ckey = "".join( d[0:ilevel] )
			ckname = ""
			cgname = d[0]
			exflag = 0
			if len( excludelist ) > 0:
				for k in range( ml ):
					dkey = d[k]
					if len( dkey ) > 0 and dkey[0] == "[":
						dkey = dkey[1:]
					if len( dkey ) > 0 and dkey[ -1 ] == "]":
						dkey = dkey[0:-1] 
					if len( dkey ) > 0 and dkey in excludelist:
						exflag = 1
						break
			if exflag == 1:
				if revex == "none":
					continue
			else:
				if len( excludelist ) and revex != "none":
					continue
			for k in range( ilevel - 1, -1, -1 ):
				if len( d[k] ) > 0:
					ckname = d[k]
					break
			if ml > 0 and len( d[1] ) > 0:
				cgname = d[1]
			if not ckey in kdict:
				kdict[ ckey ] = cm
				kdnames[ cm ] = ckname
				kgnames[ cm ] = cgname
				cm = cm + 1
				
#print >> sys.stderr, 'pca_sp 1'		
adist = np.array( edata, dtype=float )
#if ptypes == "species":
#	adist = adist.transpose()
nadist = preprocessing.scale( adist )
pca = decomposition.PCA( n_components = ncomp )
pca.fit( nadist )
coords = np.asarray( pca.transform( nadist ), order='C' )
pcavar = np.asarray( pca.explained_variance_[ 0 : ncomp ], order='C' )
pcacomp = np.asarray( pca.components_[ 0 : ncomp ], order='C' )

#print >> sys.stderr, 'pca_sp 2'		
tadist = nadist.transpose()
tpca = decomposition.PCA( n_components = ncomp )
tpca.fit( tadist )
tcoords = np.asarray( tpca.transform( tadist ), order='C' )
tpcavar = np.asarray( tpca.explained_variance_[ 0 : ncomp ], order='C' )
tpcacomp = np.asarray( tpca.components_[ 0 : ncomp ], order='C' )
	
alabels = [ "(%d) %4.2f" % ( int( pc1 ), pcavar[ int( pc1 ) - 1 ] ), "(%d) %4.2f" % ( int( pc2 ), pcavar[ int( pc2 ) - 1 ] ) ]
splist = []
sglist = []
svlist = []
ssvlist = []
sspec = []
for i in range( ncomp ):
	olist = {}
	for k in range( len( kdict ) ):
		olist[ pow( pcacomp[i][ k ], 2 )] = k
	solist = collections.OrderedDict( sorted( olist.items(), reverse=True ) )
	csplist = []
	csglist = []
	csvlist = []
	cssvlist = []
	for ik, ij in solist.iteritems():
		if len( csplist ) < 20:
			if len( csplist ) < 15 and not ij in sspec:
				sspec.append( ij )
			csplist.append( kdnames[ knorder[ ij] ] )
			csglist.append( kgnames[ knorder[ ij] ] )
			csvlist.append( ik )
			cssvlist.append( pcacomp[i][ ij ] )
	splist.append( csplist )
	sglist.append( csglist )
	svlist.append( csvlist )
	ssvlist.append( cssvlist )
	
spcoords = [ [] for x in xrange( ncomp ) ]
spnames = []
spgroups = []
for k in sspec:
	spnames.append( kdnames[ knorder[k] ] )
	spgroups.append( kgnames[ knorder[k] ] )
	for i in range( ncomp ):
		spcoords[i].append( pcacomp[ i ][ k ] )
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


	
if resolution == "high":
	print "<table><tr><td><svg width=\"1200\" height=\"1200\" id=\"ssamples\"></svg>"
	print "<td><svg width=\"1200\" height=\"1200\" id=\"sspecies\"></svg></table><br>"
	print "<svg width=\"1800\" height=\"1200\" id=\"slists\"></svg>"
else:
	print "<table><tr><td><svg width=\"600\" height=\"600\" id=\"ssamples\"></svg>"
	print "<td><div id=\"tt\" class=\"tooltip\" style=\"opacity:0;\"></div><svg width=\"600\" height=\"600\" id=\"sspecies\"></svg></table><br>"
	print "<svg width=\"1200\" height=\"600\" id=\"slists\"></svg>"

print "<script>"
print "var pc1 = %d; " % int( pc1 )
print "var pc2 = %d; " % int( pc2 )
print "var ncomp = %d; " % ncomp
print "var dataxy = [" 
for i in range( ncomp ):
	print "%s %s" % ( json.dumps( coords[:,i].tolist() ), "," if i + 1 < ncomp else "" )
print "];";
#print "var datay =  %s;" % json.dumps( coords[:, int( pc2 ) - 1 ].tolist() )
print "var datac =  %s;" % json.dumps( mtags[ color ] )
print "var datal =  %s;" % json.dumps( mtags[ labels ] )
print "var datas =  %s;" % json.dumps( mtags[ shape ] )
print "var alabels = %s;" % json.dumps( alabels ) 
print "var atitle = 'PC ';"
print "var tdataxy = [" 
for i in range( ncomp ):
	print "%s %s" % ( json.dumps( spcoords[i] ), "," if i + 1 < ncomp else "" )
print "];";
print "var tdatal = %s; " % ( json.dumps( spnames ) )
print "var tdatac = %s; " % ( json.dumps( spgroups ) )
#print "var tdatax =  %s;" % json.dumps( tcoords[:, int( pc1 ) - 1 ].tolist() )
#print "var tdatay =  %s;" % json.dumps( tcoords[:, int( pc2 ) - 1 ].tolist() )


print "var spdata = [ "
for i in range( len( splist ) ):
	print "{ names: %s, values: %s, signs: %s, pc: \"%5.2f\", ind: %d }" %  ( json.dumps( splist[i] ), json.dumps( svlist[i] ), json.dumps( ssvlist[i] ), pcavar[i], i )
	if i + 1 < len( splist ):
		print ","
print "];"

if False:
	slabels = [ "" ] * len( kdict )
	clabels = [ "" ] * len( kdict )
	nlabels = [ "" ] * len( kdict )
	for ckey in kdict:
		cind = kdict[ ckey ]
		if labels != "none":
			nlabels[ cind ] = kdnames[ cind ]
		if color != "none":
			clabels[ cind ] = kgnames[ cind ]
		if shape != "none":
			slabels[ cind ] = kgnames[ cind ]
	print "var datac =  %s;" % json.dumps( clabels )# var data = [[5,3], [10,17], [15,4], [2,8]];
	print "var datal =  %s;" % json.dumps( nlabels )# var data = [[5,3], [10,17], [15,4], [2,8]];
	print "var datas =  %s;" % json.dumps( slabels )# var data = [[5,3], [10,17], [15,4], [2,8]];


print """
    var pdata = [];
    var i;
    for ( i = 0; i < dataxy[0].length; i++ )
    {
		pdata.push( [ dataxy[ pc1 - 1 ][i], dataxy[ pc2 - 1 ][i], datac[i], datal[i], datas[i] ] );
	}
	
	var tpdata = [];
    for ( i = 0; i < tdataxy[0].length; i++ )
    {
		if ( spdata[ pc1 - 1 ].names.indexOf( tdatal[i] ) == -1 && spdata[ pc2 - 1 ].names.indexOf( tdatal[i] ) == -1 ) continue;
		tpdata.push( [ tdataxy[ pc1 - 1 ][i], tdataxy[ pc2 - 1 ][i], tdatac[i], "", "", tdatal[i] ] );
	}
	
	
"""

d3bf.print_drawscatter()

print "drawscatter( \"#ssamples\", pdata );"
print "drawscatter( \"#sspecies\", tpdata );"
print """
function drawbars()
{
	var chart = d3.select( \"#slists\" )
			.attr('background-color', 'white' );
	var margin = {top: 20, right: 15, bottom: 60, left: 60};
	var c10 = d3.scale.category10();
	var bheight = 40;
	var width = +chart.attr( "width" ) - margin.right - margin.left;
	var gwidth = width / 5;
	var fs = width / 110.; 
	for ( i = 0; i < spdata.length; i++ )
	{
		var bclass = ".dodo" + i;
		var tclass = ".coco" + i;
		var cspdata = [];
		var k;
		var dmax = 0;
		for ( k = 0; k < spdata[i].values.length; k++ )
		{
			var cv = +spdata[i].values[k];
			var signv = ( + spdata[i].signs[k] > 0 ) ? 0 : 1;
			cspdata.push( [ k, cv, spdata[i].names[k], signv ] );
			if ( k == 0 ) dmax = cv;
		}
		chart.selectAll( bclass )
			.data( cspdata )
			.enter().append("rect")
				.attr("class", bclass )
				.attr("x", margin.left + i * gwidth + gwidth * 0.47 )
				.attr("y", function(d) { return bheight + 60 + d[0] * 15; } )
				.attr("width", function(d) { return 80 * d[1] / dmax; } )
				.attr("height", 10 )
				.style("fill", function(d) { if ( d[3] != 0 ) return c10(1); else return c10(0); } );
		chart.selectAll( tclass )
			.data( cspdata )
      .		enter().append("text")
				.attr("class", tclass )
				.attr("x", margin.left + i * gwidth + gwidth * 0.45 )
				.attr("y", function(d) { return bheight + 70 + d[0] * 15; } )
				.attr("text-anchor", "end")
				.text(function(d) { return d[2].substring(0,20);})
				.style("font-size", fs + "px");  
	}
	chart.selectAll( ".bobo" )
			.data( spdata )
      .		enter().append("text")
				.attr("class", ".bobo" )
				.attr("x", function(d) { return margin.left + d.ind * gwidth + gwidth * 0.47;} )
				.attr("y", function(d) { return bheight + 50; } )
				//.attr("text-anchor", "end")
				.text(function(d) { return "PC " + ( d.ind + 1 ) + " (" + d.pc + ")";})
				.style("font-size", fs * 1.2 + "px"); 
}

drawbars();
"""

print "</script>"	
sys.exit( 0 )
		
if False:
	print """
	<a href="emap.php?id=%s">back</a><br>
	<a href="emapsort.php?id=%s&level=%d">table data</a> (to set species filter)<br>
	""" % ( cgi.escape( id ), cgi.escape( id ), ilevel )
	print """
	<form method="get" action="emapd3pca_sp.py">
	<input type=hidden name="id" value="%s"/>
	<input type=hidden name="exclude" value="%s"/>
	""" % ( cgi.escape( id ), cgi.escape( excludestr ) )
	if len( excludelist ):
		print "Filtered species: %s <input type=checkbox name=noexclude value=off>Remove species filter<br>" % ( cgi.escape( "".join( "{0};".format( v ) for v in excludelist ) ) )
	selstr = '<select name="{0}">\n{1}</select>\n'
	optstr = '<option value="{0}" {1} >{0}</option>\n'
	rbstr = '<input type=radio name="{0}" value="{1}" {2}>{1}\n'
	print "<br>Level "
	print selstr.format( "level", ''.join( optstr.format(v, "selected" if v == ilevel else "" ) for v in range( 1, ml + 2 ) ) )
	print "<br>Filter "
	print selstr.format( "filter", ''.join( optstr.format(v, "selected" if v == dfilter else "" ) for v in tkeys ) )
	print "<br>Color "
	print selstr.format( "color", ''.join( optstr.format(v, "selected" if v == color else "" ) for v in tkeys ) )
	print "<br>Shape "
	print selstr.format( "shape", ''.join( optstr.format(v, "selected" if v == shape else "" ) for v in tkeys ) )
	print "<br>Labels "
	print selstr.format( "labels", ''.join( optstr.format(v, "selected" if v == labels else "" ) for v in tkeys ) )
	print "<br>PC1 "
	print ''.join( rbstr.format( 'pc1', v, "checked" if v == int( pc1 ) else "" ) for v in range( 1, 6 ) )
	print "<br>PC2 "
	print ''.join( rbstr.format( 'pc2', v, "checked" if v == int( pc2 ) else "" ) for v in range( 1, 6 ) )
	print "<br>Reverse Exculuded <input type=checkbox name=revex value=ok %s>" % ( "checked" if revex != "none" else "" )
	print "<br>Point Types "
	print selstr.format( "ptypes", ''.join( optstr.format(v, "selected" if v == ptypes else "" ) for v in [ "samples", "species" ] ) )

print "var datax =  %s;" % json.dumps( coords[:, int( pc1 ) - 1 ].tolist() )# var data = [[5,3], [10,17], [15,4], [2,8]];
print "var datay =  %s;" % json.dumps( coords[:, int( pc2 ) - 1 ].tolist() )# var data = [[5,3], [10,17], [15,4], [2,8]];
if ptypes == "species":
	slabels = [ "" ] * len( kdict )
	clabels = [ "" ] * len( kdict )
	nlabels = [ "" ] * len( kdict )
	for ckey in kdict:
		cind = kdict[ ckey ]
		if labels != "none":
			nlabels[ cind ] = kdnames[ cind ]
		if color != "none":
			clabels[ cind ] = kgnames[ cind ]
		if shape != "none":
			slabels[ cind ] = kgnames[ cind ]
	print "var datac =  %s;" % json.dumps( clabels )# var data = [[5,3], [10,17], [15,4], [2,8]];
	print "var datal =  %s;" % json.dumps( nlabels )# var data = [[5,3], [10,17], [15,4], [2,8]];
	print "var datas =  %s;" % json.dumps( slabels )# var data = [[5,3], [10,17], [15,4], [2,8]];
else:
	print "var datac =  %s;" % json.dumps( mtags[ color ] )# var data = [[5,3], [10,17], [15,4], [2,8]];
	print "var datal =  %s;" % json.dumps( mtags[ labels ] )# var data = [[5,3], [10,17], [15,4], [2,8]];
	print "var datas =  %s;" % json.dumps( mtags[ shape ] )# var data = [[5,3], [10,17], [15,4], [2,8]];
print "var atitle = 'PC ';"
print "var alabels = %s;" % json.dumps( alabels ) 
print "var spdata = [ "
for i in range( len( splist ) ):
	print "{ names: %s, values: %s, signs: %s, pc: \"%5.2f\", ind: %d }" %  ( json.dumps( splist[i] ), json.dumps( svlist[i] ), json.dumps( ssvlist[i] ), pcavar[i], i )
	if i + 1 < len( splist ):
		print ","
print "];"

print """
    var data = [];
    for ( i = 0; i < datax.length; i++ )
    {
		data.push( [ datax[i], datay[i], datac[i], datal[i], datas[i] ] );
	}
    //[5,3], [10,17], [15,4], [2,8]];
    var margin = {top: 20, right: 15, bottom: 20, left: 60}
      , width = 1200 - margin.left - margin.right
      , height = 1100 - margin.top - margin.bottom, bheight = height - 400;
    
    var xscale = d3.scale.linear().range([ 0, width ]);
    
    var yscale = d3.scale.linear().range([ bheight, 0 ]);
	
	var xAxis = d3.svg.axis().scale(xscale).orient('bottom').ticks(2);

	var yAxis = d3.svg.axis().scale(yscale).orient('left').ticks(2);
	
	var xValue = function(d) { return d[0];}; 

	var yValue = function(d) { return d[1];};
	
	var cValue = function(d) { return d[2];};
	var color = d3.scale.category10();
	
	var vmax = Math.max( -d3.min(data, xValue), d3.max(data, xValue), -d3.min(data, yValue), d3.max(data, yValue), 0.6  ) + 0.6;
	xscale.domain([-vmax, vmax]);	
	//xscale.domain([-vmax, vmax]);	
	
	//yscale.domain([d3.min(data, yValue)-0.2, d3.max(data, yValue)+0.2]);	
	yscale.domain([-vmax, vmax]);	
	
	xMap = function(d) { return xscale(xValue(d)); };
 
	yMap = function(d) { return yscale(yValue(d)); };
	
	var sValue = function(d) { return d[4];};
	var slabels = d3.set(data.map(function(d) { return d[4];})).values();	
	var shapeScale = d3.scale.ordinal()
            .domain(slabels)
            .range([ "circle", "cross", "diamond", "square", "triangle-down", "triangle-up" ]);
            
var chart = d3.select('body')
	.append('svg')
	.attr('width', width + margin.right + margin.left)
	.attr('height', height + margin.top + margin.bottom)
	.append("g")
		.attr("transform", "translate(" + margin.left + "," + margin.top + ")");
   

    chart.append('g')
	.attr('transform', 'translate(0,' + bheight / 2 + ')')
	.attr('class', 'x axis')
	.call(xAxis)
	.append("text")
      .attr("class", "label")
      .attr("x", width)
      .attr("y", -6)
      .style("text-anchor", "end")
      .text( atitle + alabels[0] );	


    chart.append('g')
    .attr('transform', 'translate(' + width / 2 + ',0)')
	.attr('class', 'y axis')
	.call(yAxis)
	.append("text")
      .attr("class", "label")
      .attr("transform", "rotate(-90)")
      .attr("y", 6)
      .attr("dy", ".71em")
      .style("text-anchor", "end")
      .text( atitle + alabels[1] );


    //var node = chart.selectAll("g")
      //          .data(data)
      //          .enter()
      //          .append("g");
                
	//node.append("circle")
     // .attr("class", "dot")
     // .attr("r", 5.)
     // .attr("cx", xMap)
     // .attr("cy", yMap)
     // .style("fill", function(d) { return color( cValue( d ) ); } );
    
    //node.append("text")
	//	.attr("x", xMap)
	//	.attr("y", yMap)
	//	.text( function(d) { return d[3]; } );
    
    chart.selectAll(".point")
      .data(data)
      .enter().append("path")
      .attr("class", "point")
      //.attr("r", 5.)
      .attr("transform", function(d) { return "translate(" + xMap(d) + "," + yMap( d ) + ")"; })
      .attr("d", function(d,i){ 
		return d3.svg.symbol().type( shapeScale( d[4] ) )();
		} )
       //.attr("cx", xMap)
      //.attr("cy", yMap)
      .style("fill", function(d) { return color( cValue( d ) ); } )
      //.on("mouseover", function(d) {
      //    tooltip.transition()
      //         .duration(200)
      //         .style("opacity", .9);
      //    tooltip.html(d["Cereal Name"] + "<br/> (" + xValue(d) 
	  //      + ", " + yValue(d) + ")")
      //         .style("left", (d3.event.pageX + 5) + "px")
      //         .style("top", (d3.event.pageY - 28) + "px");
      //})
      //.on("mouseout", function(d) {
      //    tooltip.transition()
      //         .duration(500)
      //         .style("opacity", 0);
      //})
      ;
	chart.selectAll(".dodo")
      .data(data)
      .enter().append("text")
      .attr("class", "dodo")
      .attr("x", xMap)
      .attr("y", yMap)
	  //.attr("dx", ".71em")
      //.attr("dy", ".35em")
	  .text(function(d) { return d[3];});      
	  //.style("fill", function(d) { return color( cValue( d ) ); } )
      //.on("mouseover", function(d) {
      //    tooltip.transition()
      //         .duration(200)
      //         .style("opacity", .9);
      //    tooltip.html(d["Cereal Name"] + "<br/> (" + xValue(d) 
	  //      + ", " + yValue(d) + ")")
      //         .style("left", (d3.event.pageX + 5) + "px")
      //         .style("top", (d3.event.pageY - 28) + "px");
      //})
      //.on("mouseout", function(d) {
      //    tooltip.transition()
      //         .duration(500)
      //         .style("opacity", 0);
      //})
    
	var l2t = 0;
      
"""

if color != "none"	:
	print """
 
   var legend = chart.selectAll(".legend")
      .data( color.domain() )
    .enter().append("g")
      .attr("class", "legend")
      .attr("transform", function(d, i) { return "translate(0," + i * 20 + ")"; });

  legend.append("rect")
      .attr("x", width - 18)
      .attr("width", 18)
      .attr("height", 18)
      .style("fill", color); 
  
  legend.append("text")
      .attr("x", width - 24)
      .attr("y", 9)
      .attr("dy", ".35em")
      .style("text-anchor", "end")
      .text(function(d) { return d;});
  
  l2t = color.domain().length;
      
"""

if shape != "none"	:
	print """
 
   var slegend = chart.selectAll(".slegend")
      .data( shapeScale.domain() )
    .enter().append("g")
      .attr("class", "slegend")
      .attr("transform", function(d, i) { return "translate(0," + ( i + l2t ) * 20 + ")"; });

  slegend.append("path")
	  .attr('stroke', 'black')
      .attr('stroke-width', 1)
	  .attr('transform', 'translate(' + ( width - 9 ) + ',' + 9 + ')')
      .attr("d", function(d,i){ 
		return d3.svg.symbol().type( shapeScale( d ) )();
		 } );
      //.attr("x", width - 18)
      //.attr("width", 18)
      //.attr("height", 18)
      //.style("fill", color); 
  
  slegend.append("text")
      .attr("x", width - 24)
      .attr("y", 9)
      .attr("dy", ".35em")
      .style("text-anchor", "end")
      .text(function(d) { return d;});
      
"""

    
print """
	var c10 = d3.scale.category10();
	for ( i = 0; i < spdata.length; i++ )
	{
		var bclass = ".dodo" + i;
		var tclass = ".coco" + i;
		var cspdata = [];
		var k;
		var dmax = 0;
		for ( k = 0; k < spdata[i].values.length; k++ )
		{
			var cv = +spdata[i].values[k];
			var signv = ( + spdata[i].signs[k] > 0 ) ? 0 : 1;
			cspdata.push( [ k, cv, spdata[i].names[k], signv ] );
			if ( k == 0 ) dmax = cv;
		}
		chart.selectAll( bclass )
			.data( cspdata )
			.enter().append("rect")
				.attr("class", bclass )
				.attr("x", margin.left + i * 200 + 80 )
				.attr("y", function(d) { return bheight + 60 + d[0] * 15; } )
				.attr("width", function(d) { return 80 * d[1] / dmax; } )
				.attr("height", 10 )
				.style("fill", function(d) { if ( d[3] != 0 ) return c10(1); else return c10(0); } );
		chart.selectAll( tclass )
			.data( cspdata )
      .		enter().append("text")
				.attr("class", tclass )
				.attr("x", margin.left + i * 200 + 75 )
				.attr("y", function(d) { return bheight + 70 + d[0] * 15; } )
				.attr("text-anchor", "end")
				.text(function(d) { return d[2].substring(0,20);})
				.style("font-size", "10px");  
	}
	chart.selectAll( ".bobo" )
			.data( spdata )
      .		enter().append("text")
				.attr("class", ".bobo" )
				.attr("x", function(d) { return margin.left + d.ind * 200 + 80 ;} )
				.attr("y", function(d) { return bheight + 50; } )
				//.attr("text-anchor", "end")
				.text(function(d) { return "PC " + ( d.ind + 1 ) + " (" + d.pc + ")";})
				.style("font-size", "12px");  



	d3.select("#generate")
		.on("click", writeDownloadLink);

	function writeDownloadLink(){
		try {
			var isFileSaverSupported = !!new Blob();
		} catch (e) {
			alert("blob not supported");
		}

		var html = '<svg width="960" height="500" title="test2" version="1.1" xmlns="http://www.w3.org/2000/svg">'
					+ d3.select("svg").node().innerHTML + '</svg>';

		var blob = new Blob([html], {type: "image/svg+xml"});
		saveAs(blob, "pca_sp.svg");
	};

    //g.selectAll("scatter-dots")
    //  .data(data)
    //  .enter().append("svg:circle")
    //      .attr("cx", function (d,i) { return x(d[0]); } )
    //      .attr("cy", function (d) { return y(d[1]); } )
    //      .attr("r", 8);
 </script>
 </body>
 </html>
"""