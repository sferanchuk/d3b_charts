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
d3bf.chdir( form.getvalue( "datapath" ) )
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
level = form.getvalue( "level" )
fmethod = form.getvalue( "fmethod", "PCA" )
color = form.getvalue( "color", "none" )
shape = form.getvalue( "shape", "none" )
labels = form.getvalue( "labels", "none" )
dsizes = form.getvalue( "dsizes", "equal" )
spsizes = form.getvalue( "spsizes", "equal" )
pc1 = form.getvalue( "pc1", "1" )
pc2 = form.getvalue( "pc2", "2" )
axis = form.getvalue( "axis", "default" )
resolution = form.getvalue( "resolution", "low" )

ilevel = int( level ) 
( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, ilevel )
( findex, mtags ) = d3bf.processtags_m( volumes, tags, dfilter )
( edata, site_ids, species_ids ) = d3bf.load_edata_m( data, ilevel, mn, ml, kdict, volumes, findex, kdnames )


aedata = np.array( edata, dtype=float )
aenorm = np.sum( aedata, axis=1 )
spenorm = np.sum( aedata, axis=0 )
aedata /= aenorm.reshape( len(edata), 1 )
ncomp = min( 5, len( volumes ) )
ptypes = "abundances"

if dsizes == "equal":
	sizes_data = [ 1. ] * len( aenorm )
else:
	mean_norm = np.mean( aenorm )
	tsdata = ( aenorm / mean_norm )
	if dsizes == "linear":
		sizes_data = tsdata.tolist()
	elif dsizes == "-(1/2)":
		sizes_data = [ math.sqrt( v ) for v in tsdata ]
	elif dsizes == "-(1/4)":
		sizes_data = [ math.sqrt( math.sqrt( v ) ) for v in tsdata ]

if spsizes == "equal":
	sp_sizes_data = [ 1. ] * len( aenorm )
else:
	mean_norm = np.mean( spenorm )
	tsdata = ( spenorm / mean_norm )
	if dsizes == "linear":
		sp_sizes_data = tsdata.tolist()
	elif dsizes == "-(1/2)":
		sp_sizes_data = [ math.sqrt( v ) for v in tsdata ]
	elif dsizes == "-(1/4)":
		sp_sizes_data = [ math.sqrt( math.sqrt( v ) ) for v in tsdata ]
				
#print >> sys.stderr, 'pca_sp 1'		
adist = np.array( edata, dtype=float )
#if ptypes == "species":
#	adist = adist.transpose()
if axis == "default":
	nadist = preprocessing.scale( adist )
elif axis == "1":
	nadist = preprocessing.scale( adist, axis = 1 )
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
sorange = 40
for i in range( ncomp ):
	olist = {}
	for k in range( len( kdict ) ):
		olist[ pow( pcacomp[i][ k ], 2 )] = k
	solist = collections.OrderedDict( sorted( list(olist.items()), reverse=True ) )
	csplist = []
	csglist = []
	csvlist = []
	cssvlist = []
	for ik, ij in solist.items():
		if len( csplist ) < sorange:
			if len( csplist ) < sorange and not ij in sspec:
				sspec.append( ij )
			csplist.append( kdnames[ knorder[ ij] ] )
			csglist.append( kgnames[ knorder[ ij] ] )
			csvlist.append( ik )
			cssvlist.append( pcacomp[i][ ij ] )
	splist.append( csplist )
	sglist.append( csglist )
	svlist.append( csvlist )
	ssvlist.append( cssvlist )
	
spcoords = [ [] for x in range( ncomp ) ]
spnames = []
spgroups = []
for k in sspec:
	spnames.append( kdnames[ knorder[k] ] )
	spgroups.append( kgnames[ knorder[k] ] )
	for i in range( ncomp ):
		spcoords[i].append( pcacomp[ i ][ k ] )

d3bf.print_popupstyle()

if resolution == "high":
	print("<svg width=\"2400\" height=\"3600\" id=\"normal\"></svg>")
else:
	print("<div id=\"tt\" class=\"tooltip\" style=\"opacity:0;\"></div><svg width=\"1200\" height=\"1800\" id=\"normal\"></svg>")

print("<script>")
print("var pc1 = %d; " % int( pc1 ))
print("var pc2 = %d; " % int( pc2 ))
print("var ncomp = %d; " % ncomp)
print("var dataxy = [") 
for i in range( ncomp ):
	print("%s %s" % ( json.dumps( coords[:,i].tolist() ), "," if i + 1 < ncomp else "" ))
print("];");
#print "var datay =  %s;" % json.dumps( coords[:, int( pc2 ) - 1 ].tolist() )
print("var datac =  %s;" % json.dumps( mtags[ color ] ))
print("var datal =  %s;" % json.dumps( mtags[ labels ] ))
print("var datas =  %s;" % json.dumps( mtags[ shape ] ))
print("var sizes_data =  %s;" % json.dumps( sizes_data ))
print("var alabels = %s;" % json.dumps( alabels )) 
print("var atitle = 'PC ';")
print("var tdataxy = [") 
for i in range( ncomp ):
	print("%s %s" % ( json.dumps( spcoords[i] ), "," if i + 1 < ncomp else "" ))
print("];");
print("var tdatal = %s; " % ( json.dumps( spnames ) ))
print("var tdatac = %s; " % ( json.dumps( spgroups ) ))
print("var sizes_tdata = %s; " % ( json.dumps( sp_sizes_data ) ))
#print "var tdatax =  %s;" % json.dumps( tcoords[:, int( pc1 ) - 1 ].tolist() )
#print "var tdatay =  %s;" % json.dumps( tcoords[:, int( pc2 ) - 1 ].tolist() )


print("var spdata = [ ")
for i in range( len( splist ) ):
	print("{ names: %s, values: %s, signs: %s, pc: \"%5.2f\", ind: %d }" %  ( json.dumps( splist[i] ), json.dumps( svlist[i] ), json.dumps( ssvlist[i] ), pcavar[i], i ))
	if i + 1 < len( splist ):
		print(",")
print("];")


print("""
    var pdata = [];
    var i;
    for ( i = 0; i < dataxy[0].length; i++ )
    {
		pdata.push( [ dataxy[ pc1 - 1 ][i], dataxy[ pc2 - 1 ][i], datac[i], datal[i], datas[i], sizes_data[i] ] );
	}
	
	var tpdata = [];
    for ( i = 0; i < tdataxy[0].length; i++ )
    {
		if ( spdata[ pc1 - 1 ].names.indexOf( tdatal[i] ) == -1 && spdata[ pc2 - 1 ].names.indexOf( tdatal[i] ) == -1 ) continue;
		tpdata.push( [ tdataxy[ pc1 - 1 ][i], tdataxy[ pc2 - 1 ][i], tdatac[i], "", "", sizes_tdata[i], tdatal[i] ] );
	}
""")

d3bf.print_drawscatter()

print("""
var chart_all = d3.select( "#normal" ).attr('background-color', 'white' );
var width_all = +chart_all.attr( "width" );
var half_width = width_all / 2;
var ssamples = chart_all.append("g").attr( "transform", "translate( 0, 0 )")
	.attr( "width", half_width ).attr( "height", half_width );
var sspecies = chart_all.append("g").attr( "transform", "translate( " + half_width + ", 0 )" )
	.attr( "width", half_width ).attr( "height", half_width );
var sbars = chart_all.append("g").attr( "transform", "translate( 0, " + half_width + " )" )
	.attr( "width", width_all ).attr( "height", width_all );

drawscatter( ssamples, pdata, "samples" );
drawscatter( sspecies, tpdata, "species" );

function drawbars( chart )
{
	var margin = {top: 20, right: 15, bottom: 60, left: 60};
	var c10 = d3.scale.category10();
	var bheight = 40;
	var width = +chart.attr( "width" ) - margin.right - margin.left;
	var gwidth = width / 5;
	var fs = width / 110.; 
	console.log( "bbb" );
	for ( var i = 0; i < spdata.length; i++ )
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
drawbars( sbars );
</script>
""")
