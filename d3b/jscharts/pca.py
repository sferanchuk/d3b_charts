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

import d3bf

form = cgi.FieldStorage()
d3bf.chdir( form.getvalue( "datapath" ) )
dmethod = form.getvalue( "dmethod", "Pearson" )
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
level = form.getvalue( "level" )
fmethod = form.getvalue( "fmethod", "PCA" )
color = form.getvalue( "color", "none" )
shape = form.getvalue( "shape", "none" )
labels = form.getvalue( "labels", "none" )
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

rev = 0
if fmethod == "PCA (sklearn)": #or fmethod == "MDS (sklearn)":
	rev = 1

cdata = d3bf.calc_distances( edata, aedata, dmethod, kdata, knorder, rev )

adist = np.array( cdata )
amax = np.amax( adist )
adist /= amax
alabels = [ "", "" ]

if fmethod.find( "MDS" ) != -1:
	mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
	results = mds.fit(adist)
	coords = results.embedding_
	lcoords = [ coords[:,0].tolist(), coords[:,1].tolist() ]
elif fmethod.find( "PCA" ) != -1:
	pca = decomposition.PCA(n_components=2)
	pca.fit(adist)
	coords = pca.transform( adist )
	alabels = [ "%4.2f" % pca.explained_variance_[0], "%4.2f" % pca.explained_variance_[1] ]
	lcoords = [ coords[:,0].tolist(), coords[:,1].tolist() ]
elif fmethod.find( "PCoA" ) != -1:
	M = pd.DataFrame( np.array( cdata ), site_ids, site_ids )
	ordres = ordination.pcoa( M )
	alabels = [ "%4.2f" % ordres.eigvals[ "PC1" ], "%4.2f" % ordres.eigvals[ "PC2" ] ]
	lcoords = [ ordres.samples[ "PC1" ].tolist(), ordres.samples[ "PC2" ].tolist() ]
elif fmethod == "CA (skbio)":
	M = pd.DataFrame( np.array( cdata ), site_ids, site_ids )
	ordres = ordination.ca( M )
	alabels = [ "%4.2f" % ordres.eigvals[ "CA1" ], "%4.2f" % ordres.eigvals[ "CA2" ] ]
	lcoords = [ ordres.samples[ "CA1" ].tolist(), ordres.samples[ "CA2" ].tolist() ]
	
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
	print "<svg width=\"2400\" height=\"2400\" id=\"normal\"></svg>"
else:
	print "<svg width=\"960\" height=\"960\" id=\"normal\"></svg>"
print "<script>"
print "var datax =  %s;" % json.dumps( lcoords[0] )
print "var datay =  %s;" % json.dumps( lcoords[1] )
print "var datac =  %s;" % json.dumps( mtags[ color ] )
print "var datal =  %s;" % json.dumps( mtags[ labels ] )
print "var datas =  %s;" % json.dumps( mtags[ shape ] )
print "var atitle = '%s';" % ( "dim" if fmethod == "MDS" else "PC" )
print "var alabels = %s;" % json.dumps( alabels ) 

print """
    var pdata = [];
    var i;
    for ( i = 0; i < datax.length; i++ )
    {
		pdata.push( [ datax[i], datay[i], datac[i], datal[i], datas[i] ] );
	}
"""

d3bf.print_drawscatter()

print "drawscatter( \"#normal\", pdata );"
print "</script>"	
 