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
#import skbio.diversity as skdiv
import pandas as pd
from skbio import DistanceMatrix
import skbio.stats as skbiostats
import scipy.stats as stats
import math

import d3bf

form = cgi.FieldStorage()
id = "emap"
d3bf.chdir( form.getvalue( "datapath" ) )
dgroup = form.getvalue( "dgroup", "none" )
dfilter = form.getvalue( "dfilter", "none" )
level = form.getvalue( "level" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
cunits = form.getvalue( "cunits", "probability" )
pmethod = form.getvalue( "pmethod", "permanova" )


if dgroup == "none" or dgroup == "name":
	sys.exit( 0 )

( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )

maxlevel = ml

dmeasures = [ "Pearson", "Kendall", "Spearman", "Euclidean", "Bray-Curtis", "Jaccard", "Morisita-Horn", "Unifrac-unweighted", "Unifrac-weighted" ]

print "<table class=\"indextable\"><tr><td class=\"columnheader\"><b>Level</b>"
for k in range( len( dmeasures ) ):
	print "<td class=\"columnheader\">" + dmeasures[k]


#for rilevel in range( maxlevel + 1 ):
#	ilevel = maxlevel + 1 - rilevel
if True:
	ilevel = int( level )
	( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, ilevel )
	( findex, mtags ) = d3bf.processtags_m( volumes, tags, dfilter )
	( edata, site_ids, species_ids ) = d3bf.load_edata_m( data, ilevel, mn, ml, kdict, volumes, findex, kdnames )
	aedata = np.array( edata, dtype=float )
	aenorm = np.sum( aedata, axis=1 )
	aedata /= aenorm.reshape( len(edata), 1 )
	
	print "<tr><td class=\"rowheader\">" + str( ilevel )
	for k in range( len( dmeasures ) ):
		print "<td>"
		cdata = d3bf.calc_distances( edata, aedata, dmeasures[k], kdata, knorder, 0 )
		#M = pd.DataFrame( np.array( cdata ), site_ids, site_ids )
		dm = DistanceMatrix( cdata, site_ids )
		if pmethod == "permanova":
			dmres = skbiostats.distance.permanova( dm, mtags[ dgroup ] )
		if pmethod == "anosim":
			dmres = skbiostats.distance.anosim( dm, mtags[ dgroup ] )
		pv = dmres.get( "p-value" ) 
		if cunits == "probability":
			print pv
		if cunits == "log-probability":
			print "%6.2f" % -math.log( float( pv ) )
		
			 
			
print "</table>"

