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

if False:
	form = cgi.FieldStorage()
	dfilter = form.getvalue( "dfilter", "none" )
	level = form.getvalue( "level" )
	fmethod = form.getvalue( "fmethod", "PCA" )
	excludestr = form.getvalue( "exclude", "" )
	excludelist = excludestr.split( "--" )
	excludelist.remove( "" )
	color = form.getvalue( "color", "none" )
	shape = form.getvalue( "shape", "none" )
	labels = form.getvalue( "labels", "none" )
	resolution = form.getvalue( "resolution", "low" )
	pc1 = form.getvalue( "pc1", "1" )
	pc2 = form.getvalue( "pc2", "2" )
	ptype = "species"


( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
ilevel = ml + 1 
excludelist = []
( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, excludelist, ilevel )

treedict = {}
nodecnt = 0
knindex = {}

for d in data[1:]:
	cs = "".join( d[0:ilevel] )
	if not cs in kdnames:
		continue
	cd = treedict
	for i in range( ilevel ):
		val = d[i]
		if not val in cd:
			cd[ val ] = { "node": "n" + str( nodecnt ) }
			nodecnt = nodecnt + 1
		cd = cd[val]
	cd[ "leaf" ] = { "title": kdnames[ cs ], "node": "n" + str( nodecnt ) }
	knindex[ "n" + str( nodecnt + 1 ) ] = cs
	nodecnt = nodecnt + 1
	cd[ "node" ] = "n" + str( nodecnt )
	nodecnt = nodecnt + 1

def printtax( d, parent ):
	nn = "n0"
	nname = "node"
	if nname in d:
		nn = d[ nname ]
	ncnt = 0
	nlen = len( d ) - 1 if ( parent != "root" and parent != "null" ) else len( d )
	for k, v in d.iteritems():
		if isinstance( v, dict ):
			kn = k if len( k ) > 0 else "n/a"
			if not "leaf" in v:
				print "{ name: \"" + v[ nname ] + "\", title: \"" + kn + "\", parent :\"" + parent + "\", children: [ "
				printtax( v, v[ nname ] )
				print "] }"
			else:
				print "{ name: \"" + v[ nname ] + "\", title: \"" + v[ "leaf" ][ "title" ] + "\", parent :\"" + parent + "\" }"
			if ncnt + 1 < nlen:
				print ","
			ncnt = ncnt + 1

if len( treedict ) > 1:
	print "[ { name: \"root\", title: \"\", parent: \"null\", children: [ "
	printtax( treedict, "root" )
	print " ] } ]"
else:
	print "[ "
	printtax( treedict, "null" )
	print " ]"
