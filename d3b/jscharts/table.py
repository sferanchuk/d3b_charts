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
from operator import itemgetter
import d3bf



form = cgi.FieldStorage()
id = "emap"
d3bf.chdir( form.getvalue( "datapath" ) )
dgroup = form.getvalue( "dgroup", "none" )
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
dnorm = form.getvalue( "dnorm", "count" )
dorder = form.getvalue( "dorder", "none" )
level = form.getvalue( "level" )
numbest = form.getvalue( "numbest", "inf" )

( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
( findex, gtags ) = d3bf.processtags( volumes, tags, dfilter, dgroup )
ilevel = int( level ) 
( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, ilevel )
edata = d3bf.load_edata( data, ilevel, ml, kdict, findex, gtags )
	
if True:
		
	aedata = np.array( edata, dtype=float )
	aenorm = np.sum( aedata, axis=1 )
	aedata /= aenorm.reshape( len(edata), 1 )
		
	inumbest = 0
	if numbest.isdigit():
		inumbest = int( numbest )
	if ( inumbest > 0 and inumbest < len( kdict ) ):
		( nedata, nkdict ) = d3bf.select_toptax( edata, kdict, inumbest )
		edata = nedata
		nknorder = []
		naedata = []#[ [] for k in xrange( len( nkdict ) ) ]
		taedata = aedata.transpose()
		for cknum in range( len( knorder ) ):
			ckey = knorder[ cknum ]
			if ckey in nkdict:
				nknorder.append( ckey )
				naedata.append( taedata[ kdict[ ckey ] ].tolist() )
		aedata = np.array( naedata ).transpose()
		kdict = nkdict
		knorder = nknorder

	if dorder == "taxonomy":
		korder = map( itemgetter(1), sorted( kdict.items() ) )
	elif dorder in gtags:
		korder = sorted(range(len(kdict)), key=lambda k: aedata[ gtags[ dorder ] ][k], reverse=True )
	else:
		korder = range( len( kdict ) )
		
	print "<table id=\"restable\" class=\"indextable\"><tr>"
	for si in range( ilevel ):
		print "<td class=\"columnheader\">" + str( si + 1 )
	for cgtag in d3bf.sorted_alnum( gtags.keys() ):
		print "<td class=\"columnheader\">" + cgtag
		
	for ind in korder:
		ckey = knorder[ ind ]
		print "<tr>"
		for si in range( ilevel ):
			print "<td class=\"speciesrow\">" + kdata[ ckey ][ si ]
		for cgtag in d3bf.sorted_alnum( gtags.keys() ):
			if dnorm == "percent":
				print "<td>" + "%.1f" % ( 100. * aedata[ gtags[ cgtag ] ][ ind ] )
			else:
				print "<td>" + str( edata[ gtags[ cgtag ] ][ kdict[ ckey ] ] )
	print "</table>"

