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
import skbio.diversity as skdiv
import scipy.stats as stats
import math

import d3bf

form = cgi.FieldStorage()
id = "emap"
dgroup = form.getvalue( "dgroup", "none" )
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )


( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
( findex, gtags ) = d3bf.processtags( volumes, tags, dfilter, dgroup )

maxlevel = ml

slist = [ "total", "undef", "gini", "loglog", "log-sqrlog", "gini_index", "ace", "chao1", "shannon", "simpson", "fisher_alpha", "observed_otus", "singles", "doubles" ]
slnames = [ "Total", "Unclassified", "Gini", "PowerLaw", "LogNormal", "SKBio Gini", "Ace", "Chao1", "Shannon", "Simpson", "Fisher Alpha", "Otu Number", "Singletons", "Doubletons" ]

def Gini( distr ):
	tsum = 0
	psum = 0
	n = len( distr )
	if n < 2:
		return 0.
	for i in range( n ):
		tsum += distr[i]
		psum += distr[i] * ( i + 1 )
	return 1. - 2. * ( float( n ) - float( psum ) / tsum ) / ( n - 1 )

def LogLog( distr ):
	datax = []
	datay = []
	xn = float( len( distr ) )
	yn = float( distr[0] )
	for i in range( len( distr ) ):
		vy = math.log( distr[i] / yn )
		vx = math.log( ( i + 1 ) / xn )
		datax.append( vx )
		datay.append( vy )
	(a_s,b_s,r,tt,stderr) = stats.linregress( datax, datay )
	if a_s != a_s:
		return 0
	return a_s

def LogSqrLog( distr ):
	datax = []
	datay = []
	for i in range( len( distr ) ):
		
		vy = math.log( distr[i] )
		vx = math.log( i + 1 ) * math.log( i + 1 )
		datax.append( vx )
		datay.append( vy )
	(a_s,b_s,r,tt,stderr) = stats.linregress( datax, datay )
	if a_s != a_s:
		return 0
	return a_s

def is_undef( s ):
	ukeys = [ "*", "unassigned", "unclassified", "uncultured" ]
	fkeys = [ "gotu_", "fatu_", "ortu_", "cltu_" ]
	if len( s ) == 0 or ( s.lower() in ukeys ):
		return True
	if s.lower().find( "unclassified" ) != -1 :
		return True
	if s[0:5] in fkeys :
		return True
	return False


for rilevel in range( maxlevel + 1 ):
	ilevel = maxlevel + 1 - rilevel
	( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, ilevel )
	undef_list = []
	for ckey in kdata:
		if is_undef( kdata[ ckey ][ -1 ] ):
			undef_list.append( kdict[ ckey ] )

	edata = d3bf.load_edata( data, ilevel, ml, kdict, findex, gtags )
	print "<p><br><br><span class=\"levellabel\">Level %d</span><br>" % ilevel

	print "<table class=\"indextable\"><tr><td>*"
	for cgtag in d3bf.sorted_alnum( gtags.keys() ):
		print "<td class=\"columnheader\">" + cgtag
		
	for k in range( len( slist ) ):
		print "<tr><td class=\"rowheader\">" + slnames[k]
		for cgtag in d3bf.sorted_alnum( gtags.keys() ):
			sitot = sorted( edata[ gtags[ cgtag ] ], reverse=True )
			if 0 in sitot:
				zind = sitot.index( 0 )
				si = sitot[ 0 : zind ]
			else:
				si = sitot
			v = 0
			ef = 0
			try:
				if k == 0:
					v = sum( si )
				elif k == 1:
					v = 0
					for tn in undef_list:
						v += edata[ gtags[ cgtag ] ][ tn ]
				elif k == 2:
					v = -Gini( si )
				elif k == 3:
					v0 = -LogLog( si )
					v = v0
				elif k == 4:
					v = -LogSqrLog( si )
				else:
					v = skdiv.alpha_diversity( slist[ k ], si ).values[0]
			except ValueError as err:
				ef = 1
			except TypeError as err:
				ef = 1
			except ZeroDivisionError as err:
				ef = 1
			if ef == 0:
				sv = "%.2f" % v
				if v == float( int( v ) ):
					sv = "%d" % int( v )
				print "<td>" + sv
			else:
				print "<td>??"
	print "</table>"

