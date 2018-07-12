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
d3bf.chdir( form.getvalue( "datapath" ) )
dgroup = form.getvalue( "dgroup", "none" )
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
cmethod = form.getvalue( "cmethod", "anova" )
cunits = form.getvalue( "cunits", "probability" )

( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
#( findex, gtags ) = d3bf.processtags( volumes, tags, dfilter, dgroup )
( findex, mtags ) = d3bf.processtags_m( volumes, tags, dfilter )

maxlevel = ml

gtags = []
for vind in range( len( mtags[ dgroup ] ) ):
	if not mtags[ dgroup ][ vind ] in gtags:
		gtags.append( mtags[ dgroup ][ vind ] )

slist = [ "gini", "loglog", "log-sqrlog", "gini_index", "ace", "chao1", "shannon", "simpson", "fisher_alpha", "observed_otus", "singles", "doubles" ]
slnames = [ "Gini", "PowerLaw", "LogNormal", "SKBio Gini", "Ace", "Chao1", "Shannon", "Simpson", "Fisher Alpha", "Otu Number", "Singletons", "Doubletons" ]

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

print "<table class=\"indextable\"><tr><td class=\"columnheader\"><b>Level</b>"
for k in range( len( slist ) ):
	print "<td class=\"columnheader\">" + slnames[k]


for rilevel in range( maxlevel + 1 ):
	ilevel = maxlevel + 1 - rilevel
	( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, ilevel )
	( edata, site_ids, species_ids ) = d3bf.load_edata_m( data, ilevel, mn, ml, kdict, volumes, findex, kdnames )
	print "<tr><td class=\"rowheader\">" + str( ilevel )
	for k in range( len( slist ) ):
		#print "<tr><td class=\"rowheader\">" + slnames[k]
		rlist = [ [] for i in range( len( gtags ) ) ]
		for vind in range( len( edata ) ):
			sitot = sorted( edata[ vind ], reverse=True )
			if 0 in sitot:
				zind = sitot.index( 0 )
				si = sitot[ 0 : zind ]
			else:
				si = sitot
			v = 0
			ef = 0
			try:
				if k == 0:
					v = -Gini( si )
				elif k == 1:
					v0 = -LogLog( si )
					v = v0
				elif k == 2:
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
				gind = gtags.index( mtags[ dgroup ][ vind ] )
				rlist[ gind ].append( v )
		ef = 0
		pvalue = 1.
		try:
			if cmethod == "anova":
				( fstat, pvalue ) = stats.f_oneway( *rlist )
			elif cmethod == "best-ttest":
				for i in range( len( rlist ) ):
					for j in range( i ):
						( tstat, cpv ) = stats.ttest_ind( rlist[ i ], rlist[ j ] )
						pvalue = min( pvalue, cpv )
			elif cmethod == "best-wilcoxon":
				for i in range( len( rlist ) ):
					for j in range( i ):
						( stat, cpv ) = stats.mannwhitneyu( rlist[ i ], y=rlist[ j ] )
						pvalue = min( pvalue, cpv )
		except ValueError as err:
			ef = 1
		except TypeError as err:
			ef = 1
		except ZeroDivisionError as err:
			ef = 1
		strval = "??"
		if ef == 0 and not math.isnan( pvalue ):
			if cunits == "log-probability":
				if pvalue == 0:
					strval = "Inf"
				else: 
					strval = "%5.2f" % -math.log( pvalue )
			else:
				strval = "%5.3f" % pvalue
		print "<td>" + strval
			
print "</table>"

