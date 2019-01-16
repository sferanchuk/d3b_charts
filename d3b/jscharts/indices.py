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
import skbio.stats as skstats
from scipy.special import comb
from scipy.optimize import fmin_powell
import math

import d3bf

form = cgi.FieldStorage()
id = "emap"
d3bf.chdir( form.getvalue( "datapath" ) )
dgroup = form.getvalue( "dgroup", "none" )
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
mmethod = form.getvalue( "mmethod", "all" )
level = form.getvalue( "level", "all" )

( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
( findex, gtags ) = d3bf.processtags( volumes, tags, dfilter, dgroup )

maxlevel = ml

slist = [ "total", "undef", "gini", "loglog", "log-sqrlog", "gini_index", "ace", "chao1", "shannon", "simpson", "fisher_alpha", "observed_otus", "singles", "doubles", "mm-fit" ]
slnames = [ "Total", "Unclassified", "Gini", "PowerLaw", "LogNormal", "SKBio Gini", "Ace", "Chao1", "Shannon", "Simpson", "Fisher Alpha", "Otu Number", "Singletons", "Doubletons", "Michaelis-Menten" ]

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
	if level != "all" and ilevel != int( level ):
		continue
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
		
	krange = range( len( slist ) - 1 ) if mmethod != "mm-fit" else [ slist.index( "mm-fit" ) ]
	for k in krange:
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
			estr = ""
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
				elif slist[k] == "mm-fit":
					
					n_indiv = sum( si )
					n_otu = len( si )
					def subsample( si, i ):
						ssi = skstats.subsample_counts( si, i )
						return np.count_nonzero( ssi )

					i_step = max( n_indiv / 200, 1 )
					num_repeats = max( 2000 / i_step, 1 ) 
					print >>sys.stderr, ( i_step, num_repeats )
					S_max_guess = n_otu
					B_guess = int( round( n_otu / 2 ) )
					params_guess = ( S_max_guess, B_guess )
					xvals = np.arange( 1, n_indiv, i_step )
					ymtx = np.empty( ( num_repeats, len( xvals ) ), dtype=int )
					for i in range( num_repeats ):
						ymtx[i] = np.asarray( [ subsample( si, n ) for n in xvals ], dtype=int )
					yvals = ymtx.mean(0)
					def errfn(p, n, y):
						return ( ( ( p[0] * n / (p[1] + n ) ) - y ) ** 2 ).sum()
						#return ( ( p[0] * ( 1. - np.exp( n / p[1] ) ) - y ) ** 2 ).sum()
					params_guess = ( n_otu,  int( round( n_otu / 2 ) ) )
					print >>sys.stderr, yvals
					print >>sys.stderr, xvals

					mparams = fmin_powell( errfn, params_guess, ftol=1e-5, args=(xvals, yvals),	disp = False )
					ef = 2
					sv = "%.2f %.2f %.2f" % ( mparams[0], mparams[1], math.sqrt( errfn( mparams, xvals, yvals) / len( xvals ) ) )
				else:
					v = skdiv.alpha_diversity( slist[ k ], si ).values[0]
			except ( ValueError, TypeError, ZeroDivisionError ) as err:
				estr = str( err )
				ef = 1
			#except TypeError as err:
			#	ef = 1
			#except ZeroDivisionError as err:
			#	ef = 1
			if ef == 0:
				sv = "%.2f" % v
				if v == float( int( v ) ):
					sv = "%d" % int( v )
				print "<td>" + sv
			elif ef == 1:
				print "<td>??"
				print >>sys.stderr, estr
			elif ef == 2:
				print "<td>" + sv
	print "</table>"

