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

import d3bf
import os

form = cgi.FieldStorage()
d3bf.chdir( form.getvalue( "datapath" ) )

if os.path.isfile( "summary.txt" ):
    with open( "summary.txt" ) as f:
	print f.read()
	sys.exit( 0 )

print >>sys.stderr, form.getvalue( "datapath" )

with open( "name" ) as f:
    name=f.read()
id = "emap"
dgroup = form.getvalue( "dgroup", "none" )
dglabels = form.getvalue( "dglabels", "no" )
dtlabels = form.getvalue( "dtlabels", "no" )
#order2 = form.getvalue( "order2", "none" )
#labels = form.getvalue( "labels", "none" )
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
dlegend = form.getvalue( "dlegend", "no" )
level = form.getvalue( "level" )
dnorm = form.getvalue( "dnorm", "percent" )
dprestype = form.getvalue( "dprestype", "bubble" )
resolution = form.getvalue( "resolution", "low" )


( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
( findex, gtags ) = d3bf.processtags( volumes, tags, dfilter, dgroup )
( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, 1 )
taxtype = "none"
if ml == 7 and "Bacteria" in kdnames.values():
    taxtype = "qiime"
rv = json.dumps( { "name": name, "maxlevel": ml + 1, "taxtype": taxtype, "numvolumes": len( volumes ), "volumes": volumes } )
with open( "summary.txt", "w" ) as f:
    f.write( rv )
print rv
