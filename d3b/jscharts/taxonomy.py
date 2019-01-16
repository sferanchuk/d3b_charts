#!/usr/bin/python

# Turn on debug mode.
import sys
import cgi
import cgitb
cgitb.enable()
import csv
import json
import os.path
import d3bf

import os
import csv

sfname = "emap_filters.txt"
reslist = { "none": [ "0", "" ] }

form = cgi.FieldStorage()
d3bf.chdir( form.getvalue( "datapath" ) )

nsfdata = form.getvalue( "newfilters", "" )

if len( nsfdata ) > 0:
	newfilters = json.loads( nsfdata )
	with open( sfname, "w" ) as f:
		for sfilter in newfilters:
			nlist = ";".join( newfilters[ sfilter ][1].strip().split( "\n" ) )
			f.write( "%s\t%s\t%s\n" % ( sfilter, newfilters[ sfilter ][ 0 ], nlist ) )

if os.path.isfile( sfname ):
	sfreader = csv.reader(open( sfname, "r"), delimiter='\t')
	sfdata = list( sfreader )
	for j in range( 0, len( sfdata ) ):
		if len( sfdata[j] ) == 3:
			reslist[ sfdata[j][0] ] = [ sfdata[j][ 1 ], sfdata[j][2].replace( ";", "\n" ) ]
print json.dumps( reslist )
