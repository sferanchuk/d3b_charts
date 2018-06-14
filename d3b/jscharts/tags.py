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

form = cgi.FieldStorage()
id = "emap"

tagdata = form.getvalue( "newtags", "" )

if len( tagdata ) > 0:
	newtags = json.loads( tagdata )
	volumes = newtags[ 'name' ]
	with open( "emap_tags.txt", "w" ) as f:
		f.write( "name\t" )
		tkeys = newtags.keys()
		for tag in tkeys:
			if tag != "name" and tag != "none":
				f.write( tag + "\t" )
		f.write( "\n" )
		for i in range( len( volumes ) ):
			f.write( volumes[i] + "\t" )
			for tag in tkeys:
				if tag != "name" and tag != "none":
					f.write( newtags[ tag ][ i ] + "\t" )
			f.write( "\n" )

( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
print json.dumps( tags )
