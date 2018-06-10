#!/usr/bin/python

# Turn on debug mode.
import sys
import cgi
import cgitb
cgitb.enable()
import csv
import json
import os.path

import os
import csv

sfname = "emap_filters.txt"
reslist = { "none": [ "0", "" ] }

if os.path.isfile( sfname ):
	sfreader = csv.reader(open( sfname, "r"), delimiter='\t')
	sfdata = list( sfreader )
	for j in range( 0, len( sfdata ) ):
		if len( sfdata[j] ) == 3:
			reslist[ sfdata[j][0] ] = sfdata[j][ 1 : ]
print json.dumps( reslist )
