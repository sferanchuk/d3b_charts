#!/usr/bin/python

# Turn on debug mode.
import sys
import cgi
import cgitb
cgitb.enable()
import csv
import numpy as np
from sklearn import manifold
from sklearn import decomposition
import scipy.stats as stats
from scipy.spatial import distance
import json
import os.path
import sqlite3
import time
import math



#form = cgi.FieldStorage()
id = "emap"
reader = csv.reader( open( "emap.txt", "r"), delimiter='\t' )
data = list(reader)
ml = 0
mn = 0
volumes = []
for i in range( 0, len( data[0] ) ):
	if len( data[0][i] ) > 0:
		mn = i
		c0 = data[0][i][0]
		if c0 >= '0' and c0 <= '9':
			ml = i
		else:
			volumes.append( data[0][i] )

fl = open( "maxlevel", "w" )
fl.write( "%d" % ml )
fl.close()
ft = open( "emap_tags.txt", "w" )
ft.write( "*\tnone\tname\t\n" )
for v in volumes:
	ft.write( "%s\t\t%s\t\n" % ( v, v ) )
ft.close()
