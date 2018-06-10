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
if len( tagdata ) == 0:
    ( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
    ( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
    print json.dumps( tags )
else:
    newtags = json.loads( newtags )
