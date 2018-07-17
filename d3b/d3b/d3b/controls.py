# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os
import time
import hashlib
from django.conf import settings
import requests

dataroot = "pool/"
scriptsroot = "jscharts/"
phantomjs = "phantomjs/"
staticjs = "d3b/static/"
env = "d3b_py2/"
tempdir = "tmp/"

def submit_job( reqf, jobname ):
	abspath = settings.BASE_DIR + '/'
	job = hashlib.md5( str( time.time() ) ).hexdigest()
	os.mkdir( abspath + dataroot + job )
	os.chdir( abspath + dataroot + job )
	f = reqf[ 'file' ]
	if len( f.name ) > 4 and f.name[ -4: ].lower() == "biom":
		with open( abspath + dataroot + job + '/emap.biom', 'wb+') as destination:
			for chunk in f.chunks():
				destination.write(chunk)
		cmd = "python " + abspath + scriptsroot + "biom2emap.py emap.biom >emap.txt 2>convert.err"
		os.system( cmd )
	else:
		with open( abspath + dataroot + job + '/emap.txt', 'wt+') as destination:
			for chunk in f.chunks():
				destination.write(chunk)
	with open( abspath + dataroot + job + '/name', 'wt+') as destination:
		destination.write( jobname )
	os.chdir( abspath )
	return job

def run_script( params, job, script ):
	abspath = settings.BASE_DIR + '/'
	jobpath = abspath + dataroot + job
	nparams = params
	if len( nparams ) > 0:
		nparams += "&" 
	nparams += "&datapath=" + jobpath
	os.putenv( "QUERY_STRING", nparams )
	if False:
		cmd = abspath + scriptsroot + "runscript.py " + jobpath + " " + script + " regular"
		res = os.popen( cmd ).read()
	cmd = "python " + abspath + scriptsroot + script + ".py" + " 2>python.err"
	res = os.popen( cmd ).read()
	return res

def render_png( params, job, script, host, jscripts ):
	abspath = settings.BASE_DIR + '/'
	jobpath = abspath + dataroot + job
	os.putenv( "QUERY_STRING", params + "&datapath=" + jobpath + "&resolution=high" )
	res = os.popen( "python " + abspath + scriptsroot + script + ".py" ).read()
	tempname = abspath + tempdir + job
	fo = open( tempname + ".html", "w" )
	fo.write( "<html><head>\n" )
	host = "file://" + abspath + staticjs
	for jscript in jscripts:
		fo.write( "<script src=\"%s\"></script>\n" % ( host + jscript ) )
	fo.write( "</head><body style=\"background-color:#f2f2f2;\">\n" )
	fo.write( res )
	fo.write( "</body></html>" )
	fo.close()
	cmd = abspath + phantomjs + "phantomjs " + abspath + phantomjs + "d3b_render.js %s.html %s.png" % ( tempname, tempname )
	print cmd
	os.system( cmd )
	pngres = ""
	if os.path.isfile( tempname + ".png" ):
		with open( tempname + ".png" ) as f:
			pngres = f.read()
		os.remove( tempname + ".png" )
	return pngres

def render_svg( params, job, script, host, jscripts ):
	abspath = settings.BASE_DIR + '/'
	jobpath = abspath + dataroot + job
	os.putenv( "QUERY_STRING", params + "&datapath=" + jobpath + "&resolution=high" )
	res = os.popen( "python " + abspath + scriptsroot + script + ".py" ).read()
	tempname = abspath + tempdir + job
	fo = open( tempname +".html", "w" )
	fo.write( "<html><head>\n" )
	host = "file://" + abspath + staticjs
	for jscript in jscripts:
		fo.write( "<script src=\"%s\"></script>\n" % ( host + jscript ) )
	fo.write( "</head><body>\n" )
	fo.write( res )
	fo.write( "</body></html>" )
	fo.close()
	cmd = abspath + phantomjs + "phantomjs " + abspath + phantomjs + "d3b_savepage.js %s.html %s_rendered.html" % ( tempname, tempname )
	print cmd
	os.system( cmd )
	svgres = ""
	if os.path.isfile( tempname + "_rendered.html" ):
		with open( tempname + "_rendered.html" ) as f:
			htmlres = f.read()
			svgbeg = htmlres.find( "<svg" )
			svgend = htmlres.rfind( "</svg>" )
			if svgbeg != -1 and svgend != -1:
				svgres = htmlres[ svgbeg : svgend + 6 ]
	return svgres
			

