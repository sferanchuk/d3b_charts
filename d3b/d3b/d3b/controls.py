# -*- coding: utf-8 -*-


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

archivepath = "/var/www/html/server/pool/"

def submit_job( reqf, jobname, transform_type ):
	abspath = settings.BASE_DIR + '/'
	job = str( hashlib.md5( str( time.time() ).encode() ).hexdigest() )
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
		with open( abspath + dataroot + job + '/emap_raw.txt', 'wb+') as destination:
			for chunk in f.chunks():
				destination.write(chunk)
		if transform_type == "none":
			os.system( "mv emap_raw.txt emap.txt" )
		else:
			os.system( "python " + abspath + scriptsroot + "transform_emap.py emap_raw.txt emap.txt \"" + transform_type + "\" >transform.out 2> tramsform.log" ) 
	with open( abspath + dataroot + job + '/name', 'wt+') as destination:
		destination.write( jobname )
	os.chdir( abspath )
	return job

def run_script( params, job, script ):
	abspath = settings.BASE_DIR + '/'
	jobpath = abspath + dataroot + job
	if not os.path.isdir( jobpath ):
		jobpath = archivepath + job
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
	#print(res)
	if os.path.isfile( "python.err" ):
		errres = open( "python.err" ).read()
		print( errres )
		
	#with open( "python.err" ) as f:
	#	res += f.read()
	#os.chdir( abspath )
	return res

def render_png( params, job, script, host, jscripts ):
	abspath = settings.BASE_DIR + '/'
	jobpath = abspath + dataroot + job
	if not os.path.isdir( jobpath ):
		jobpath = archivepath + job
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
	print( cmd )
	os.system( cmd )
	pngres = ""
	if os.path.isfile( tempname + ".png" ):
		with open( tempname + ".png", "rb" ) as f:
			pngres = f.read()
		os.remove( tempname + ".png" )
	return pngres

def render_svg( params, job, script, host, jscripts ):
	abspath = settings.BASE_DIR + '/'
	jobpath = abspath + dataroot + job
	if not os.path.isdir( jobpath ):
		jobpath = archivepath + job
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
	print( cmd )
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
			

ss = """

<?php

include "std_include.php";

$name = $_POST[ 'name' ];
$datatype = $_POST['datatype'];
if (!isset( $_POST['input'] ) ) exit();
$input = $_POST['input'];
$job = md5( microtime() );
system( "mkdir pool/$job" );
system( "chmod 777 pool/$job" );
write_name( $job, $name );
$tag = substr( $input, 1, 3 );
$tag1 = substr( $input, 0, 1 );
chdir( "pool/$job" );
if ( $tag == "HDF" || $tag1 == "{" || !is_numeric( $tag1 ) )
{
	$fd = fopen( "input.biom", "w" );
	fwrite( $fd, $input );
	fclose( $fd );
	system( "python ../../aux/biom2emap.py input.biom >emap.txt" );
}
else
{
	$fd = fopen( "emap.txt", "w" );
	$sinp = explode( "\n", $input );
	foreach ( $sinp as $k=>$v )
	{
		fwrite( $fd, rtrim( $v ) . "\t\n" );
	}
	fclose( $fd );
}
exec( "python ../../emap/init_tags.py" );
echo "Job id: ".$job;
?>
"""