# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os
import time
import hashlib
from django.conf import settings

dataroot = "pool/"
scriptsroot = "jscharts/"
phantomjs = "phantomjs/"
staticjs = "d3b/static/"


def submit_job( fn, jobname ):
	abspath = os.path.abspath( "." ) + "/"
	job = hashlib.md5( str( time.time() ) ).hexdigest()
	os.mkdir( abspath + dataroot + job )
	with open( abspath + dataroot + job + '/emap.txt', 'wt+') as destination:
		for chunk in f.chunks():
			destination.write(chunk)
	with open( dataroot + job + '/name', 'wt+') as destination:
		destination.write( jobname )
		
	return job

def run_script( params, job, script ):
	abspath = settings.BASE_DIR + '/'
	os.chdir( abspath + dataroot + job )
	os.putenv( "QUERY_STRING", params )
	cmd = "python " + abspath + scriptsroot + script + ".py 2>python.err"
	print params
	print cmd
	res = os.popen( cmd ).read()
	with open( "python.err" ) as f:
		res += f.read()
	os.chdir( abspath )
	return res

def render_png( params, job, script, host, jscripts ):
	abspath = settings.BASE_DIR + '/'
	os.chdir( abspath + dataroot + job )
	os.putenv( "QUERY_STRING", params + "&resolution=high" )
	res = os.popen( "python " + abspath + scriptsroot + script + ".py" ).read()
	fo = open( "lastres.html", "w" )
	fo.write( "<html><head>\n" )
	host = "file://" + abspath + staticjs
	for jscript in jscripts:
		fo.write( "<script src=\"%s\"></script>\n" % ( host + jscript ) )
	fo.write( "</head><body>\n" )
	fo.write( res )
	fo.write( "</body></html>" )
	fo.close()
	cmd = abspath + phantomjs + "phantomjs " + abspath + phantomjs + "d3b_render.js lastres.html lastres.png >phantomjs.log 2>phantomjs.err"
	print cmd
	os.system( cmd )
	pngres = ""
	if os.path.isfile( "lastres.png" ):
		with open( "lastres.png" ) as f:
			pngres = f.read()
		os.remove( "lastres.png" )
	os.chdir( abspath )
	return pngres

def render_svg( params, job, script, host, jscripts ):
	abspath = settings.BASE_DIR + '/'
	os.chdir( abspath + dataroot + job )
	os.putenv( "QUERY_STRING", params + "&resolution=high" )
	res = os.popen( "python " + abspath + scriptsroot + script + ".py" ).read()
	fo = open( "lastres.html", "w" )
	fo.write( "<html><head>\n" )
	host = "file://" + abspath + staticjs
	for jscript in jscripts:
		fo.write( "<script src=\"%s\"></script>\n" % ( host + jscript ) )
	fo.write( "</head><body>\n" )
	fo.write( res )
	fo.write( "</body></html>" )
	fo.close()
	cmd = abspath + phantomjs + "phantomjs " + abspath + phantomjs + "d3b_savepage.js lastres.html lastres_rendered.html >phantomjs.log 2>phantomjs.err"
	print cmd
	os.system( cmd )
	svgres = ""
	if os.path.isfile( "lastres_rendered.html" ):
		with open( "lastres_rendered.html" ) as f:
			htmlres = f.read()
			svgbeg = htmlres.find( "<svg" )
			svgend = htmlres.rfind( "</svg>" )
			if svgbeg != -1 and svgend != -1:
				svgres = htmlres[ svgbeg : svgend + 6 ]
	os.chdir( abspath )
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