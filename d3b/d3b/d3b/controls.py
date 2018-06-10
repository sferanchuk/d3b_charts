# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os
import time
import hashlib

dataroot = "pool/"
scriptsroot = "../../jscharts/"

def submit_job( fn, jobname ):
	job = hashlib.md5( str( time.time() ) ).hexdigest()
	os.mkdir( dataroot + job )
	with open( dataroot + job + '/emap.txt', 'wt+') as destination:
		for chunk in f.chunks():
			destination.write(chunk)
	with open( dataroot + job + '/name', 'wt+') as destination:
		destination.write( jobname )
		
	return job

def run_script( params, job, script ):
	os.chdir( dataroot + job )
	os.putenv( "QUERY_STRING", params )
	res = os.popen( "python " + scriptsroot + script + ".py" ).read()
	os.chdir( "../.." )
	return res

			

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