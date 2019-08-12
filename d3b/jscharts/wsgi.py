
import os
from cgi import parse_qs, escape


def application(environ, start_response):
	qs = environ['QUERY_STRING']
	d = parse_qs( qs )
	script = d[ 'script' ]
	job = d[ 'job' ]
	base_path = "/home/sergey/work/d3b/pool/"
	scriptsroot =  "/home/sergey/work/d3b/jscharts/"
	pythonpath = "/home/sergey/work/d3b/d3b_py2/bin/python2"
	os.chdir( base_path + job )
	os.putenv( "QUERY_STRING", params )
	response = os.popen( pythonpath + " " + scriptsroot + script + ".py" )
	start_response('200 OK', [('Content-Type', 'text/plain')])
	yield response