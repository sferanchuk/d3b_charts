#!/usr/bin/python

import os
import sys

path = sys.argv[1]
script = sys.argv[2]
task = sys.argv[3]
location = os.path.abspath(__file__)
scriptsroot = location[ : -len( "runscript.py" ) ]
pythonpath = scriptsroot[ : -len( "jscharts/" ) ] + "d3b_py2/bin/python2"

os.chdir( path )
os.system( pythonpath + " " + scriptsroot + script + ".py 2>python.err" )
