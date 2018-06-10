#!/usr/bin/python

# Turn on debug mode.
import sys
import os
import cgi
import cgitb
cgitb.enable()
import csv
import numpy as np
import json
import os.path
import collections
from operator import itemgetter
import hashlib
import time
import d3bf



form = cgi.FieldStorage()
id = "emap"
dgroup = form.getvalue( "dgroup", "none" )
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
dnorm = form.getvalue( "dnorm", "count" )
dorder = form.getvalue( "dorder", "none" )
athreshold = form.getvalue( "athreshold", "none" )
task = form.getvalue( "task", "biom" )
#level = form.getvalue( "level" )
fl = open( "maxlevel" )
maxlevel = int( fl.read() )
fl.close()

##OTU ID T2      T5      T7      L6      OV5     L3      Lm6     L1      L7      OV2     OV7     T3      OV8     L5      L4      T1      v1      v2      v3      T4      T6      Lm1     Lm10   Lm13     Lm15    Lm2     OV4     Lm3     OV6     W       OV1     OV3     Lm8     L2      v3d     Lm7     Lm12    Lm14    Lm9     L8      Lm4     Lm11    Lm16    Lm5     taxonomy
#973124  2.0     1.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0    0.0      0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     k__Bacteria; p__Bacteroidetes; c__[Saprospirae]; o__[Saprospirales]; f__Chitinophagaceae; g__; s__

( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
( findex, gtags ) = d3bf.processtags( volumes, tags, dfilter, dgroup )
ilevel = maxlevel + 1
( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, ilevel )
edata = d3bf.load_edata( data, ilevel, ml, kdict, findex, gtags )
mincount = [ "none", "singletons", "sinletons/doubletons" ].index( athreshold )

def tempname():
	rv = hashlib.md5( str( time.time() ) ).hexdigest()
	return rv[:6]

tprefix = [ "k", "p", "c", "o", "f", "g", "s" ]
fpr = tempname()
tf = open( fpr + ".txt", "w" )
tf.write( "#OTU ID\t" )
gsorted = d3bf.sorted_alnum( gtags.keys() )
for cgtag in gsorted:
	tf.write( cgtag + "\t" )
tf.write( "taxonomy\n" )
for ind in range( len( knorder ) ):
	ckey = knorder[ ind ]
	vs = ""
	zflag = 1
	for cgtag in gsorted:
		vc = edata[ gtags[ cgtag ] ][ kdict[ ckey ] ]
		if vc > mincount:
			vs += "%.1f\t" % float( vc )
			zflag = 0
		else:
			vs += "0.0\t"
	if zflag == 1:
		continue
	tf.write( kdata[ ckey ][ -1 ] + "\t" )
	tf.write( vs )
	for si in range( ilevel - 1 ):
		ac = chr( ord('a' ) + si )
		if ilevel - 1 == 7:
			ac = tprefix[ si ]
		tf.write( ac + "__" + kdata[ ckey ][ si ] )
		if si + 1 < ilevel - 1:
			tf.write( "; " )
	tf.write( "\n" )
tf.close()
os.system( "../../programs/biom convert -i %s.txt -o %s.biom --to-json --table-type=\"OTU table\" --process-obs-metadata taxonomy" % ( fpr, fpr ) )
#os.system( "../../programs/biom convert -i %s.txt -o %s.biom --to-hdf5 --table-type=\"OTU table\" --process-obs-metadata taxonomy" % ( fpr, fpr ) )
if task == "biom":
	bf = open( fpr + ".biom" )
	bdata = bf.read()
	sys.stdout.write( bdata )

if task == "picrust-kegg":
	ppath = "/home/sergey/soft/picrust-1.1.2/"
	os.putenv( "PYTHONPATH", ppath )
	os.system( ppath + "scripts/normalize_by_copy_number.py -i %s.biom -o %s_n.biom" % ( fpr, fpr ) )
	os.system( ppath + "scripts/predict_metagenomes.py -i %s_n.biom -o %s_p.biom" % ( fpr, fpr ) )
	os.system( ppath + "scripts/categorize_by_function.py -i %s_p.biom -c KEGG_Pathways -l 3 -o %s_f.biom" % ( fpr, fpr ) )
	os.system( "../../programs/biom convert -i %s_f.biom -o %s_c.biom --to-json" % ( fpr, fpr ) )
	bf = open( fpr + "_c.biom" )
	bdata = bf.read()
	sys.stdout.write( bdata )


#os.system( "rm %s*.biom %s.txt" % ( fpr, fpr ) )
sys.exit( 0 )

#    normalize_by_copy_number.py -i picrust_input.biom -o normalized_otus.biom
#    predict_metagenomes.py -i normalized_otus.biom -o picrust.biom
