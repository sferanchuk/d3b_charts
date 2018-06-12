#!/usr/bin/python

# Turn on debug mode.
import sys
tfp = '/home/sferanchuk/.local/lib/python2.7/site-packages'
if tfp in sys.path:
	sys.path.remove( tfp )
import cgi
import cgitb
cgitb.enable()
import csv
import numpy as np
#from sklearn import manifold
#from sklearn import decomposition
#import scipy.stats as stats
#from scipy.spatial import distance
import json
import os.path
import collections
import math

import d3bf

form = cgi.FieldStorage()
id = "emap"
dgroup = form.getvalue( "dgroup", "none" )
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
resolution = form.getvalue( "resolution", "low" )
ptype = form.getvalue( "ptype", "proportional-presence" )
samples = form.getlist( "samples" )
ilevel = int( form.getvalue( "level" ) )

( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
( findex, gtags ) = d3bf.processtags( volumes, tags, dfilter, dgroup )
( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, ilevel )
edata = d3bf.load_edata( data, ilevel, ml, kdict, findex, gtags )

ptnum = [ "proportional-presence", "proportional-abundance", "interactive-presence" ].index( ptype )

if ( ptnum == 0 or ptnum == 1 ) and not ( len( samples ) == 2 or len( samples ) == 3 ):
	print "2 or 3 samples should be selected for proportional Venn diagram"
	sys.exit( 0 )

if ( ptnum == 0 or ptnum == 1 ) and ( len( samples ) == 2 or len( samples ) == 3 ):
	if ptnum == 1:
		aedata = np.array( edata, dtype=float )
		aenorm = np.sum( aedata, axis=1 )
		aedata /= aenorm.reshape( len(edata), 1 )
	npair = len( samples ) * ( len( samples ) - 1 ) / 2
	sscores = [ 0. ] * len( samples )
	pscores = [ 0. ] * npair
	psets = []
	jscore = 0.
	jcounts = edata[ gtags[ samples[0] ] ][:]
	jsums = [ 0. ] * len( jcounts )
	vscale = 500
	for i in range( len( samples ) ):
		if ptnum == 0:
			sitot = sorted( edata[ gtags[ samples[i] ] ], reverse=True )
			v = len( sitot )
			if 0 in sitot:
				v = sitot.index( 0 )
		else:
			v = 0
			for vv in aedata[ gtags[ samples[i] ] ].tolist():
				v += vv * vscale
		sscores[ i ] = int( v )
		for j in range( i ):
			v = 0
			for k in range( len( kdict ) ):
				if edata[ gtags[ samples[ i ] ] ][ k ] != 0 and edata[ gtags[ samples[ j ] ] ][ k ] != 0:
					if ptnum == 0:
						v += 1
					else:
						v += 0.5 * ( aedata[ gtags[ samples[ i ] ] ][ k ] + aedata[ gtags[ samples[ j ] ] ][ k ] ) * vscale
			pn = i * ( i - 1 ) / 2 + j
			pscores[ pn ] = int( v )
			psets.append( json.dumps( [ samples[i], samples[j] ] ) )
		for k in range( len( kdict ) ):
			jcounts[k] = min( jcounts[k], edata[ gtags[ samples[i] ] ][k] )
			if ptnum == 1:
				jsums[k] += aedata[ gtags[ samples[i] ] ][k]
	if ptnum == 0:
		sitot = sorted( jcounts, reverse=True )
		v = len( sitot )
		if 0 in sitot:
			v = sitot.index( 0 )
		jscore = v
	else:
		v = 0
		for k in range( len( kdict ) ):
			if jcounts[k] != 0:
				v = jsums[k] / len( samples )		
		jscore = int ( v * vscale )

	if resolution == "high":
		print "<svg width=\"1200\" height=\"1200\" id=\"normal\"></svg>"
	else:
		print "<svg width=\"600\" height=\"600\" id=\"normal\"></svg>" 
	print "<script>"
	print "var vareas = [ "
	for i in range( len( samples ) ):
		print "{ sets: [ \"" + samples[i] + "\" ], size: " + str( sscores[i] ) + "}, "
	for i in range( npair ):
		print "{ sets: " + psets[i]  + ", size: " + str( pscores[i] ) + "}"
		if len( samples ) == 3:
			print ","
	if len( samples ) == 3:
		print "{ sets: " + json.dumps( samples ) + ", size: " + str( jscore ) + "}"
	print "];"
	print """


	//function getSetIntersections() {
	//	areas = d3.selectAll(".venn_area")[0].map(
	//		function (element) { 
	//			return { sets: element.id.split(","), 
	//					size: parseFloat(element.value)};} );
	//	return areas;
	//}

	function getSetIntersections() {
		areas = d3.map( vareas,	function( d ) { return { sets: d[ sets ], size: d[ size ] }; } );
		return areas;
	}

	var svg0 = d3.select( "#normal" ),
		diameter = +svg0.attr("width");

	// draw the initial set
	var chart = venn.VennDiagram()
					.width( diameter - 50 )
					.height( diameter - 50 );

	svg0.datum( vareas ).call(chart);

	// redraw the sets on any change in input
	//d3.selectAll("input").on("change", function() {
	//    d3.select("#dynamic").datum(getSetIntersections()).call(chart);
	//});

	</script>
	"""
	sys.exit( 0 )
elif ptnum == 2 and len( samples ) >= 2 and len( samples ) < 7:
	print "<script>"
	print "$(document).ready(function(){"
	print "$('#jvenn-container').jvenn({"
	print "series: [";
	for i in range( len( samples ) ):
		print "{ name: '" + samples[ i ] + "', \ndata: [";
		for k in range( len( kdict ) ):
			if edata[ gtags[ samples[i] ] ][k] > 0:
				print "'" + kdnames[ knorder[k] ] + "',"
		print "] } "
		if i + 1 < len( samples ):
			print ","
	print """
	],
	fnClickCallback: function() {
				var value = "";
				var vtitle = "";
				if (this.listnames.length == 1) {
					vtitle += "Elements only in ";
				} else {
					vtitle += "Common elements in ";
				}
				for (name in this.listnames) {
					vtitle += this.listnames[name] + " ";
				}
				vtitle += ":";
				for (val in this.list) {
					value += this.list[val] + "\\n";
				}
				$("#vtitle").val( vtitle );
				$("#names").val( value );
			}
		});
	});
	</script>
	<div class="row-fluid">
	<div id="jvenn-container"></div>
	</div>
	<div class="row-fluid">
		<div>
		<p> Click on a venn diagram figure to display the linked elements: </p>
		<input type="text" id="vtitle" size="50" style="background-color:#eee; color: black;" readonly><br><br>
		<textarea readonly id="names" cols="60" style="background-color:#ddd; color: black;" wrap="off" rows="10"></textarea>
	</div>

	"""
	sys.exit(0)
else:
	sys.exit(0)
print """
<?php
echo "</head><body>";
echo "<a href=\"emap.php?id=$id\">back</a><br>\n";


echo "<form action=\"emapjvenn.php\" method=GET><input type=hidden name=id value=\"$id\">\n";
echo "<select name=level>";
for ( $k = 1; $k <= $nh; $k++ ) 
{
	$s = "";
	if ( $k == $level ) $s = "selected";
	echo "<option value=$k $s>$k</option>";
}
echo "</select><br>\n";
$vs0="selected";
$vs1="";
if ( $type == 1 )
{
	$vs0 = "";
	$vs1 = "selected";
}
#echo "<select name=type><option value=0 $vs0>by expression</option><option value=1 $vs1>by presence</option></select><br>\n";
for ( $k = 0; $k < sizeof( $header ) - $nh; $k++ )
{
	$v = $header[ $k + $nh ];
	$cv = "";
	if ( isset( $sets[ $v ] ) ) $cv = "checked";
	echo "$v<input type=checkbox name=sets[$v] value=1 $cv><input type=text name=title[$v] value=\"". $title[$v] ."\"><br>\n";
}
echo "<input type=submit name=\"ok\"></form>\n";
if ( sizeof( $slist ) < 2 )
{
	echo "2 or 3 samples should be selected\n";
	exit();
}
?>
"""
	