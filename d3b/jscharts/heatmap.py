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
import copy

import d3bf

form = cgi.FieldStorage()
id = "emap"
order1 = form.getvalue( "order1", "none" )
order2 = form.getvalue( "order2", "none" )
labels = form.getvalue( "labels", "none" )
taxlabels = form.getvalue( "taxlabels", "no" )
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
datatype = form.getvalue( "dtype", "count" )
level = form.getvalue( "level" )
numbest = form.getvalue( "numbest", "inf" )
numhighlight = form.getvalue( "numhighlight", "all" )
spshow = form.getvalue( "spshow", "none" )
customtax = form.getvalue( "spcustom", "[]" )
resolution = form.getvalue( "resolution", "low" )
splist = customtax if spshow == "custom" else json.dumps( d3bf.loadfilters( "emap_filters.txt", spshow )[0] )

( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
#( findex, gtags ) = d3bf.processtags( volumes, tags, dfilter, dgroup )
ilevel = int( level ) 
( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, ilevel )
valhighlight = 0. if numhighlight == "all" else 0.01 * float( numhighlight[:-1] )

( findex, mtags ) = d3bf.processtags_m( volumes, tags, dfilter )
( edata, site_ids, species_ids ) = d3bf.load_edata_m( data, ilevel, mn, ml, kdict, volumes, findex, kdnames )

aedata = np.array( edata, dtype=float )
aenorm = np.sum( aedata, axis=1 )
maxdata = np.amax( aedata )
aedata /= aenorm.reshape( len(edata), 1 )
ataxsum = np.sum( aedata, axis=0 )

stags = {}
for j in range( len( edata ) ):
	key = ""
	if order1 != "none":
		key = key + mtags[ order1 ][ j ]
	if order2 != "none":
		key = key + mtags[ order2 ][ j ]
	if labels != "none":
		key = key + mtags[ labels ][ j ]
	key = key + mtags[ "name" ][j]
	stags[ key ] = j
	
ftlist = []
ftindex = []
ostags = collections.OrderedDict( sorted( stags.items() ) )
for k, j in ostags.iteritems():
	ftlist.append( mtags[ labels ][ j ] if labels != "none" else mtags[ "name" ][j] )
	ftindex.append( site_ids[j] )
	

aesnorm = [ 0 ] * len( aedata )
for i in range( len( aenorm ) ):
	ii = site_ids.index( ftindex[i] )
	aesnorm[ i ] = aenorm[ ii ]

nzdict = {}
for i in findex:
	for d in data[1:]:
		ckey = "".join( d[0:ilevel] )
		if ckey in kdict:
			v = float( d[ i + ml + 1 ] )
			if v > 0:
				nzdict[ ckey ] = 1

inumbest = 0
if numbest.isdigit():
	inumbest = int( numbest )
if ( inumbest > 0 and inumbest < len( kdict ) ):
	( nedata, nkdict ) = d3bf.select_toptax( edata, kdict, inumbest )
	edata = nedata
	nknorder = []
	naedata = [ [] for k in xrange( len( nkdict ) ) ]
	taedata = aedata.transpose()
	for cknum in range( len( knorder ) ):
		ckey = knorder[ cknum ]
		if ckey in nkdict:
			nknorder.append( ckey )
			naedata[ nkdict[ ckey ] ] = taedata[ kdict[ ckey ] ].tolist()
	aedata = np.array( naedata ).transpose()
	kdict = nkdict
	knorder = nknorder

		
if datatype == "z-score":
	zedata = []
	rsum = [ 0. ] * len( aedata[0] )
	rsumsq = [ 0. ] * len( aedata[0] )
	rcnt = len( aedata )
	for i in range( len( aedata ) ):
		for k in range( len( aedata[i] ) ):
			v = float( aedata[i][k] )
			rsum[k] += v
			rsumsq[k] += v * v
	rmean = []
	rmdisp = []
	for k in range( len( rsum ) ):
		crmean = rsum[k] / rcnt
		crmdisp = math.sqrt( rsumsq[k] / rcnt - crmean * crmean )
		rmean.append( crmean )
		rmdisp.append( crmdisp )
	for i in range( len( aedata ) ):
		zrow = []
		for k in range( len( aedata[i] ) ):
			v = 0.
			if rmdisp[k] > 0:
				v = ( float( aedata[i][k] ) - rmean[k] ) / rmdisp[k]
			zrow.append( v )
		zedata.append( zrow )
			

treedict = {}
nodecnt = 0
knindex = {}
kdused = []
kdabund = {}
for d in data[1:]:
	cs = "".join( d[0:ilevel] )
	if not cs in kdnames:
		continue
	if not cs in nzdict:
		continue
	if not cs in kdict:
		continue
	if cs in kdused:
		continue
	cd = treedict
	for i in range( ilevel ):
		val = d[i]
		if not val in cd:
			cd[ val ] = { "node": "n" + str( nodecnt ) }
			nodecnt = nodecnt + 1
		cd = cd[val]
	cd[ "leaf" ] = { "title": kdnames[ cs ], "node": "n" + str( nodecnt ) }
	kdabund[ "n" + str( nodecnt ) ] = ataxsum[ kdict[ cs ] ] 
	knindex[ "n" + str( nodecnt + 1 ) ] = cs
	kdused.append( cs )
	nodecnt = nodecnt + 1
	cd[ "node" ] = "n" + str( nodecnt )
	nodecnt = nodecnt + 1
	
def taxabund( d ):
	rv = 0.
	if isinstance( d, dict ):
		if "leaf" in d:
			return kdabund[ d[ "leaf" ][ "node" ] ]
		for k, v in d.iteritems():
			if isinstance( v, dict ):
				if not "leaf" in v:
					rv += taxabund( v )
				else:
					rv += kdabund[ v[ "leaf" ][ "node" ] ]
	return rv

def printtax( d, parent ):
	nn = "n0"
	nname = "node"
	if nname in d:
		nn = d[ nname ]
	ncnt = 0
	nlen = len( d ) - 1 if ( parent != "root" and parent != "null" ) else len( d )
	sklist = []
	for k, v in d.iteritems():
		sklist.append( ( k, v ) )
	sklist.sort( reverse = True, key = lambda k: taxabund( k[1] ) )
	for k, v in sklist:
		if isinstance( v, dict ):
			kn = k if len( k ) > 0 else "n/a"
			abund = taxabund( v )
			if not "leaf" in v:
				print "{ name: \"" + v[ nname ] + "\", title: \"" + kn + "\", parent :\"" + parent + "\", abund : \"" + str( abund ) + "\",  children: [ "
				printtax( v, v[ nname ] )
				print "] }"
			else:
				print "{ name: \"" + v[ nname ] + "\", title: \"" + v[ "leaf" ][ "title" ] + "\", parent :\"" + parent +  "\", abund : \"" + str( abund ) + "\" }"
			if ncnt + 1 < nlen:
				print ","
			ncnt = ncnt + 1
				
maxidsize = 0
for cid in kdnames.values():
	if len( cid ) > maxidsize:
		maxidsize = len( cid )
smallrect = []

print """
<style>

.main text {
    font: 10px sans-serif;	
}

.axis line, .axis path {
    shape-rendering: crispEdges;
    stroke: black;
    fill: none;
}

.node circle {
	  fill: #899;
	  stroke: #344;
	  stroke-width: 1px;
	}
</style>
"""

print "<div id=\"tt\" class=\"tooltip\" style=\"opacity:0; position: absolute; text-align: center; width: 100px; height: 16px; padding: 2px; font: 12px sans-serif; background: #f4f2ec; border: 1px solid #7E713D; border-radius: 8px; pointer-events: none;\"></div>"
if resolution == "high":
	print "<svg width=\"2400\" height=\"%d\" id=\"normal\"></svg>" % ( len( kdict ) * 40 + 600 )
else:
	print "<svg width=\"1200\" height=\"%d\" id=\"normal\"></svg>" % ( len( kdict ) * 20 + 300 )
print "<script>"
print "var ilevel = " + str( ilevel ) + ";"
print "var maxdata = " + str( int( maxdata ) ) + ";"
print "var nsamples = " + str( len( edata) ) + ";"
print "var maxidsize = " + str( maxidsize ) + ";"
print "var taxlabels = \"" + taxlabels + "\";"
print "var datal =  %s;" % json.dumps( ftlist )
print "var datatype =  \"%s\";" % datatype
print "var htnorm =  %s;" % json.dumps( aesnorm )
if resolution == "high":
	print "var minrectsize = 30;"
else:
	print "var minrectsize = 15;"

if len( treedict ) > 1:
	print "var treeData = [ { name: \"root\", title: \"\", parent: \"null\", children: [ "
	printtax( treedict, "root" )
	print " ] } ];"
	print "var multipleroots = 1;"
else:
	print "var treeData = [ "
	printtax( treedict, "null" )
	print " ];"
	print "var multipleroots = 0;"
print "var hmData = { "
ikcnt = 0
for ikey in knindex.keys():
	
	inum = kdict[ knindex[ ikey ] ]
	olst = []
	atsum = 0
	for i in range( len( edata ) ):
		ii = site_ids.index( ftindex[i] )
		if datatype == "z-score":
			v = zedata[ ii ][ inum ]
		elif datatype == "percent":
			v = aedata[ ii ][ inum ] * 100.
		else:
			v = edata[ ii ][ inum ]
		olst.append( edata[ ii ][ inum ] )
		atsum += aedata[ ii ][ inum ]
	print "\"" + ikey + "\": " + json.dumps( olst )
	if atsum / len( edata ) < valhighlight:
		smallrect.append( ikey )
	if ikcnt + 1 < len( knindex ):
		print ","
	ikcnt = ikcnt + 1
print "};"

if datatype == "z-score":
	print "var cdomain = [ -10., -1.2, -0.5, 0., 0.9, 2.2, 10. ]; var cdscale = 1;"
elif datatype == "percent":
	print "var cdomain = [ 0., 0.01, 0.25, 1., 5., 25., 100. ]; var cdscale = 1;"
else:
	print "var cdomain = [ 0, 1, 3, 10, 50, 200, 50000 ]; var cdscale = 0;"
	
print "var smallrect = %s;" % json.dumps( smallrect )
print "var splist = %s;" % splist

print """
// ************** Generate the tree diagram	 *****************
//var margin = {top: 100, right: 15, bottom: 100, left: 60}
//      , width = 1200 - margin.left - margin.right;
//     , height = 500 - margin.top - margin.bottom;
//var margin = {top: 20, right: 120, bottom: 20, left: 120},
//	width = 960 - margin.right - margin.left,
//	height = 500 - margin.top - margin.bottom;

var i = 0;

var val_id = "#normal";
var svg0 = d3.select( val_id ),
    //.attr('id', "mainsvg" ),
    diameter = +svg0.attr("width");
if ( diameter > 1500 )
{
	svg0.append("rect")
		.attr("width", "100%")
		.attr("height", "100%")
		.attr("fill", d3.hsl( 0, 0, 0.95 ) );
}
var fs = diameter / 100;
var margin = { top: 10 * fs, right: 1.5 * fs, bottom: 10 * fs, left: 10 * fs }
    
var width = diameter - margin.left - margin.right;
var height = +svg0.attr( "height" ) - margin.top - margin.bottom;

var ctreewidth = Math.min( ilevel * fs * 2, width / 3 );
var treewidth = ctreewidth;
var maxlevel = ilevel - 1;
var rootshift = 0;
if ( multipleroots == 1 ) 
{
	rootshift = ctreewidth / ( 2 * ilevel );
	treewidth -= rootshift;
	maxlevel = ilevel;
}

var tree = d3.layout.tree()
	.size([height, ctreewidth ])
	.separation(function(a, b) {
			return 1;
           });
           
var diagonal = d3.svg.diagonal()
	.projection(function(d) { return [d.y, d.x]; });


var svg = svg0.append("g")
	.attr("transform", "translate(" + margin.left + "," + margin.top + ")");

root = treeData[ treeData.length - 1 ];
  
function clickcircle( d ) 
{
	if (d.title) {
		//console.log( " mm " + d.title );
		d._title = d.title;
		d.title = null;
		return "";
	} else {
		//console.log( " nn " + d._title );
		d.title = d._title;
		d._title = null;
		return d.title;
	}
}

function update(source) {

  // Compute the new tree layout.
  var nodes = tree.nodes(root).reverse(),
	  links = tree.links(nodes);
  
  var nindex = {};
  var ndata = [];
  var taxroot = 0;
   
  
  // Normalize for fixed-depth.
  nodes.forEach(function(d) { 
	if ( d.depth == 0 ) 
	{ 
		d.y = 0 
	} else 
	{ 
		d.y = d.depth * ctreewidth / maxlevel - rootshift; 
	} 
	if ( d.depth == maxlevel )
	{
		nindex[ d.name ] = d.x;
		if ( d.name in hmData )
		{
			var crow = hmData[ d.name ];
			var ii;
			var mean;
			var mdisp;
			var sfactor = ( smallrect.indexOf( d.name ) != -1 ) ? 0.7 : 1.;
			if ( datatype == "z-score" )
			{
				var sum = 0.;
				var sumsq = 0.;
				for ( ii = 0; ii < crow.length; ii++ )
				{
					var v = crow[ii] / htnorm[ii];
					sum += v;
					sumsq += v * v;
				}
				mean = sum / crow.length;
				mdisp = Math.sqrt( sumsq / crow.length - mean * mean );
			}
			for ( ii = 0; ii < crow.length; ii++ )
			{
				var v = crow[ii];
				if ( datatype == "percent" ) v = ( 100 * crow[ii] ) / htnorm[ii];
				else if ( datatype == "z-score" ) v = ( ( crow[ii] / htnorm[ii] )  - mean ) / mdisp;
				ndata.push( [ d.x, ii, v, "", crow[ii], sfactor ] );
			}
		}
		else
		{
			ndata.push( [ d.x, 0, "ff", d.name, crow[ii], 1 ] );
		}

	}
	if ( d.title == "Bacteria" ) taxroot = d.depth;
  });
  
  var chue = 0.3;
  var csat = 0.3;
  var labelSize = 0;

  // Declare the nodes
  var node = svg.selectAll("g.node")
	  .data(nodes, function(d) { return d.id || (d.id = ++i); });

  // Enter the nodes.
  var nodeEnter = node.enter().append("g")
	  //.attr("class", function(d) { return "node" + (d.children ? " node--internal" : " node--leaf"); })
	  .attr("transform", function(d) { 
		  return "translate(" + d.y + "," + d.x + ")"; });

 nodeEnter.append("circle")
	  .attr("r", fs * 0.3 )
	  .style("fill", d3.hsl( chue, 0.45, 0.3 ) )
		.on("click", function()
		{
			d3.select(this.parentNode ).select( "text" ).text(function(d) { 
			  if (d.title) {
					//console.log( " mm " + d.title );
					if ( typeof removecustomtax == "function" ) removecustomtax( d.title );
					d._title = d.title;
					d.title = null;
					return "";
				} else {
					//console.log( " nn " + d._title );
					d.title = d._title;
					d._title = null;
					if ( typeof addcustomtax == "function" ) addcustomtax( d.title );
					return d.title;
				}
			} );
		} );
		

  nodeEnter.append("text")
	  //.attr("x", function(d) { 
	//	  return d.children || d._children ? -13 : 13; })
	  .attr("dy", function(d) { return d.children ? fs : 0.5 * fs; } )
	  .attr("text-anchor", function(d) { return d.children ? "end" : "start"; } )
	  .attr( "dx", function(d) { return d.children ? -0.2 * fs : 0.5 * fs; } ) 
      .style("font", fs * 1.2 + "px sans-serif" )
	  .text(function(d) { 
		if ( d.children )
		{
			if ( splist.indexOf( d.title ) != -1 )
			{
				if ( typeof addcustomtax == "function" ) addcustomtax( d.title );
				return d.title;
			}
			else
			{
				d._title = d.title;
				d.title = null;
				return ""; 
			}
		} else return d.title;
		}
		)
	  .style("fill-opacity", 1)
	  .each(function(d,i) {
        var thisWidth = this.getComputedTextLength();
        if ( labelSize < thisWidth ) labelSize = thisWidth;
		});
	  
	   var nodeUpdate = node.transition()
      .duration(750);
      //attr("transform", function(d) { return "translate(" + d.y + "," + d.x + ")"; });

	nodeUpdate.select("text")
		.style("fill-opacity", 1);


	  
	  

  // Declare the links
  var link = svg.selectAll("path.link")
	  .data(links, function(d) { return d.target.id; });

  // Enter the links.
  link.enter().insert("path", "g")
	  .attr("class", "link")
	  .attr("d", diagonal)
	.attr("stroke", "black")
    .attr("stroke-width", fs * 0.17 +"px")
    .attr("shape-rendering", "auto")
    .attr("fill", "none");   

	var recth = minrectsize;
	var hmbeg = treewidth + labelSize + 10;
	var rectw = ( width - hmbeg ) / nsamples;
	if ( rectw > minrectsize ) rectw = minrectsize;
	
	var cgenerator = d3.scale.linear()
		.domain([0, 3, cdomain.length-1 ])
		.range([
			d3.hsl( 0, 0, 1),
			d3.hsl( 0, 0, 0.8),
			d3.hsl( 0, 0, 0)]
		);
		//.interpolate(d3.interpolateCubehelix);
		
	var crange = d3.range( cdomain.length ).map(cgenerator);

	var chue = 0.3;
	var csat = 0.3;
	var colorScale = d3.scale.threshold()
        .domain( cdomain )
		.range([
			d3.hsl( chue, csat, 1),
			d3.hsl( chue, csat, 1),
			d3.hsl( chue, csat, 0.9),
			d3.hsl( chue, csat, 0.8),
			d3.hsl( chue, csat, 0.6),
			d3.hsl( chue, csat, 0.3),
			d3.hsl( chue, csat, 0.1)]
		);

   var div = d3.select("#tt");
 	
	svg.selectAll(".bobo")
      .data(ndata)
      .enter().append("rect")
      .attr("class", "bobo")
      .attr("x", function(d) { return d[1] * rectw + hmbeg + 0.16 * rectw + ( 0.5 - 0.5 * d[5] ) * 0.96 * rectw; })
      .attr("y", function(d) { return d[0] - recth * 0.5 + ( 0.5 - 0.5 * d[5] ) * recth;})
	  .attr("width", function(d) { return rectw * 0.96 * d[5]; } )
      .attr("height", function(d) { return recth * d[5]; } )
      .style("stroke", "black" )
      .style("stroke-width", 0.02 * rectw )
      .style("fill", function(d) { return colorScale( d[2] ); } )
      
		.on("mouseover", function(d) {		
			div.transition()		
				.duration(200)		
				.style("opacity", .9);		
			div	.html( d[4] + " [" + ( 100 * d[4] / htnorm[ d[1] ] ).toFixed(1) + "%]" )	
				.style("left", (d3.event.pageX) + "px")		
				.style("top", (d3.event.pageY - 28) + "px");	
			})					
		.on("mouseout", function(d) {		
			div.transition()		
				.duration(500)		
				.style("opacity", 0);	
		});
	  
	var vldata = [];
	for ( i = 0; i < datal.length; i++ )
	{
		vldata.push( [ i, datal[i] ] );
	}
	svg.selectAll(".dodo")
      .data( vldata )
      .enter().append("text")
      .attr("class", "dodo")
      .attr("y", function(d) { return d[0] * rectw + hmbeg + rectw - 1; })
      .attr("x", 10 )
      .attr("transform", "rotate(-90)")
      .style("font", fs * 1.2 + "px sans-serif" )
	  .text(function(d) { return d[1];});      

	if ( taxlabels == "yes" )
	{
		var vtdata = [ "Phylum", "Class", "Order", "Family", "Genus", "Species", "Otu" ];
		var vtaltdata = [ "Phylum", "Class", "Order", "Family", "Genus", "Species", "", "", "Reference", "", "Otu" ];
		var cvtdata = vtdata.slice( 0, ilevel - 1 );
		if ( ilevel > vtdata.length ) cvtdata = vtaltdata.slice( 0, ilevel - 1 );
		svg.selectAll(".toto")
		.data( cvtdata )
		.enter().append("text")
		.attr("class", "toto")
		.attr("y", function( d, i ) { return ( i + 1 + taxroot ) * ctreewidth / maxlevel - rootshift + 0.3 * fs; })
		.attr("x", 0.5 * fs )
		.attr("transform", "rotate(-90)")
		.style("font", fs * 1.2 + "px sans-serif" )
		.style("font-style", "oblique" )
		.text(function(d) { return d;});      
	}

	var vcdata = [];
	var vcrect = 4 * fs;
	var vcheight = 1.4 * fs;
	for ( i = 0; i < cdomain.length - 1; i++ )
	{
		var md = cdomain[ i + 1 ];
		if ( i + 2 == cdomain.length && cdscale == 0 )
		{
			md = maxdata + 1;
		}
		vcdata.push( [ i - 1, cdomain[i], md ] );
	}
	svg.selectAll(".coco1")
      .data( vcdata )
      .enter().append("text")
      .attr("class", "coco1")
      .attr("x", function(d) { return d[0] * vcrect; })
      .attr("y", height + vcheight * 2 )
      .style("font", fs * 1.2 + "px sans-serif" )
	  .text(function(d) { return d[1];});      

	svg.selectAll(".coco2")
      .data( vcdata )
      .enter().append("text")
    .attr("class", "coco2")
      .attr("x", function(d) { return d[0] * vcrect; })
      .attr("y", height + 4.3 * vcheight )
      .style("font", fs * 1.2 + "px sans-serif" )
	  .text(function(d) { return ( cdscale == 0 ) ? ( d[2] - 1 ) : d[2]; });      
	
    
	svg.selectAll(".coco3")
      .data( vcdata )
      .enter().append("rect")
      .attr("class", "coco3")
      .attr("x", function(d) { return d[0] * vcrect; })
      .attr("y", height + vcheight * 2.3 )
	  .attr("width", vcrect )
      .attr("height", vcheight )
	  //.attr("dx", ".71em")
      //.attr("dy", ".35em")
      .style("stroke", "black")
      .style("fill", function(d) { return colorScale( d[1] ); } );


}


function click(d) {
  console.log( "bb " );
  if (d.title) {
	console.log( " mm " + d.title );
	d._title = d.title;
    d.title = null;
    if ( typeof removecustomtax == "function" ) removecustomtax( d.title );
  } else {
	console.log( " nn " + d._title );
    d.title = d._title;
    d._title = null;
    if ( typeof addcustomtax == "function" ) addcustomtax( d.title );
  }
  update( d );
}


update(root);

 </script>
"""