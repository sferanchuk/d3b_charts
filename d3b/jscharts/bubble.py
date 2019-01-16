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
d3bf.chdir( form.getvalue( "datapath" ) )
dgroup = form.getvalue( "dgroup", "none" )
dglabels = form.getvalue( "dglabels", "no" )
dtlabels = form.getvalue( "dtlabels", "no" )
#order2 = form.getvalue( "order2", "none" )
#labels = form.getvalue( "labels", "none" )
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
dlegend = form.getvalue( "dlegend", "no" )
level = form.getvalue( "level" )
dnorm = form.getvalue( "dnorm", "percent" )
dprestype = form.getvalue( "dprestype", "bubble" )
resolution = form.getvalue( "resolution", "low" )


( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
( findex, gtags ) = d3bf.processtags( volumes, tags, dfilter, dgroup )
ilevel = int( level ) 

( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, ilevel )

edata = d3bf.load_edata( data, ilevel, ml, kdict, findex, gtags )

aedata = np.array( edata, dtype=float )
aenorm = np.sum( aedata, axis=1 )
aedata /= aenorm.reshape( len(edata), 1 )
if dnorm == "percent-quantile":
	aedata = d3bf.quantileNormalize( aedata )

minsizes = []
for rnum in range( len( edata ) ):
	csrow = sorted( edata[ rnum ], reverse=True )
	rlen = len( csrow )
	if 0 in csrow:
		rlen = csrow.index( 0 )
	if rlen < 70:
		minsize = 0
	else:
		if dprestype == "treemap" and False:
			rind = min( rlen - 1, 70 )
		else:
			rind = min( rlen - 1, max( 70, rlen / 3 ) )
		minsize = csrow[ rind ] 
	minsizes.append( minsize )

print """
<style>
circle {
  fill: rgb(31, 119, 180);
  fill-opacity: .25;
  stroke: rgb(31, 119, 180);
  stroke-width: 1px;
}

.leaf circle {
  fill: #ff7f0e;
  fill-opacity: 1;
}

text {
  font: 10px sans-serif;
  text-anchor: middle;
}

.node {
  font: 10px sans-serif;
  line-height: 12px;
  overflow: hidden;
  position: absolute;
  text-indent: 2px;
}
</style>
"""

if dprestype == "bubble" or dprestype == "bars":
	if resolution == "high":
		print "<svg width=\"2400\" height=\"2400\" id=\"normal\"></svg>"
	else:
		print "<svg width=\"960\" height=\"960\" id=\"normal\"></svg>"
#	print "<svg width=2400 height=2400 id=highres></svg>"
else:
	if resolution == "high":
		print "<div width=\"2400\" height=\"2400\" id=\"normal\"></div>"
	else:
		print "<div width=\"960\" height=\"960\" id=\"normal\"></div>"
#	print "<div width=960 height=500 id=normal></div>"
#	print "<div width=2000 height=1200 id=highres></div>"

print "<script type=\"text/javascript\">"
print "var height = 600;"
print "var outfn = \"bubble_%s_%d_%s\";" % ( id, ilevel, dgroup ) 
print "var tfont = \"12px sans-serif\";"
print "var ilevel = " + str( ilevel ) + ";"
print "var dglabels = \"" + dglabels + "\";"
print "var dtlabels = \"" + dtlabels + "\";"
print "var dlegend = \"" + dlegend + "\";"
#print "var maxdata = " + str( int( aenorm ) ) + ";"
#print "var nsamples = " + str( len( edata) ) + ";"
#print "var maxidsize = " + str( maxidsize ) + ";"
print "var rroot =  { name: \"root\", children: ["
gsize = len( gtags )
gscnt = 0
ksize = len( edata[0] )
ordlist = []
for gnum in sorted( gtags.values() ):
	gkey = gtags.keys()[ gtags.values().index( gnum ) ]
	print "{ sname: \"%s\", children: [ " % gkey
	kcnt = 0
	others_sum = 0.
	for ckey in kdnames:
		csize = edata[ gscnt ][ kdict[ ckey ] ]
		if csize <= minsizes[ gscnt ]:
			others_sum += csize
		else:
			if dnorm == "percent" or dnorm == "percent-quantile":
				csize = int( 10000 * aedata[ gscnt ][ kdict[ ckey ] ] )
			dname = kdnames[ ckey ]
			dratio = aedata[ gscnt ][ kdict[ ckey ] ]
			percent = "%4.1f %%" % ( dratio * 100. )
			#pdig = [ch.isdigit() for ch in dname].index( True )
			print "{ name: \"%s\", sname: \"%s\", order: \"%s\", size: %d, percent: \"%s\" } " % ( kdnames[ ckey ] + "\\n" + percent, kdnames[ ckey ], kgnames[ ckey ], csize, percent  )
			if kcnt + 1 < ksize or others_sum > 0:
				print ","
			if not kgnames[ ckey ] in ordlist and dratio > 0.01:
				ordlist.append( kgnames[ ckey ] )
		kcnt = kcnt + 1
	if others_sum > 0:
		if dnorm == "percent":
			csize = int( 10000 * others_sum / aenorm[ gscnt ] )
		else:
			csize = others_sum
		percent = "%4.1f %%" % ( others_sum * 100. / aenorm[ gscnt ] )
		#pdig = [ch.isdigit() for ch in dname].index( True )
		print "{ name: \"Others\\n%s\", sname: \"Others\", order: \"Mixed\", size: %d, percent: \"%s\" } " % ( percent, csize, percent  )
	print "] }"
	if gscnt + 1 < gsize:
		print ","
	gscnt = gscnt + 1
print "] };"
print "var ordlist = %s;" % json.dumps( ordlist )
print """

//var margin = {top: 100, right: 15, bottom: 100, left: 60};

root = d3.hierarchy(rroot)
      .sum(function(d) { return d.size; })
      .sort(function(a, b) { return b.value - a.value; });

var color = d3.scaleOrdinal(d3.schemeCategory20);

function draw( val_id ) {

"""
if dprestype == "bubble":
	print """

	var svg0 = d3.select( val_id );
	var diameter = +svg0.attr("width") * 0.75;
	var fs = diameter / 100;
	margin = { top: 8 * fs, right: 2 * fs, bottom: 16 * fs, left: 8 * fs };

	var svg = svg0.append("g")
		.attr("transform", "translate(" + margin.left + "," + margin.top + ")");

	var g = svg.append("g").attr("transform", "translate(2,2)"),
	format = d3.format(",d");

	var pack = d3.pack()
		.size([diameter - 4, diameter - 4]);

	var maxbottom = 0;

	var node = g.selectAll(".node")
		.data(pack(root).descendants())
		.enter().append("g")
			.attr("class", function(d) { return d.children ? "node" : "leaf node"; })
			.attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; });

	node.append("title")
		.style("font", tfont )
		.style( "text-anchor", "middle" )
		.text(function(d) { return d.data.name; });

	node.append("circle")
		.attr("r", function(d) { return d.r; })
		.style( "stroke", "rgb(31, 119, 180)" )
		.style( "stroke-width", function(d) 
		{ 
			if ( d.data.name == "root" ) 
			{
				return "0px";
			}
			else if ( d.children )
			{
				return "1px";
			}
			else 
			{
				maxbottom = Math.max( maxbottom, d.y + d.r ); 				
				return "1px"; 
		} } )
		.style( "fill-opacity", function(d) 
		{ 
			if ( d.data.name == "root" ) 
			{
				return "0.0";
			}
			else if ( d.children )
			{
				return "0.0";
			}
			else return "0.5"; 
		})
		.style("fill", function(d) 
		{ 
			if ( d.data.name == "root" ) 
			{
				return d3.hsl( 0.3, 0.3, 0.1 );
			}
			else if ( d.children )
			{
				return d3.hsl( 0.3, 0.3, 0.2 );
			}
			else return color(d.data.order); 
		});
			
	var smap = d3.scaleOrdinal().domain( [ 1, 2, 3, 4, 5 ] ).range( [ 1., 1.2, 1.25, 1.3, 1.33 ] );
	function ssmap( v ) {	if ( v <= 5 ) return smap( v );	return 1000; }
	var maxfs = fs * 0.2 * Math.sqrt( 8. / fs );
	var minfs = maxfs * 0.5;
  
	function scaledsize( radius, textlength, dimparam )
	{
		var curfs = maxfs * Math.sqrt( radius );
		var lenfs = 5 * dimparam * radius / textlength;
		return Math.max( dimparam * curfs, Math.min( curfs, lenfs ) ); 
	}

	node.filter(function(d) { return !d.children; }).append("text")
		.attr( "dy", function(d) { return -0.2 * scaledsize( d.r, 1, 0.7 ) + "px"; } ) 
		.style( "font", function(d) { return scaledsize( d.r, d.data.sname.length, 0.7 ) + "px sans-serif"; } )
		.style( "text-anchor", "middle" )
		.text( function(d) { return d.data.sname; } )
		.each(function(d,i) 
		{
			d.drawtext = true;
			var thisWidth = this.getComputedTextLength();
			if ( thisWidth > 1.95 * d.r )
			{
				this.remove();
				d.drawtext = false;
			}
		}
		);

  
	if ( dglabels == "yes" )
	{
		node.filter(function(d) { return d.children && d.data.name != "root"; }).append("text")
			.attr("dy", function( d ) { return ( d.r + fs * 2 ); } )
			.attr("dx", fs )
			.style("font", fs * 2. + "px sans-serif" )
			.style("font-weigth", "bolder" )
			.style( "text-anchor", "middle" )
			.text(function(d) { return d.data.sname; } 
			);
	}
	if ( dtlabels == "yes" )
	{
		node.filter(function(d) { return !d.children; }).append("text")
			.attr( "dy", function(d) { return 0.8 * scaledsize( d.r, 1, 0.7 ) + "px"; } ) 
			.style( "font", function(d) { return scaledsize( d.r, 1, 0.7 ) + "px sans-serif"; } )
			.style( "font-weigth", "bolder" )
			.style( "text-anchor", "middle" )
			.text(function(d) { return d.drawtext ? d.data.percent : ""; } 
			);
	}

	
	if ( dlegend == "yes" )
	{
		var cellwidth = fs * 11;
		var cellheight = fs * 1.3;
		var height = Math.min( diameter + cellheight * 4, maxbottom + cellheight * 6 );
		svg.selectAll(".coco1")
			.data( ordlist )
			.enter().append("text")
			.attr("class", "coco1")
			.attr("x", function(d,i) { return i * cellwidth; } )
			.attr("y", function(d,i) { return height + 5 * cellheight + ( i % 2 ) * fs * 1.8; } )
			.style("font", 1.6 * fs + "px sans-serif" )
			.style("text-anchor", "start" ) 
			.text(function(d) { return d;}
			);      

		svg.selectAll(".coco3")
			.data( ordlist )
			.enter().append("rect")
			.attr("class", "coco3")
			.attr("x", function(d,i) { return i * cellwidth; })
			.attr("y", height + 2 * cellheight )
			.attr("width", cellwidth )
			.attr("height", cellheight )
			.style("stroke",  "black")
			.style("stroke-width", "0.5px" )
			.style("fill-opacity", 0.5 )
			.style("fill", function(d) { return color( d ); } 
			);
	}

      
"""
elif dprestype == "bars":
	print """

	var svg0 = d3.select( val_id );
	var diameter = +svg0.attr("width") * 0.75;
	var fs = diameter / 100;
	var margin = { top: 8 * fs, right: 2 * fs, bottom: 16 * fs, left: 4 * fs };

	var svg = svg0.append("g")
		.attr("transform", "translate(" + margin.left + "," + margin.top + ")");

	var rectdata = [];
	var labeldata = [];
	var x = margin.left;
	var labelsize = diameter * 0.5;
	var ysize = diameter * 0.8;
	var xbsize = ( diameter * 0.6 ) / rroot.children.length;
	console.log( rroot );
	for ( var i = 0; i < rroot.children.length; i++ )
	{
	
		var sample = rroot.children[ i ];
		var sum = 0;
		for ( var j = 0; j < sample.children.length; j++ )
		{
			sum += sample.children[ j ].size;
		}
		//console.log( sum );
		//console.log( root[ "children" ][ i ][ "children" ] );
		//console.log( sample );
		var y = margin.top;
		for ( var j = 0; j < sample.children.length; j++ )
		{
			var width = sample.children[j].size * ysize / sum; 
			rectdata.push( [ x, width, y, sample.children[ j ].order, sample.children[ j ].percent ] );
			y += width;
		}
		labeldata.push( sample.sname );
		x += xbsize;
	}
	
	console.log( rectdata.length );
	console.log( labeldata.length );
		
	svg.selectAll(".momo1")
		.data( rectdata )
		.enter().append("rect")
		.attr("class", "momo1")
		.attr("x", function(d) { return d[0]; } )
		.attr("y", function(d) { return d[2]; } )
		.attr("height", function(d) { return d[1]; } )
		.attr("width", xbsize * 0.8 )
		.style("stroke",  "black")
		.style("stroke-width", "0.5px" )
		.style("fill-opacity", 0.5 )
		.style("fill", function(d) { return color( d[3] ); } 
		);

	svg.selectAll(".momo1m")
		.data( rectdata )
		.enter().append("text")
		.attr("class", "momo1m")
		.attr("x", function(d) { return d[0]; } )
		.attr("y", function(d) { return d[2] + 1.4 * fs; } )
		.style("font", 1.6 * fs + "px sans-serif" )
		.style("text-anchor", "start" ) 
		.text(function(d) { return d[4];}
		);

	svg.selectAll(".momo2")
		.data( labeldata )
		.enter().append("text")
		.attr("class", "momo2")
		.attr("x", function(d,i) { return xbsize * i + margin.left; } )
		.attr("y", function(d,i) { return ysize + 2.5 * fs + margin.top; } )
		.style("font", 2 * fs + "px sans-serif" )
		.style("text-anchor", "start" ) 
		.text(function(d) { return d;}
		);      
	
	
	if ( dlegend == "yes" )
	{
		var cellheight = ysize / ordlist.length;
		var cellwidth = fs * 1.3;
		var bwidth = diameter * 0.7;
		svg.selectAll(".coco1")
			.data( ordlist )
			.enter().append("text")
			.attr("class", "coco1")
			.attr("x", function(d,i) { return bwidth + cellwidth + fs; } )
			.attr("y", function(d,i) { return margin.top + i * cellheight + 2 * fs; } )
			.style("font", 1.6 * fs + "px sans-serif" )
			.style("text-anchor", "start" ) 
			.text(function(d) { return d;}
			);      

		svg.selectAll(".coco3")
			.data( ordlist )
			.enter().append("rect")
			.attr("class", "coco3")
			.attr("x", function(d,i) { return bwidth; })
			.attr("y", function(d,i) { return margin.top + i * cellheight; } )
			.attr("width", cellwidth )
			.attr("height", cellheight )
			.style("stroke",  "black")
			.style("stroke-width", "0.5px" )
			.style("fill-opacity", 0.5 )
			.style("fill", function(d) { return color( d ); } 
			);
	}

      
"""
elif dprestype == "treemapslice":
	print """

	var svg0 = d3.select( val_id );
	var diameter = +svg0.attr("width") * 0.75;
	var fs = diameter / 100;
	var margin = { top: 8 * fs, right: 2 * fs, bottom: 16 * fs, left: 8 * fs };

	var svg = svg0.append("g")
		.attr("transform", "translate(" + margin.left + "," + margin.top + ")");

	var g = svg.append("g").attr("transform", "translate(2,2)"),
	format = d3.format(",d");

	var hroot = d3.hierarchy(root);
	//hroot.sum(d => d.value);
	
	var treemap = d3.treemap()
		.size([diameter - 4, diameter - 4]);
	treemap.tile( d3[ 'treemapSliceDice' ] );
	var tdata = treemap(root);

	var maxbottom = 0;

	var node = g.selectAll(".node")
		.data(tdata.descendants())
		.enter().append("g")
			.attr("class", function(d) { return d.children ? "node" : "leaf node"; })
			.attr("transform", function(d) { return "translate(" + d.x0 + "," + d.y0 + ")"; });


	var border = fs;
	node.append("rect")
		.attr('width', d => Math.max(0, d.x1 - d.x0 - 2 * border ) + "px" )
		.attr('height', d => d.y1 - d.y0)
		.attr('dx', border )
		.style( "stroke", "rgb(31, 119, 180)" )
		.style( "stroke-width", function(d) 
		{ 
			if ( d.data.name == "root" ) 
			{
				return "0px";
			}
			else if ( d.children )
			{
				return "1px";
			}
			else 
			{
				//maxbottom = Math.max( maxbottom, d.y + d.r ); 				
				return "1px"; 
		} } )
		.style( "fill-opacity", function(d) 
		{ 
			if ( d.data.name == "root" ) 
			{
				return "0.0";
			}
			else if ( d.children )
			{
				return "0.0";
			}
			else return "0.5"; 
		})
		.style("fill", function(d) 
		{ 
			/*
			if ( d.data.name == "root" ) 
			{
				return d3.hsl( 0.3, 0.3, 0.1 );
			}
			else if ( d.children )
			{
				return d3.hsl( 0.3, 0.3, 0.2 );
			}
			else 
			*/
			return color(d.data.order); 
		});
			
	/*
	var smap = d3.scaleOrdinal().domain( [ 1, 2, 3, 4, 5 ] ).range( [ 1., 1.2, 1.25, 1.3, 1.33 ] );
	function ssmap( v ) {	if ( v <= 5 ) return smap( v );	return 1000; }
	var maxfs = fs * 0.2 * Math.sqrt( 8. / fs );
	var minfs = maxfs * 0.5;
  
	function scaledsize( radius, textlength, dimparam )
	{
		var curfs = maxfs * Math.sqrt( radius );
		var lenfs = 5 * dimparam * radius / textlength;
		return Math.max( dimparam * curfs, Math.min( curfs, lenfs ) ); 
	}
	*/

	node.filter(function(d) { return !d.children; }).append("text")
		//.attr( "dy", function(d) { return -0.2 * scaledsize( d.r, 1, 0.7 ) + "px"; } ) 
		.style( "font", function(d) { return 0.2 * fs + "px sans-serif"; } )
		.style( "text-anchor", "right" )
		.text( function(d) { return d.data.sname; } );
		/*
		.each(function(d,i) 
		{
			d.drawtext = true;
			var thisWidth = this.getComputedTextLength();
			if ( thisWidth > 1.95 * d.r )
			{
				this.remove();
				d.drawtext = false;
			}
		}
		);
		*/

  
	if ( dglabels == "yes" )
	{
		node.filter(function(d) { return d.children && d.data.name != "root"; }).append("text")
			.attr("dy", function( d ) { return ( d.y1 + fs * 2 ); } )
			.attr("dx", fs )
			.style("font", fs * 2. + "px sans-serif" )
			.style("font-weigth", "bolder" )
			.style( "text-anchor", "start" )
			.text(function(d) { return d.data.sname; } 
			);
	}
	
	if ( dtlabels == "yes" )
	{
		node.filter(function(d) { return !d.children; }).append("text")
			.attr( "dx", function(d) { return 1.2 * border + "px"; } ) 
			.attr( "dy", function(d) { return 1.5 * fs + "px"; } ) 
			.style( "font", function(d) { fs + "px sans-serif"; } )
			.style( "font-weigth", "bolder" )
			.style( "text-anchor", "start" )
			.text(function(d) { return d.data.percent; } 
			);
	}
	
	
	if ( dlegend == "yes" )
	{
		var cellwidth = fs * 11;
		var cellheight = fs * 1.3;
		var height = diameter + cellheight * 4;
		svg.selectAll(".coco1")
			.data( ordlist )
			.enter().append("text")
			.attr("class", "coco1")
			.attr("x", function(d,i) { return i * cellwidth; } )
			.attr("y", function(d,i) { return height + 5 * cellheight + ( i % 3 ) * fs * 1.8; } )
			.style("font", 1.6 * fs + "px sans-serif" )
			.style("text-anchor", "start" ) 
			.text(function(d) { return d;}
			);      

		svg.selectAll(".coco3")
			.data( ordlist )
			.enter().append("rect")
			.attr("class", "coco3")
			.attr("x", function(d,i) { return i * cellwidth; })
			.attr("y", height + 2 * cellheight )
			.attr("width", cellwidth )
			.attr("height", cellheight )
			.style("stroke",  "black")
			.style("stroke-width", "0.5px" )
			.style("fill-opacity", 0.5 )
			.style("fill", function(d) { return color( d ); } 
			);
	}

      
"""
else:
	print """
       
  var div = d3.select( val_id )
    //.attr('id', "mainsvg" )
    .style("position", "relative")
    //.style("width", (width + margin.left + margin.right) + "px")
    //.style("height", (height + margin.top + margin.bottom) + "px")
    //.style("left", margin.left + "px")
    //.style("top", margin.top + "px");
    
  width = +div.attr( "width" );
  height = +div.attr( "height" );

  var treemap = d3.treemap().size([width, height]).padding( 2 );
  //treemap.tile( d3[ 'treemapSliceDice' ] )
  var tree = treemap(root);
  
    
  var node = div.datum(root).selectAll(".node")
      .data(tree.leaves())
    .enter().append("div")
      .attr("class", "node")
      .style("left", (d) => d.x0 + "px")
      .style("top", (d) => d.y0 + "px")
      .style("width", (d) => Math.max(0, d.x1 - d.x0 - 1) + "px")
      .style("height", (d) => Math.max(0, d.y1 - d.y0  - 1) + "px")
      .style("background", (d) => color(d.data.order))
      .style("border", "solid 1px white")
      .text((d) => d.data.name);
      
	d3.selectAll("input").on("change", function change() {
    const value = ( 1 == 0 )
        ? (d) => { return d.size ? 1 : 0;}
        : (d) => { return d.size; };

    const newRoot = d3.hierarchy(data, (d) => d.children)
      .sum(value);

    node.data(treemap(newRoot).leaves())
      .transition()
        .duration(1500)
        .style("left", (d) => d.x0 + "px")
        .style("top", (d) => d.y0 + "px")
        .style("width", (d) => Math.max(0, d.x1 - d.x0 - 1) + "px")
        .style("height", (d) => Math.max(0, d.y1 - d.y0  - 1) + "px")
  });





"""

    
print """


}


draw( "#normal" );

function click(d) {
  console.log( "bb " );
  if (d.title) {
	console.log( " mm " + d.title );
	d._title = d.title;
    d.title = null;
  } else {
	console.log( " nn " + d._title );
    d.title = d._title;
    d._title = null;
  }
  update( d );
}


	
 </script>
"""