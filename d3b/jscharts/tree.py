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
from sklearn import manifold
from sklearn import decomposition
import scipy.stats as stats
from scipy.spatial import distance
import json
import os.path
import sqlite3
import time
import math

import d3bf

form = cgi.FieldStorage()
id = "emap"
labels = form.getvalue( "labels", "name" )
dfilter = form.getvalue( "dfilter", "none" )
spfilter = d3bf.loadfilters( "emap_filters.txt", form.getvalue( "spfilter", "none" ) )
dmethod = form.getvalue( "dmethod", "Pearson" )
tmethod = form.getvalue( "tmethod", "upgma" )
lmethod = form.getvalue( "lmethod", "lines" )
level = form.getvalue( "level" )
resolution = form.getvalue( "resolution", "low" )
jobname = form.getvalue( "jobname", "noname" )

ilevel = int( level ) 
( data, volumes, mn, ml ) = d3bf.loaddata( "emap.txt" )
( tags, tkeys ) = d3bf.loadtags( "emap_tags.txt", volumes )
( kdict, kdnames, kgnames, knorder, kdata ) = d3bf.loadtaxonomy( data, ml, spfilter, ilevel )
( findex, mtags ) = d3bf.processtags_m( volumes, tags, dfilter )
( edata, site_ids, species_ids ) = d3bf.load_edata_m( data, ilevel, mn, ml, kdict, volumes, findex, kdnames )

aedata = np.array( edata, dtype=float )
aenorm = np.sum( aedata, axis=1 )
aedata /= aenorm.reshape( len(edata), 1 )

rev = 0
cdata = d3bf.calc_distances( edata, aedata, dmethod, kdata, knorder, rev )

if resolution == "high":
	print "<svg width=\"2400\" height=\"%d\" id=\"normal\"></svg>" % ( len( kdict ) * 40 )
else:
	print "<svg width=\"1200\" height=\"%d\" id=\"normal\"></svg>" % ( len( kdict ) * 20 )
print "<script>"
print "var tmethod = \"" + tmethod + "\";"
print "var lmethod = \"" + lmethod + "\";"
print "var jobname = \"" + jobname + "\";"
print "var dist = ["
for i in range( len( cdata ) ):
	ch = "," if i + 1 < len( cdata ) else "];"
	print json.dumps( cdata[i] ) + ch
print "var datal =  %s;" % json.dumps( mtags[ labels ] )

print """

console.log( datal.length )

  function clustering( pairwiseDistances, linkage ) 
   {
        var clusters, minDistance, clusterA, clusterB, distance, distanceA,
            distanceB, candidates, mergedCluster, i, j, k;
         clusters = [];
        for (i = 0; i < pairwiseDistances.length; i++) {
            clusters.push({
                name: i,
                key: i,
                index: i,
                size: 1,
                label: datal[i],
                //parent: "root",
                distance: 0.
            });
        }
        k = i;
        while (true) {
         // Stop if all clusters have been merged into a single cluster.
            if (clusters.length === 1) {
                delete clusters[0].index;
                delete clusters[0].key;
                break;
            }
         // Find closest clusters.
            minDistance = Infinity;
            for (i = 0; i < clusters.length; i++) {
                clusterA = clusters[i];
                for (j = 0; j < clusters.length; j++) {
                    if (i !== j) {
                        clusterB = clusters[j];
                        distance = pairwiseDistances[clusterA.key][clusterB.key];
                        if (distance < minDistance) {
                            minDistance = distance;
                            candidates = [clusterA, clusterB];
                        }
                    }
                }
            }
         // Merge clusters.
            mergedCluster = {
                children: candidates,
                key: candidates[0].key,
                distance: minDistance,
                branch_length: 0,
                size: candidates[0].size + candidates[1].size,
                label: "",
                name: k
                //,parent: "root"
            };
            candidates[0].branch_length = minDistance - candidates[0].distance;
            candidates[1].branch_length = minDistance - candidates[1].distance;
            //candidates[0].parent = candidates[1].parent = k;
            k = k + 1;
         // Replace first cluster with merged cluster in list of clusters.
            clusters[candidates[0].index] = mergedCluster;
         // Remove second cluster from list of clusters.
            clusters.splice(candidates[1].index, 1);
         // Recompute distances from merged cluster to all other clusters.
            for (i = 0; i < clusters.length; i++) {
                if (clusters[i].key === candidates[0].key) {
                    distance = 0;
                } else {
                    distanceA = pairwiseDistances[candidates[0].key][clusters[i].key];
                    distanceB = pairwiseDistances[candidates[1].key][clusters[i].key];
                    switch (linkage) {
                        case "single":
                            if (distanceA < distanceB) {
                                distance = distanceA;
                            } else {
                                distance = distanceB;
                            }
                            break;
                        case "complete":
                            if (distanceA > distanceB) {
                                distance = distanceA;
                            } else {
                                distance = distanceB;
                            }
                            break;
                        case "average":
                            distance = ((distanceA * candidates[0].size) + (distanceB * candidates[1].size)) / (candidates[0].size + candidates[1].size);
                            break;
                    }
                }
                pairwiseDistances[candidates[0].key][clusters[i].key] = distance;
                pairwiseDistances[clusters[i].key][candidates[0].key] = distance;
            }
         // Remove column of second cluster from pairwise distance matrix.
            for (i = 0; i < pairwiseDistances.length; i++) {
                pairwiseDistances[i].splice(candidates[1].key, 1);
            }
         // Remove row of second cluster from pairwise distance matrix.
            pairwiseDistances.splice(candidates[1].key, 1);
         // Update keys of clusters to reflect removal of the column.
            for (i = candidates[1].key; i < clusters.length; i++) {
                clusters[i].key--;
            }
         // Remove obsolete key and index of merged clusters.
            delete candidates[0].key;
            delete candidates[0].index;
            delete candidates[1].key;
            delete candidates[1].index;
         // Reindex clusters.
            for (i = 0; i < clusters.length; i++) {
                clusters[i].index = i;
            }
        }
        return clusters;
    };


   var margin = {top: 20, right: 15, bottom: 100, left: 100}
      , width = 1200 - margin.left - margin.right
      , height = 800 - margin.top - margin.bottom;

var svg0 = d3.select( "#normal" );
   width = +svg0.attr( "width" ) * 0.9;
   var fs = width / 50;
   height = +svg0.attr( "height" ) * 0.9;
	var tsize = datal.length;
	if ( fs * tsize < height ) height = fs * tsize;
	if ( fs * tsize * 1.2 + 200 < width ) width = fs * tsize * 1.2 + 200;
   margin.top = margin.right = fs;
   margin.left = margin.bottom = 5 * fs;
 

var svg = svg0.append("g")
    .attr("transform",
          "translate(" + margin.left + "," + margin.top + ")");    


var tpos = width - 200;

var cluster = d3.layout.cluster()
    .size([height, tpos - 10 ])
    .separation(function(a, b) { return 1; });
    
var treeData = clustering( dist, tmethod );
	console.log( JSON.stringify( treeData ) );

var diagonal = d3.svg.diagonal()
	.projection(function(d) { return [d.y, d.x]; });

//var line = d3.svg.line()
//    .x(function(d) { return x(d.x); })
//    .y(function(d) { return y(d.y); });

root = treeData[0];

var nodes = cluster.nodes( root );
var links = cluster.links( nodes );
console.log( nodes.length )
console.log( links.length )

//console.log( JSON.stringify( nodes ) );

var link = svg.selectAll(".link")
      .data( links )
    .enter().append("path")
      .attr("class", "link")
      .style( "fill", "none" )
	  .style( "stroke", "#555" )
	  .style( "stroke-width", "4px" )
      .attr("d", function( d ) { return ( lmethod == "bezier" ) ? diagonal( d ) : "M" + d.source.y + "," + d.source.x + "L" + d.source.y + "," + d.target.x + "L" + d.target.y + "," + d.target.x; } );
      //.attr("d", diagonal);

  var node = svg.selectAll(".node")
      .data(nodes)
    .enter().append("g")
      .attr("class", "node")
      .attr("transform", function(d) { return "translate(" + d.y + "," + d.x + ")"; })

  node.append("circle")
      .attr("r", 4.);

  node.append("text")
      .attr("dx", function(d) { return d.children ? -8 : 8; })
      .attr("dy", 3)
      .attr("text-anchor", function(d) { return d.children ? "end" : "start"; })
      .style("font", fs + "px sans-serif" )
      .text(function(d) { return d.label; });

d3.select("#gpng")
		.on("click", writeDownloadPng);
d3.select("#gnewick")
		.on("click", exportNewick);
		
		
	function nested(nest){
		var subtree = "";

		if(nest.hasOwnProperty('children')){
			var children = [];
			nest.children.forEach(function(child){
				var subsubtree = nested(child);
				children.push(subsubtree);
			});
      var substring = children.join();
      if(nest.hasOwnProperty('label')){
        subtree = "("+substring+")" + nest.label;
      }
      if(nest.hasOwnProperty('branch_length')){
        subtree = subtree + ":"+nest.branch_length;
      }
		}
		else{
      var leaf = "";
      if(nest.hasOwnProperty('label')){
        leaf = nest.label;
      }
      if(nest.hasOwnProperty('branch_length')){
        leaf = leaf + ":"+nest.branch_length;
      }
      subtree = subtree + leaf;
		}
		return subtree;
	}

	
	function writeDownloadPng(){
		var element = document.getElementById('svg1');
		element.style.background="white";
		html2canvas(element, {
			onrendered: function(canvas) {
				canvas.toBlob(function(blob) {
					saveAs(blob, "tree.png" );
				}, "image/png");
		} } );
	}
	function exportNewick(){
		var str = nested( treeData[0] );
		var blob = new Blob([ str ], {type: "text/plain;charset=utf-8"});
		saveAs(blob, jobname + "_tree.nwk" );
	}


</script>
"""