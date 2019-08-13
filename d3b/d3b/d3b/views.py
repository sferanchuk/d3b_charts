
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.template import loader
from django import forms
import controls
import dforms
import json
import urllib.request, urllib.parse, urllib.error
import re

jvenn_site = "http://jvenn.toulouse.inra.fr/app/"
jvenn_css = ""
jvenn_js = "\n".join( [	"<script type=\"text/javascript\" src=\"" + jvenn_site + "js/" + js + ".js\"></script>" for js in [ "jquery.min", "bootstrap.min", "bootstrap-colorpicker.min", "canvas2svg", "jvenn.min" ] ] )

def handle_uploaded_file(f):
	with open( 'pool/d3btestdata.txt', 'wb+') as destination:
		for chunk in f.chunks():
			destination.write(chunk)

def new(request):
	if request.method == 'POST':
		form = dforms.UploadFile(request.POST, request.FILES)
		if form.is_valid():
			job = controls.submit_job( request.FILES, form.data[ 'name' ], form.data[ 'transform' ] )
			return HttpResponseRedirect( '/summary/' + job )
	else:
		form = dforms.UploadFile()
	template = loader.get_template('new.html')
	context = { 'header': "Starting page", 'form' : form, 'title' : 'submit'  }
	return HttpResponse( template.render( context, request )  )

def format_params( form ):
	params = ""
	print(form.data)
	for field in list( form.declared_fields ):
		values = form.data.getlist( field )
		for value in values:
			if len( value ) > 0:
				params += field + "=" + value + "&"
	return params[:-1]

def run_script( form, job, script ):
	return controls.run_script( format_params( form ), job, script )

def job_name( job ):
	summary = json.loads( controls.run_script( "", job, "summary" ) )
	return summary[ 'name' ]

def summary(request,job):
	summary = json.loads( controls.run_script( "", job, "summary" ) )
	tags = json.loads( controls.run_script( "", job, "tags" ) )
	jobtitle = summary[ 'name' ]
	template = loader.get_template('summary.html')
	vollist = ", ".join( summary[ 'volumes' ] )
	result = '<br>'.join( [ "maxlevel: " + str( summary[ 'maxlevel' ] ), "taxonomy type: " + summary[ 'taxtype' ], "num_volumes: " + str( summary[ 'numvolumes' ] ), "num_records: " + str( summary[ 'numrecords' ] ), "volumes: " + vollist ] )
	context = { 'job': job, 'title' : jobtitle, 'service' :  'summary', 'servicename' : 'summary of the dataset', 'result' : result  }
	return HttpResponse( template.render( context, request ) )

def generic_view( request, job, service, formclass, **kwargs ):
	jobtitle = job_name( job )
	jscripts = kwargs.get( 'jscripts', [] )
	servicename = kwargs.get( 'servicename', service )
	script = kwargs.get( 'script', service )
	template = kwargs.get( 'template', "default_service" if len( jscripts ) == 0 else "default_chart" )
	outname = jobtitle + "_" + service
	result = "select parameters and press <i>Submit</i><br>"
	tags = controls.run_script( "", job, "tags" )
	if request.method == 'POST':
		form = formclass( request.POST, request.FILES, job_id = job, tags_dict = tags )
		if not service in [ "indices", "summary", "anova" ]:
			outname += "_" + form.data[ 'level' ]
		if "command" in form.data and form.data[ "command" ] != "Submit":
			host = "http://secure.bri-shur.com:%s/static" % request.META[ 'SERVER_PORT' ]
			print(form.data[ "command" ])
			#host = "http://bri-shur.com/javascripts"
			if form.data[ "command" ] == "Download as PNG":
				pngres = controls.render_png( format_params( form ), job, script, host, jscripts )
				if len( pngres ) > 0:
					response = HttpResponse( pngres, content_type="image/png")
					response[ 'Content-Disposition' ] = 'attachment; filename="%s.png"' % outname.replace( "\n", "_" )
					return response
			if form.data[ "command" ] == "Download as SVG":
				svgres = controls.render_svg( format_params( form ), job, script, host, jscripts )
				if len( svgres ) > 0:
					response = HttpResponse( svgres, content_type="image/svg+xml")
					response[ 'Content-Disposition' ] = 'attachment; filename="%s.svg"' % outname.replace( "\n", "_" )
					return response
		print(script)
		result = run_script( form, job, script  )
	else:
		form = formclass( job_id=job, tags_dict = tags )
	template = loader.get_template( template + '.html' )
	context = { 'job': job, 'title' : jobtitle, 'service' : service, 'servicename' : servicename, 'result' : result, 'form' : form , 'outname' : outname, 'jscripts' : jscripts, 'tags' : tags }
	if service == "venn" and request.method == 'POST' and form.data[ "ptype" ] == "interactive-presence":
		context[ 'extracss' ] = jvenn_css
		context[ 'extrascripts' ] = jvenn_js
	return HttpResponse( template.render( context, request ) )

def tags(request,job):
	jobtitle = job_name( job )
	service = 'tags'
	servicename = 'subsets of samples'
	tags = controls.run_script( "", job, "tags" )
	if request.method == 'POST':
		re_pattern = re.compile('[^\u0000-\uD7FF\uE000-\uFFFF]', re.UNICODE)
		ntstr = re_pattern.sub( '_', ''.join( [ urllib.parse.unquote_plus( request.POST[ 'jstags' ] ) ] ) )
		ntstr = "".join( [ uc if uc < "\u00ff" else "_" for uc in urllib.parse.unquote_plus( request.POST[ 'jstags' ] ) ] )
		print(ntstr)
		newtags = controls.run_script( "newtags=" + ntstr, job, "tags" )
		tags = newtags.replace( "u'", "'" )
	ptags = json.loads( tags )
	ptkeys = list(ptags.keys())
	ptkeys.remove( "name" )
	ptkeys.remove( "none" )
	volumes = ptags[ 'name' ]
	template = loader.get_template( 'tags.html' )
	context = { 'job': job, 'title' : jobtitle, 'service' : service, 'servicename' : servicename, 'tags' : tags, 'ptkeys' : ptkeys, 'volumes' : json.dumps( volumes ), 'pvolumes' : volumes }
	return HttpResponse( template.render( context, request ) )

def taxonomy(request,job):
	jobtitle = job_name( job )
	service = 'taxonomy'
	servicename = 'taxonomic groups'
	taxfilters = controls.run_script( "", job, "taxonomy" )
	treedata = controls.run_script( "", job, "get_taxonomy" )
	if request.method == 'POST':
		print(request.POST[ 'jsfilters' ])
		newfilters = controls.run_script( "newfilters=" + urllib.parse.unquote_plus( request.POST[ 'jsfilters' ] ), job, "taxonomy" )
		taxfilters = newfilters.replace( "u'", "'" )
	pfilters = json.loads( taxfilters )
	pfkeys = list(pfilters.keys())
	pfkeys.remove( "none" )
	template = loader.get_template( 'taxonomy.html' )
	context = { 'job': job, 'title' : jobtitle, 'service' : service, 'servicename' : servicename, 'taxfilters' : taxfilters, 'pfkeys' : pfkeys, 'treedata' : treedata, 'jscripts' : [ "d3.v3.min.js" ] }
	return HttpResponse( template.render( context, request ) )

def table(request,job):
	return generic_view( request, job, "table", dforms.Table, template = "table", servicename = "matrix of abundance values" )

def indices(request,job):
	return generic_view( request, job, "indices", dforms.Indices, servicename = "alpha-diversity indices" )

def anova(request,job):
	return generic_view( request, job, "anova", dforms.Anova, servicename = "variance of alpha-diversity indices" )

def permanova(request,job):
	return generic_view( request, job, "permanova", dforms.Permanova, servicename = "variance of distances" )

def pca(request,job):
	return generic_view( request, job, "pca", dforms.PCA, jscripts = [ "d3.v3.min.js" ], servicename = "PCA / constrained ordination"  )

def tree(request,job):
	return generic_view( request, job, "tree", dforms.Tree, jscripts = [ "d3.v3.min.js" ], template = "tree", servicename = "dendrogramm of distances"  )

def heatmap(request,job):
	return generic_view( request, job, "heatmap", dforms.Heatmap, template="heatmap", jscripts = [ "d3.v3.min.js" ], servicename = "heatmap presentation"  )

def venn(request,job):
	return generic_view( request, job, "venn", dforms.Venn, template = "mchoice_chart", jscripts = [ "d3.v3.min.js", "venn.js" ], servicename = "venn diagramm"  )

def ternary(request,job):
	return generic_view( request, job, "ternary", dforms.Ternary, template = "mchoice_chart", jscripts = [ "d3.v3.min.js" ], servicename = "ternary chart"  )

def bubbles(request,job):
	return generic_view( request, job, "bubbles", dforms.BubbleChart, script="bubble", jscripts = [ "d3.v4.min.js" ], servicename = "bubble chart"  )

def whittaker(request,job):
	return generic_view( request, job, "whittaker", dforms.Whittaker, jscripts = [ "d3.v4.min.js" ], servicename = "whittaker chart / rarefaction curve"  )

def volcano(request,job):
	return generic_view( request, job, "volcano", dforms.Volcano, template = "mchoice_chart", jscripts = [ "d3.v4.min.js" ], servicename = "volcano chart / MD-plot"  )

def pca_sp(request,job):
	return generic_view( request, job, "pca_sp", dforms.PCA2P, jscripts = [ "d3.v3.min.js" ], servicename = "two-panel PCA"  )

def network(request,job):
	return generic_view( request, job, "network", dforms.Network, jscripts = [ "d3.v3.min.js" ], servicename = "network of species"  )
	
if False:
	jobtitle = job_name( job )
	result = ""
	if request.method == 'POST':
		form = dforms.BubbleChart(request.POST, request.FILES)
		result = run_script( form, job, "bubble"  )
	else:
		form = dforms.BubbleChart()
	template = loader.get_template('bubbles.html')
	context = { 'job': job, 'title' : jobtitle, 'service' :  'bubbles', 'servicename' : 'Bubble chart', 'result' : result, 'form' : form  }
#	return HttpResponse( template.render( context, request ) )
