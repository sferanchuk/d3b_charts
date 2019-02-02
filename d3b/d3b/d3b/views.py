# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.template import loader
from django import forms
from . import controls
from . import forms as dforms
import json
import urllib
import re

# Create your views here.
jvenn_site = "http://jvenn.toulouse.inra.fr/app/"
jvenn_css = ""
#"\n".join( [ "<link href=\"" + jvenn_site + "css/" + css + ".css\" rel=\"stylesheet\">" for css in [ "bootstrap",  "bootstrap-responsive", "bootstrap-colorpicker.min" ] ] )
#jvenn_css = "\n".join( [ "<link href=\"" + jvenn_site + "css/" + css + ".css\" rel=\"stylesheet\">" for css in [ "bootstrap", "prettify", "bootstrap-responsive", "bootstrap-colorpicker.min" ] ] )
jvenn_js = "\n".join( [	"<script type=\"text/javascript\" src=\"" + jvenn_site + "js/" + js + ".js\"></script>" for js in [ "jquery.min", "bootstrap.min", "bootstrap-colorpicker.min", "canvas2svg", "jvenn.min" ] ] )

def handle_uploaded_file(f):
	with open( 'pool/d3btestdata.txt', 'wb+') as destination:
		for chunk in f.chunks():
			destination.write(chunk)

def new(request):
	if request.method == 'POST':
		form = dforms.UploadFile(request.POST, request.FILES)
		if form.is_valid():
			job = controls.submit_job( request.FILES, form.data[ 'name' ] )
			return HttpResponseRedirect( '/summary/' + job )
	else:
		form = dforms.UploadFile()
	template = loader.get_template('new.html')
	context = { 'header': "Submit new matrix", 'form' : form, 'title' : 'submit'  }
	return HttpResponse( template.render( context, request )  )

def format_params( form ):
	params = ""
	print form.data
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
	jobtitle = summary[ u'name' ]
	template = loader.get_template('summary.html')
	vollist = ", ".join( summary[ 'volumes' ] )
	result = '<br>'.join( [ "maxlevel: " + str( summary[ 'maxlevel' ] ), "taxonomy type: " + summary[ 'taxtype' ], "num_volumes: " + str( summary[ 'numvolumes' ] ), "volumes: " + vollist ] )
	context = { 'job': job, 'title' : jobtitle, 'service' :  'summary', 'servicename' : 'Summary of the dataset', 'result' : result  }
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
			print form.data[ "command" ]
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
		print script
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
	servicename = 'Define tags'
	tags = controls.run_script( "", job, "tags" )
	if request.method == 'POST':
		re_pattern = re.compile(u'[^\u0000-\uD7FF\uE000-\uFFFF]', re.UNICODE)
		ntstr = re_pattern.sub( '_', u''.join( [ urllib.unquote_plus( request.POST[ 'jstags' ] ) ] ) )
		ntstr = "".join( [ uc if uc < u"\u00ff" else "_" for uc in urllib.unquote_plus( request.POST[ 'jstags' ] ) ] )
		print ntstr
		newtags = controls.run_script( "newtags=" + ntstr, job, "tags" )
		tags = newtags.replace( "u'", "'" )
	ptags = json.loads( tags )
	ptkeys = ptags.keys()
	ptkeys.remove( "name" )
	ptkeys.remove( "none" )
	volumes = ptags[ 'name' ]
	template = loader.get_template( 'tags.html' )
	context = { 'job': job, 'title' : jobtitle, 'service' : service, 'servicename' : servicename, 'tags' : tags, 'ptkeys' : ptkeys, 'volumes' : json.dumps( volumes ), 'pvolumes' : volumes }
	return HttpResponse( template.render( context, request ) )

def taxonomy(request,job):
	jobtitle = job_name( job )
	service = 'taxonomy'
	servicename = 'Define taxonomy filters'
	taxfilters = controls.run_script( "", job, "taxonomy" )
	treedata = controls.run_script( "", job, "get_taxonomy" )
	if request.method == 'POST':
		print request.POST[ 'jsfilters' ]
		newfilters = controls.run_script( "newfilters=" + urllib.unquote_plus( request.POST[ 'jsfilters' ] ), job, "taxonomy" )
		taxfilters = newfilters.replace( "u'", "'" )
	pfilters = json.loads( taxfilters )
	pfkeys = pfilters.keys()
	pfkeys.remove( "none" )
	template = loader.get_template( 'taxonomy.html' )
	context = { 'job': job, 'title' : jobtitle, 'service' : service, 'servicename' : servicename, 'taxfilters' : taxfilters, 'pfkeys' : pfkeys, 'treedata' : treedata, 'jscripts' : [ "d3.v3.min.js" ] }
	return HttpResponse( template.render( context, request ) )

def table(request,job):
	return generic_view( request, job, "table", dforms.Table, template = "table" )

def indices(request,job):
	return generic_view( request, job, "indices", dforms.Indices )

def anova(request,job):
	return generic_view( request, job, "anova", dforms.Anova, servicename = "Variance of alpha-diversity indices" )

def permanova(request,job):
	return generic_view( request, job, "permanova", dforms.Permanova, servicename = "Variance of distances" )

def pca(request,job):
	return generic_view( request, job, "pca", dforms.PCA, jscripts = [ "d3.v3.min.js" ]  )

def tree(request,job):
	return generic_view( request, job, "tree", dforms.Tree, jscripts = [ "d3.v3.min.js" ] )

def heatmap(request,job):
	return generic_view( request, job, "heatmap", dforms.Heatmap, template="heatmap", jscripts = [ "d3.v3.min.js" ] )

def venn(request,job):
	return generic_view( request, job, "venn", dforms.Venn, template = "mchoice_chart", jscripts = [ "d3.v3.min.js", "venn.js" ] )

def ternary(request,job):
	return generic_view( request, job, "ternary", dforms.Ternary, template = "mchoice_chart", jscripts = [ "d3.v3.min.js" ] )

def bubbles(request,job):
	return generic_view( request, job, "bubbles", dforms.BubbleChart, script="bubble", jscripts = [ "d3.v4.min.js" ] )

def whittaker(request,job):
	return generic_view( request, job, "whittaker", dforms.Whittaker, jscripts = [ "d3.v4.min.js" ] )

def rarefunction(request,job):
	return generic_view( request, job, "whittaker", dforms.Rarefunction, jscripts = [ "d3.v4.min.js" ] )

def pca_sp(request,job):
	return generic_view( request, job, "pca_sp", dforms.PCA2P, servicename = "Two-panel PCA", jscripts = [ "d3.v3.min.js" ] )

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
