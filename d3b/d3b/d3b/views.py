# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.template import loader
from django import forms
from . import controls
from . import forms as dforms
import json

# Create your views here.

def handle_uploaded_file(f):
	with open( 'pool/d3btestdata.txt', 'wb+') as destination:
		for chunk in f.chunks():
			destination.write(chunk)

def new(request):
	if request.method == 'POST':
		form = dforms.UploadFileForm(request.POST, request.FILES)
		if form.is_valid():
			job = controls.submit_job( request.FILES['file'], form.data[ 'name' ] )
			return HttpResponseRedirect( '/summary/' + job )
	else:
		form = dforms.UploadFile()
	template = loader.get_template('new.html')
	context = { 'header': "Submit new matrix", 'form' : form, 'title' : 'submit'  }
	return HttpResponse( template.render( context, request )  )


def run_script( form, job, script ):
	params = ""
	for field in list( form.declared_fields ):
		value = form.data[ field ]
		params += field + "=" + value + "&"
	return controls.run_script( params, job, script )

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
	template = kwargs.get( 'template', "default_service" )
	outname = jobtitle + "_" + service
	result = "select parameters and press <i>Submit</i><br>"
	if request.method == 'POST':
		form = formclass( request.POST, request.FILES, job_id=job )
		result = run_script( form, job, script  )
		if not service in [ "indices", "summary", "anova" ]:
			outname += "_" + form.data[ 'level' ]
	else:
		form = formclass( job_id=job )
	template = loader.get_template( template + '.html' )
	context = { 'job': job, 'title' : jobtitle, 'service' : service, 'servicename' : servicename, 'result' : result, 'form' : form , 'outname' : outname, 'jscripts' : jscripts }
	return HttpResponse( template.render( context, request ) )

def tags(request,job):
	return generic_view( request, job, "tags", dforms.GenericForm )

def taxonomy(request,job):
	return generic_view( request, job, "taxonomy", dforms.GenericForm )

def table(request,job):
	return generic_view( request, job, "table", dforms.Table, template = "table" )

def indices(request,job):
	return generic_view( request, job, "indices", dforms.Indices )

def anova(request,job):
	return generic_view( request, job, "anova", dforms.Anova, servicename = "Variance of alpha-diversity indices" )

def permanova(request,job):
	return generic_view( request, job, "permanova", dforms.Permanova, servicename = "Variance of distances" )

def pca(request,job):
	return generic_view( request, job, "pca", dforms.GenericForm )

def tree(request,job):
	return generic_view( request, job, "tree", dforms.GenericForm )

def heatmap(request,job):
	return generic_view( request, job, "heatmap", dforms.Heatmap, jscripts = [ "d3.v3.min.js" ], template = 'heatmap' )

def venn(request,job):
	return generic_view( request, job, "venn", dforms.GenericForm )

def ternary(request,job):
	return generic_view( request, job, "ternary", dforms.GenericForm )

def bubbles(request,job):
	return generic_view( request, job, "bubbles", dforms.BubbleChart, script="bubble" )

def whittaker(request,job):
	return generic_view( request, job, "pca", dforms.GenericForm )


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
