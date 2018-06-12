
import json
from django import forms
from . import controls

class UploadFile( forms.Form ):
	name = forms.CharField(max_length=50, label="Dataset name")
	file = forms.FileField( label="File (tab-delimited / biom format)" )
	
	
class VChoiceField(forms.ChoiceField):
	def valid_value(self, value):
		return True
	def validate(self, value):
		return super( VChoiceField, self ).validate( value )
		
		
class GenericForm( forms.Form ):
	requested_asset = None
	
	def __init__(self,*args,**kwargs):
		job_id = kwargs.pop( 'job_id')
		super(GenericForm,self).__init__(*args,**kwargs)
		summary = json.loads( controls.run_script( "", job_id, "summary" ) )
		tags = json.loads( controls.run_script( "", job_id, "tags" ) )
		taxfilters = json.loads( controls.run_script( "", job_id, "taxonomy" ) )
		cvolumes = [ ( v, v ) for v in summary[ 'volumes' ] ]
		levelnames = { "qiime" : [ "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU" ] }
		if summary[ "taxtype" ] != "none": 
			clevels = [ ( str(k), levelnames[ summary[ "taxtype" ] ][ k - 1 ] + " (" + str( k ) + ")" ) for k in range( 1, summary[ 'maxlevel' ] + 2 ) ]
		else:
			clevels = [ ( str(k), str(k) ) for k in range( 1, summary[ 'maxlevel' ] + 1 ) ]
		ctags = [ ( v, v ) for v in tags.keys() ]
		ctaxfilters = [ ( v, v ) for v in taxfilters.keys() ]
		if 'level' in self.fields:
			self.fields[ 'level' ].widget.choices += clevels 
		
		for fieldname in [ 'dorder', 'volume' ]:
			if fieldname in self.fields:
				self.fields[ fieldname ].widget.choices += cvolumes
		for fieldname in [ 'dfilter', 'color', 'shape', 'dgroup', 'order1', 'order2', 'labels' ]:
			if fieldname in self.fields:
				self.fields[ fieldname ].widget.choices += ctags
		for fieldname in [ 'spfilter', 'spshow' ]:
			if fieldname in self.fields:
				self.fields[ fieldname ].widget.choices += ctaxfilters
				

class Table( GenericForm ):
	level = VChoiceField()
	dtype = forms.ChoiceField( label = "Units", choices = [ ( "count", "Counts" ), ( "percent", "Percents" ) ] )
	dorder = VChoiceField( label="Sort order", choices = [ ( "none", "Unchanged" ), ( "taxonomy", "Taxonomy" ) ] )
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dgroup = VChoiceField( label="Samples grouping" )
	numbest = forms.ChoiceField( label = "Number of best species / units", choices = [ ( "all", "all" ), ( "10", "10" ), ( "25", "25"), ( "50", "50" ), ( "100", "100" ) ] )

class Indices( GenericForm ):
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dgroup = VChoiceField( label="Samples grouping" )
	
class Anova( GenericForm ):
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dgroup = VChoiceField( label="Samples grouping" )
	cmethod = forms.ChoiceField( label="Method", choices = [ ( "anova", "Anova" ), ( "best-ttest", "T-test (best pair)" ), ( "best-wilcoxon", "Wilcoxon (best pair)" ) ] )
	cunits = forms.ChoiceField( label = "Units", choices = [ ( "probability", "p-value" ), ( "log-probability", "-log( p-value )" ) ] )

class Permanova( GenericForm ):
	level = VChoiceField()
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dgroup = VChoiceField( label="Samples grouping" )
	pmethod = forms.ChoiceField( label="Method", choices = [ ( "permanova", "Permanova" ), ( "anosim", "Anosim" ) ] )
	cunits = forms.ChoiceField( label = "Units", choices = [ ( "probability", "p-value" ), ( "log-probability", "-log( p-value )" ) ] )

class Heatmap( GenericForm ):
	level = VChoiceField()
	dtype = forms.ChoiceField( label = "Units", choices = [ ( "count", "Counts" ), ( "percent", "Percents" ), ( "z-score", "Z-Score" ) ] )
	spfilter = VChoiceField( label="Taxonomy filter" )
	order1 = VChoiceField( label="Sort order (primary)" )
	order2 = VChoiceField( label="Sort order (secondary)" )
	labels = VChoiceField( label="Labels for samples" )
	numbest = forms.ChoiceField( label = "Number of best species / units", choices = [ ( "all", "all" ), ( "10", "10" ), ( "25", "25"), ( "50", "50" ), ( "100", "100" ) ] )
	numhighligtht = forms.ChoiceField( label = "Highlight best species/units", choices = [ ( v, v ) for v in [ "all", "20%", "10%", "7.5%", "5%", "2%" ] ] )
	spshow = VChoiceField( label="Show taxonomy", choices = [ ( "custom", "Custom" ) ] )
	taxlabels = forms.ChoiceField( label = "Labels for taxonomy", choices = [ ( "no", "No" ), ( "yes", "Yes" ) ] )
																		  
class BubbleChart( GenericForm ):
	level = VChoiceField()
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dgroup = VChoiceField( label="Samples grouping" )
	dnorm = forms.ChoiceField( label = "Units", choices = [ ( "percent", "Percents" ), ( "count", "Counts" ) ] )
	dglabels = forms.ChoiceField( label = "Group labels", choices = [ ( "no", "No" ), ( "yes", "Yes" ) ] )
	dplabels = forms.ChoiceField( label = "Percent labels", choices = [ ( "no", "No" ), ( "yes", "Yes" ) ] )
	dlegend = forms.ChoiceField( label = "Show legend", choices = [ ( "no", "No" ), ( "yes", "Yes" ) ] )
	dprestype = forms.ChoiceField( label = "Presentation type", choices = [ ( "bubble", "Bubble chart" ), ( "treemap", "Treemap" ) ] )