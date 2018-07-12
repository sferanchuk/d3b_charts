
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
		

class VMChoiceField(forms.MultipleChoiceField):
	def valid_value(self, value):
		return True
	def validate(self, value):
		return super( VMChoiceField, self ).validate( value )


class GenericForm( forms.Form ):
	requested_asset = None
	
	def __init__(self,*args,**kwargs):
		job_id = kwargs.pop( 'job_id' )
		tags = json.loads( kwargs.pop( 'tags_dict' ) )
		super(GenericForm,self).__init__(*args,**kwargs)
		summary = json.loads( controls.run_script( "", job_id, "summary" ) )
		#tags = json.loads( controls.run_script( "", job_id, "tags" ) )
		taxfilters = json.loads( controls.run_script( "", job_id, "taxonomy" ) )
		cvolumes = [ ( v, v ) for v in summary[ 'volumes' ] ]
		levelnames = { "qiime" : [ "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU" ] }
		if summary[ "taxtype" ] != "none": 
			clevels = [ ( str(k), levelnames[ summary[ "taxtype" ] ][ k - 1 ] + " (" + str( k ) + ")" ) for k in range( 1, summary[ 'maxlevel' ] + 2 ) ]
		else:
			clevels = [ ( str(k), str(k) ) for k in range( 1, summary[ 'maxlevel' ] + 1 ) ]
		ctags = [ ( v, v ) for v in tags.keys() ]
		ctaxfilters = [ ( v, v ) for v in taxfilters.keys() ]
		distmeasures = [ "Pearson", "Kendall", "Spearman", "Euclidean", "Bray-Curtis", "Jaccard", "Morisita-Horn", "Unifrac-unweighted", "Unifrac-weighted" ]
		cdmeasures = [ ( v, v ) for v in distmeasures ]
		
		if 'level' in self.fields:
			self.fields[ 'level' ].widget.choices += clevels 
		for fieldname in [ 'dorder', 'volume' ]:
			if fieldname in self.fields:
				self.fields[ fieldname ].widget.choices += cvolumes
		for fieldname in [ 'dmethod' ]:
			if fieldname in self.fields:
				self.fields[ fieldname ].widget.choices += cdmeasures
		for fieldname in [ 'dfilter', 'color', 'shape', 'dgroup', 'order1', 'order2', 'labels' ]:
			if fieldname in self.fields:
				self.fields[ fieldname ].widget.choices += ctags
		for fieldname in [ 'spfilter', 'spshow' ]:
			if fieldname in self.fields:
				self.fields[ fieldname ].widget.choices += ctaxfilters
		if 'samples' in self.fields:
			ctag = 'name'
			if len( args ) > 0:
				post = args[ 0 ]
				if 'dgroup' in post:
					ctag = post[ 'dgroup' ]
					if ctag == 'none':
						ctag = 'name'
			snames = set( tags[ ctag ] )
			print snames
			self.fields[ 'samples' ].widget.choices = [ ( v, v ) for v in snames ]
			
					
			
				

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


class PCA( GenericForm ):
	level = VChoiceField()
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dmethod = VChoiceField( label="Distances" )
	fmethod = forms.ChoiceField( label = "Method", choices = [ ( v, v ) for v in [ "CA (skbio)", "PCoA (skbio)", "PCA (sklearn)", "MDS (sklearn)" ] ] )
	color = VChoiceField( label="Color" )
	shape = VChoiceField( label="Shape" )
	labels = VChoiceField( label="Labels" )

class Tree( GenericForm ):
	level = VChoiceField()
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dmethod = VChoiceField( label="Distances" )
	tmethod = forms.ChoiceField( label = "Clustering method", choices = [ ( v, v ) for v in [ "average", "single", "complete" ] ] )
	lmethod = forms.ChoiceField( label = "Drawing links", choices = [ ( v, v ) for v in [ "lines", "bezier" ] ] )
	labels = VChoiceField( label="Labels" )

class Heatmap( GenericForm ):
	level = VChoiceField()
	dtype = forms.ChoiceField( label = "Units", choices = [ ( "count", "Counts" ), ( "percent", "Percents" ), ( "z-score", "Z-Score" ) ] )
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	order1 = VChoiceField( label="Sort order (primary)" )
	order2 = VChoiceField( label="Sort order (secondary)" )
	labels = VChoiceField( label="Labels for samples" )
	numbest = forms.ChoiceField( label = "Number of best species / units", choices = [ ( "all", "all" ), ( "10", "10" ), ( "25", "25"), ( "50", "50" ), ( "100", "100" ) ] )
	numhighligtht = forms.ChoiceField( label = "Highlight best species/units", choices = [ ( v, v ) for v in [ "all", "20%", "10%", "7.5%", "5%", "2%" ] ] )
	spshow = VChoiceField( label="Show taxonomy", choices = [ ( "custom", "Custom" ) ] )
	taxlabels = forms.ChoiceField( label = "Labels for taxonomy", choices = [ ( "no", "No" ), ( "yes", "Yes" ) ] )
	spcustom = forms.CharField( widget=forms.HiddenInput() )
	
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

class Ternary( GenericForm ):
	level = VChoiceField()
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dgroup = VChoiceField( label="Samples grouping", widget=forms.Select( attrs={ 'onchange': 'settagvalues();' } ) )
	numbest = forms.ChoiceField( label = "Number of best species / units", choices = [ ( "10", "10" ), ( "25", "25"), ( "50", "50" ), ( "100", "100" ) ] )
	samples = VMChoiceField( widget=forms.SelectMultiple, choices = [] )
	
class Whittaker( GenericForm ):
	level = VChoiceField()
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dgroup = VChoiceField( label="Samples grouping" )
	xscale = forms.ChoiceField( label = "Ranks scale", choices = [ ( v, v ) for v in [ "linear", "logarithmic", "sqr-log" ] ] )
	yscale = forms.ChoiceField( label = "Abundance scale", choices = [ ( v, v ) for v in [ "linear", "logarithmic" ] ] )
	dmarks = forms.ChoiceField( label = "Groups separation", choices = [ ( v, v ) for v in [ "none", "color", "shape", "both" ] ] )
	regression = forms.ChoiceField( label = "Include regression line", choices = [ ( v, v ) for v in [ "no", "yes", "2 parts" ] ] )

class Venn( GenericForm ):
	level = VChoiceField()
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dgroup = VChoiceField( label="Samples grouping", widget=forms.Select( attrs={ 'onchange': 'settagvalues();' } ) )
	ptype = forms.ChoiceField( label = "Presentation", choices = [ ( v, v ) for v in [ "proportional-presence", "proportional-abundance", "interactive-presence" ] ] )
	samples = VMChoiceField( widget=forms.SelectMultiple, choices = [] )
	
class PCA2P( GenericForm ):
	level = VChoiceField()
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	color = VChoiceField( label="Color" )
	shape = VChoiceField( label="Shape" )
	labels = VChoiceField( label="Labels" )
	pc1 = forms.ChoiceField( label = "PC1", choices = [ ( str( v ), str( v ) ) for v in range( 1, 5 ) ] )
	pc2 = forms.ChoiceField( label = "PC2", choices = [ ( str( v ), str( v ) ) for v in range( 2, 6 ) ] )

						 