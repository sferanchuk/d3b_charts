
import json
from django import forms
import controls

class UploadFile( forms.Form ):
	name = forms.CharField(max_length=50, label="Dataset name")
	#transform = forms.ChoiceField( label="Transform", choices = [ ( "none", "none" ), ( "scale 1000", "scale from 0 to 1000" ), ( "erf 1000", "erf 1000" ) ], initial = [ "none" ] ) 
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
			clevels = [ ( str(k), levelnames[ summary[ "taxtype" ] ][ k - 1 ] + " (" + str( k ) + ")" ) for k in range( 1, summary[ 'maxlevel' ] + 1 ) ]
		else:
			clevels = [ ( str(k), str(k) ) for k in range( 1, summary[ 'maxlevel' ] + 1 ) ]
		ctags = [ ( v, v ) for v in list(tags.keys()) ]
		ctaxfilters = [ ( v, v ) for v in list(taxfilters.keys()) ]
		distmeasures = [ "Pearson", "Kendall", "Spearman", "Euclidean", "Bray-Curtis", "Jaccard", "Morisita-Horn", "Unifrac-unweighted", "Unifrac-weighted" ]
		cdmeasures = [ ( v, v ) for v in distmeasures ]
		
		if 'level' in self.fields:
			self.fields[ 'level' ].widget.choices += clevels 
			self.initial[ 'level' ] = [ clevels[ min( summary[ 'maxlevel' ] - 1, 3 ) ][0] ]
			#self.fields[ 'level' ].widget.initial = [ clevels[ min( summary[ 'maxlevel' ], 3 ) ][1] ]
		for fieldname in [ 'dorder', 'volume' ]:
			if fieldname in self.fields:
				self.fields[ fieldname ].widget.choices += cvolumes
		for fieldname in [ 'dmethod' ]:
			if fieldname in self.fields:
				self.fields[ fieldname ].widget.choices += cdmeasures
		for fieldname in [ 'dfilter', 'color', 'shape', 'dgroup', 'order1', 'order2', 'labels', 'shisttag' ]:
			if fieldname in self.fields:
				self.fields[ fieldname ].widget.choices += ctags
				if fieldname in [ 'labels' ]:
					self.initial[ fieldname ] = [ 'name' ]
				else:
					self.initial[ fieldname ] = [ 'none' ]
		for fieldname in [ 'spfilter', 'spshow' ]:
			if fieldname in self.fields:
				self.initial[ fieldname ] = [ 'none' ]
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
			print(snames)
			self.fields[ 'samples' ].widget.choices = [ ( v, v ) for v in snames ]
			self.initial[ 'dgroup' ] = [ 'name' ]

class Table( GenericForm ):
	level = VChoiceField()
	dnorm = forms.ChoiceField( label = "Units", choices = [ ( "count", "Counts" ), ( "percent", "Percents" ) ] )
	dorder = VChoiceField( label="Sort order", choices = [ ( "none", "Unchanged" ), ( "taxonomy", "Taxonomy" ) ] )
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dgroup = VChoiceField( label="Samples grouping" )
	numbest = forms.ChoiceField( label = "Number of best species / units", choices = [ ( "all", "all" ), ( "10", "10" ), ( "15", "15" ), ( "25", "25"), ( "50", "50" ), ( "100", "100" ) ] )

class Indices( GenericForm ):
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dgroup = VChoiceField( label="Samples grouping" )
	mmethod = forms.ChoiceField( label="Subset of estimators", choices = [ ( "all", "Major subset" ), ( "mm-fit", "Rarefaction" ) ] )
	level = VChoiceField( choices = [ ( "all", "All" ) ] )
	
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
	permutations = forms.ChoiceField( label = "Permutations", choices = [ ( v, v ) for v in [ "999", "3999", "9999", "39999" ] ] )

class PCA( GenericForm ):
	level = VChoiceField()
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dmethod = VChoiceField( label="Distances" )
	fmethod = forms.ChoiceField( label = "Method", choices = [ ( v, v ) for v in [ "CA (skbio)", "PCoA (skbio)", "PCA (sklearn)", "MDS (sklearn)" ] ] )
	color = VChoiceField( label="Color" )
	shape = VChoiceField( label="Shape" )
	labels = VChoiceField( label="Labels" )
	dsizes = forms.ChoiceField( label = "Size of points", choices = [ ( v, v ) for v in [ "equal", "linear", "-(1/2)", "-(1/4)" ] ] )

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
	dtype = forms.ChoiceField( label = "Units", choices = [ ( "count", "Counts" ), ( "percent", "Percents" ), ( "z-score", "Z-Score" ), ( "percent-quantile", "Percents-Quantile" ), ( "raw", "Raw" ) ] )
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dgroup = VChoiceField( label="- Samples grouping - excludes sorting (*)" )
	order1 = VChoiceField( label="* Sort order (primary)" )
	order2 = VChoiceField( label="* Sort order (secondary)" )
	labels = VChoiceField( label="* Labels for samples" )
	numbest = forms.ChoiceField( label = "Number of best species / units", choices = [ ( "all", "all" ), ( "10", "10" ), ( "15", "15" ), ( "25", "25"), ( "50", "50" ), ( "100", "100" ) ] )
	percentbest = forms.ChoiceField( label = "Minimal average abundance", choices = [ ( "all", "all" ), ( "1%", "1%" ), ( "2%", "2%" ), ( "3%", "3%" ), ( "5%", "5%"), ( "10%", "10%" ) ] )
	shist = forms.ChoiceField( label = "Histogram of variations", choices = [ ( "none", "none" ), ( "best-fisher", "Exact Fisher (best pair)" ), ( "anova", "Anova" ), ( "best-wilcoxon", "Ranked (best pair)" ), ( "best-ttest", "T-Test (best pair)" ), ( "best-chisquare", "Chi-square (best pair)" ), ( "best-logfc", "Fold Count (best pair)" ) ] )
	shisttag = VChoiceField( label = "Groups for histogram" )
	spshow = VChoiceField( label="Show taxonomy", choices = [ ( "custom", "Custom" ) ] )
	cscale = forms.ChoiceField( label = "Color scale", choices = [ ( "threshold", "threshold" ), ( "linear", "linear" ) ] ) 
	cpalette = forms.ChoiceField( label = "Color palette", choices = [ ( "pink/brown", "pink/brown" ),  ( "blue/green", "blue/green" ), ( "red/black/green", "red/black/green (two-sided)" ), ( "red/white/blue", "red/white/blue (two-sided)" ) ] ) 
	numhighlight = forms.ChoiceField( label = "Highlight best species/units", choices = [ ( v, v ) for v in [ "all", "20%", "10%", "7.5%", "5%", "2%" ] ] )
	taxlabels = forms.ChoiceField( label = "Labels for taxonomy", choices = [ ( "no", "No" ), ( "yes", "Yes" ) ] )
	spcustom = forms.CharField( widget=forms.HiddenInput() )
	
class BubbleChart( GenericForm ):
	level = VChoiceField()
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dgroup = VChoiceField( label="Samples grouping" )
	dnorm = forms.ChoiceField( label = "Units", choices = [ ( "percent", "Percents" ), ( "count", "Counts" ), ( "percent-quantile", "Percents-Quantile" ) ] )
	dglabels = forms.ChoiceField( label = "Group labels", choices = [ ( "no", "No" ), ( "yes", "Yes" ) ] )
	dtlabels = forms.ChoiceField( label = "Percent labels", choices = [ ( "no", "No" ), ( "yes", "Yes" ) ] )
	dlegend = forms.ChoiceField( label = "Show legend", choices = [ ( "no", "No" ), ( "yes", "Yes" ) ] )
	dprestype = forms.ChoiceField( label = "Presentation type", choices = [ ( "bubble", "Bubble chart" ), ( "treemap", "Treemap" ), ( "bars", "Bars" ) ] )

class Ternary( GenericForm ):
	level = VChoiceField()
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dgroup = VChoiceField( label="Samples grouping", widget=forms.Select( attrs={ 'onchange': 'settagvalues();' } ) )
	numbest = forms.ChoiceField( label = "Number of best species / units", choices = [ ( "10", "10" ), ( "15", "15" ), ( "25", "25"), ( "50", "50" ), ( "100", "100" ) ] )
	samples = VMChoiceField( widget=forms.SelectMultiple, choices = [] )
	
class Whittaker( GenericForm ):
	level = VChoiceField()
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dgroup = VChoiceField( label="Samples grouping" )
	#xscale = forms.ChoiceField( label = "Ranks scale", choices = [ ( v, v ) for v in [ "linear", "logarithmic", "sqr-log" ] ] )
	#yscale = forms.ChoiceField( label = "Abundance scale", choices = [ ( v, v ) for v in [ "linear", "logarithmic" ] ] )
	ptype = forms.ChoiceField( label = "Presentation", choices = [ ( v, v ) for v in [ "log-log", "log-sqrlog", "log-linear", "lorentz", "rarefaction" ] ] ) 
	dmarks = forms.ChoiceField( label = "Groups separation", choices = [ ( v, v ) for v in [ "none", "color", "shape", "both" ] ] )
	color = VChoiceField( label="Color" )
	shape = VChoiceField( label="Shape" )
	regression = forms.ChoiceField( label = "Include regression line", choices = [ ( v, v ) for v in [ "no", "yes", "2 parts" ] ] )

class Volcano( GenericForm ):
	level = VChoiceField()
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	dgroup = VChoiceField( label="Samples grouping", widget=forms.Select( attrs={ 'onchange': 'settagvalues();' } ) )
	samples = VMChoiceField( widget=forms.SelectMultiple, choices = [] )
	ptype = forms.ChoiceField( label = "Presentation", choices = [ ( "volcano", "Volcano chart" ), ( "mdplot", "Mean-Difference chart" ) ] ) 
	cmethod = forms.ChoiceField( label = "Criterion", choices = [ ( "anova", "Anova" ), ( "wilcoxon", "Ranked" ), ( "ttest", "T-Test" ) ] )

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
	dsizes = forms.ChoiceField( label = "Size of points for samples", choices = [ ( v, v ) for v in [ "equal", "linear", "-(1/2)", "-(1/4)" ] ] )
	spsizes = forms.ChoiceField( label = "Size of points for species", choices = [ ( v, v ) for v in [ "equal", "linear", "-(1/2)", "-(1/4)" ] ] )
	pc1 = forms.ChoiceField( label = "PC1", choices = [ ( str( v ), str( v ) ) for v in range( 1, 5 ) ] )
	pc2 = forms.ChoiceField( label = "PC2", choices = [ ( str( v ), str( v ) ) for v in range( 2, 6 ) ] )
	axis = forms.ChoiceField( label = "Axis to scale", choices = [ ( str( v ), str( v ) ) for v in [ "default", "1" ] ] )
	
class Network( GenericForm ):
	level = VChoiceField()
	dfilter = VChoiceField( label="Samples filter" )
	spfilter = VChoiceField( label="Taxonomy filter" )
	numbest = forms.ChoiceField( label = "Number of best species / units", choices = [ ( "all", "all" ), ( "10", "10" ), ( "15", "15" ), ( "25", "25"), ( "50", "50" ), ( "100", "100" ) ] )
	cmethod =  forms.ChoiceField( label = "Measure of correlation", choices = [ ( str( v ), str( v ) ) for v in [ "Pearson", "Spearman", "Kendall", "Binary" ] ] )
	cthreshold = forms.ChoiceField( label = "Threshold of correlation", choices = [ ( v, v ) for v in [ "0.5", "0.7", "0.9" ] ] )
	psizes = forms.ChoiceField( label = "Size of points", choices = [ ( v, v ) for v in [ "equal", "propotional" ] ] )

	

						 
