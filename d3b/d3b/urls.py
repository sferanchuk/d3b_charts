"""d3b URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.11/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import url
from django.contrib import admin
#import sys
#sys.path.insert(0, './d3b')
import views

urlpatterns = [
	
	url( r'^$', views.new, name='new' ),
	url( r'^summary/(?P<job>[0-9a-f]+)/$', views.summary, name='summary' ),
	url( r'^tags/(?P<job>[0-9a-f]+)/$', views.tags, name='tags' ),
	url( r'^taxonomy/(?P<job>[0-9a-f]+)/$', views.taxonomy, name='taxonomy' ),
	url( r'^table/(?P<job>[0-9a-f]+)/$', views.table, name='table' ),
	url( r'^indices/(?P<job>[0-9a-f]+)/$', views.indices, name='indices' ),
	url( r'^anova/(?P<job>[0-9a-f]+)/$', views.anova, name='anova' ),
	url( r'^permanova/(?P<job>[0-9a-f]+)/$', views.permanova, name='permanova' ),
	url( r'^bubbles/(?P<job>[0-9a-f]+)/$', views.bubbles, name='bubbles' ),
	url( r'^heatmap/(?P<job>[0-9a-f]+)/$', views.heatmap, name='heatmap' ),
	url( r'^tree/(?P<job>[0-9a-f]+)/$', views.tree, name='tree' ),
	url( r'^pca/(?P<job>[0-9a-f]+)/$', views.pca, name='pca' ),
	url( r'^venn/(?P<job>[0-9a-f]+)/$', views.venn, name='venn' ),
	url( r'^ternary/(?P<job>[0-9a-f]+)/$', views.ternary, name='ternary' ),
	url( r'^whittaker/(?P<job>[0-9a-f]+)/$', views.whittaker, name='whittaker' ),
	url( r'^volcano/(?P<job>[0-9a-f]+)/$', views.volcano, name='volcano' ),
	url( r'^pca_sp/(?P<job>[0-9a-f]+)/$', views.pca_sp, name='pca_sp' ),
]

