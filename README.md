# d3b_charts

Data presentation tools for microbial ecology.

## Description

The online services available on this site have been developed for the extensive analysis of tables with estimates of bacteria abundance levels in environmental samples. In most cases the microbial ecology-specific functionality is implemented by the scikit-bio Python package, together with the other Python packages intended for analysis of big datasets. The interactive visualisation tools are implemented by the D3.js software library; therefore, the software project is named D3b. The source codes are available at github as an installable package.

## Installation

The installation process may require the skills of a system administrator, if there are incompatibilities between the python3 environment and the software dependencies.

Suggested installation via conda:

1. create virtual environment:

`conda create -n d3b python=3.6 anaconda`

`source activate d3b`

2. run setup:

`git clone https://github.com/sferanchuk/d3b_charts.git`
`cd d3b_charts/d3b`
`pip install .`

3. download phantomjs package, extract the executable file from the archive and move it to 'phantomjs' subdirectory. Suggested URLs to obtain phantomjs:

Windows: https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-windows.zip

Linux x86_64: https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-linux-x86_64.tar.bz2

Linux i686: https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-linux-i686.tar.bz2

Mac: https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-macosx.zip


This step is optional; phantomjs is used to render graphics on a server side and is required to export PNG and SVG images. Other packages external to python environment, like d3.js, are already included into the distribution.

4. Run django server (runserver.sh) in console. Point the browser to URL 'http://localhost:8000'. Modify ALLOWED_HOSTS variable in d3b/settings.py to resolve the problems with a remote access to pages.

## Getting started

Use links to the sample datasets on a start page to explore the abilities of the d3b system. On each link, the menu provides tools to analyse a dataset in various representations.

For each tool, the result becomes available in the browser after pressing the 'Submit' button. The input parameters should be specified before submitting the task. In order to become familiar with the design and conventions of the d3b system, it is advisable to first try out the default settings on the sample datasets, 

For this purpose, datasets composed of abundances of annotated OTUs in microbial communities are provided. Four of the sample datasets are estimates of amplicon sequencing followed by the QIIME1 pipeline, and one dataset represents a selection from a custom ecological community.

All the datasets are represented as matrices of abundance counts. To allow the specification of traits and taxonomic groups, the descriptors for the rows and columns of the matrix can be modified and saved. The tools which allow to specify both kinds of descriptors (Subsets of samples / Taxonomic groups) are also available in the sidebar menu for any dataset.

## Usage

To apply the tools from the d3b system to custom data, a dataset in the form of a matrix should be loaded to the server side of the system. After that, a permanent link to the processed dataset becomes available. It should be kept and can be accessed at any time when the whole system is up and running.

Two formats of input matrix are accepted. The first is tab-delimited format (TSV) and the second is BIOM format. In both formats, abundance counts for OTU are combined with taxonomic annotation of the OTU and identifiers of the samples.

A tab-delimited file exported from a spreadsheet using packages like MS Excel or LibreOfiice Calc is suitable for the system, provided the conventions for the contents of rows and columns are followed. And the processed dataset can be exported back to tab-delimited format, as provided for in the tool 'matrix of abundances'.

The output obtained from that tool applied to one of the sample datasets can provide an example of a correctly composed tab-delimited file. In that file, each line begins from a complete taxonomic annotation for the OTU, in a fixed number of cells, followed by absolute counts of abundances for this OTU in all of the samples. In the first line, cells with numbers from 1 to 2, 3, etc., are used to specify fixed levels of hierarchy in the taxonomic system, and they are followed by cells with identifiers of the samples.

Data in BIOM format ("biological observation matrix") is also suitable for import into d3b. This format is often used to keep intermediate output results in software systems for microbiologists. But due to the flexibility of the BIOM format specification, the import from a particular BIOM-formatted file may fail, or may give wrong results. In this case, the python script 'biom2emap.py' from the source codes of d3b system can be used as a template for a user-defined script, to convert any biom-formatted dataset to a correct tab-delimited file.

## Overview of tools

The tabular presentations of the data include a table of abundances, just as it is loaded into the system. There are various options for sorting the rows, merging the columns, reducing the level of taxonomic hierarchy, etc.. Most of the input forms include a choice of level of taxonomic hierarchy, the possibility to restrict the analysis or data presentation to certain taxonomic groups, and the possibility to merge samples into pre-defined groups.

The tabular presentations also include: 
1. the values of alpha-diversity calculated using several of the most informative estimators, 
2. the significance of differences for alpha-diversity values between several groups
3. the significance of differences between groups of samples calculated using the distances between samples, with several alternative measures of distance.

The graphical presentations, implemented with the use of the d3.js library, currenly include the following charts: 
1. A bubble chart and heatmap, to represent absolute/relative abundances. 
2. 2D scatter charts, to represent the results of several data ordination methods, such as principal component analysis (PCA) or multi-dimensional scaling (MDS). The choice of several measures is available here to calculate distances between the samples. 
3. A dendrogram (tree) to represent the degree of proximity between samples. 
4. A Venn diagram to represent the unique and shared taxonomic units for the samples, implemented with the use of the jVenn plugin and the venn.js library. 
5. Two kinds of diagrams to present distributions which describe a sample or a group of samples: a rank-abundance chart (Whittaker plot) to represent the distribution of relative species abundance, and a rarefaction curve to estimate the effect of insufficient coverage and sample size. 
6. A ternary chart, to represent the relative abundances of bacterial phylotypes for the three samples, or groups of samples. 
7. A Volcano chart and mean-distance plot, to represent distribution of abundances and differentiation between traits. 
8. Two combined 2D charts, to represent the results of PCA decomposition applied directly to a non-square matrix of abundances. One chart is for the samples in the survey and the second adjacent chart is for the bacterial species in the rows of the submitted matrix.










