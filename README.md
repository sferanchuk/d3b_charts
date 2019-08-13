# d3b_charts

Data presentation tools for microbial ecology.

## Description

The software is implemented as a web appication based on Django framework.
The format of input data is compatible with Biological Observation Matrix (BIOM) format.
The interactive vizulalization tools are implemented with the use of D3.js software library.

Most of the functionality specific to microbial ecology is implemented with the use of scikit-bio python package, 
in the conjunction with other python packages intended to be used in big data analysis.

## Installation

An installation process may require skills of system administrator, in a case of incompatibilities in the python3 environment and the software dependencies.

Suggested installation via conda:

1. create virtual environment:

conda create -n d3b python=3.6 anaconda

source activate d3b

2. run setup.py

3. download phantomjs package, extract the executable file from the archive and move it to 'phantomjs' subdirectory. Suggested URLs to obtain phantomjs:

Windows: https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-windows.zip

Linux x86_64: https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-linux-x86_64.tar.bz2

Linux i686: https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-linux-i686.tar.bz2

Mac: https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-macosx.zip


This step is optional; phantomjs is used to render graphics on a server side and is required to export PNG and SVG images. Other packages external to python environment, like d3.js, are already included into the distribution.

4. Run django server (runserver.sh) in console. Point the browser to URL 'http://localhost:8000'. Modify ALLOWED_HOSTS variable in d3b/settings.py to resolve the problems with a remote access to pages.

## Getting started

Use two sample datasets included into the distribution to explore the abilities of d3b system. The links to the sample datasets are provided in the start page. On each link, the tools listed in the sidebar menu allow to analyse a dataset in various representations. 

For each tool, the result becomes available in the browser after pressing 'Submit' button. The input parameters should be specified before submitting the task. For the sample datasets, the default settings may be often used on first steps to become familiar with the design and conventions of the d3b system.

By an intention, a dataset is composed from abundances of annotated OTUs in microbial comminities. The two sample datasets are estimates by 454 amplicon sequencing followed by QIIME1 pipeline, and by a whole metagenome seqiencing followed by a custom pipeline. 

Anyway, all the datasets are represented as matrices of abundance counts. To allow the specification of traits and taxonomic groups, the descriptors for rows and columns of the matrix could be modified and saved. The tools which allow to specify both kinds of descriptors (Subsets of samples / Taxonomic groups) also available in the sidebar menu for any dataset.

## Usage

To apply the tools from d3b system to a custom data, a dataset in a form of a matrix should be loaded to the server side of the system. After that, a permanent link to the processed dataset becomes available. It should be kept and could be accessed at any time when a whole system is up-running.

The two formats of input matrix are accepted. First is tab-delimited format (TSV) and second is BIOM format. In both formats, abundance counts for an OTU are combined with taxonomic annotation of the OTU and identifiers of samples. 

A tab-delimited file exported from spreadsheet using packages like MS Excel or LibreOfiice Calc is suitable for the system, if the conventions for a content of rows and columns are followed. And the processed dataset could be exported back to tab-delimited format, as it is provided in the tool 'matrix of abundances'. 

The output obtained from that tool applied to one of the sample datasets can provide an example of a correctly composed tab-delimited file. In that file, each line begins from a complete taxonomic annotation for the OTU, in a fixed number of cells, followed by absolute counts of abundances for this OTU in all of the samples. In the first line, cells wuth numbers from 1 to 2, 3, etc., used to specify fixed levels of hierarchy in taxonomic system, and they are followed by cells with identifiers of the samples.

The data in BIOM format ("biological observation matrix") is also suitable for an import into d3b. This format is often used to keep the intermediate output results in software systems for microbiologists. But due to flexibility of BIOM format specification, the import from a paticilar BIOM-formatted file may fail or may give wrong results. In that case, the python script 'biom2emap.py' from the source codes of d3b system could be used as a template for a user-defined script, to convert any biom-formatted dataset to a correct tab-delimited file.











