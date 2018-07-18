#!/bin/bash

pip2 install virtualenv
virtualenv d3b_py2
source d3b_py2/bin/activate
pip install requests
pip install .
pip install biom-format
pip install h5py

deactivate

for file in d3.v3.min.js d3.v4.min.js
do
	path=d3b/static/$file
	if [ ! -f $path ]
	then
		wget -O $path "https://d3js.org/$file" 
	fi
done

if [ ! -f d3b/static/venn.js ]
then
	wget -O d3b/static/venn.js "https://raw.githubusercontent.com/benfred/venn.js/master/venn.js"
fi

fpath=phantomjs/phantomjs
if [ ! -f $fpath ]
then
	fdname="phantomjs-2.1.1-linux-x86_64"
	faname="${fdname}.tar.bz2"
	if [ ! -f $faname ]
	then
		url="https://bitbucket.org/ariya/phantomjs/downloads/$faname"
		wget $url
	fi
	tar jxvf $faname ${fdname}/bin/phantomjs
	mv $fdname/bin/phantomjs $fpath
fi
