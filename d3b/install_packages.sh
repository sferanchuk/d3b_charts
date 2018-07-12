#!/bin/bash
#pip2 install virtualenv
#virtualenv d3b_py2
source d3b_py2/bin/activate
pip install "django<2"
#pip install numpy
#pip install scipy
#pip install scikit-learn
pip install "scikit-bio<0.5"
pip install biom-format
pip install h5py
pip install wsgi
deactivate
#wget https://d3js.org/d3.v3.min.js
#wget https://d3js.org/d3.v4.min.js
#wget https://github.com/benfred/venn.js/blob/master/venn.js