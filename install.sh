#!/bin/bash
pip2 install virtualenv
virtualenv d3b_py2
source d3b_py2/bin/activate
pip install .
deactivate
