
from setuptools import setup
import requests
import platform
import os

setup(name='d3b_charts',
      version='0.7',
      description='Data presentation tools for microbial ecology',
      author='Sergey Feranchuk',
      author_email='feranchuk@gmail.com',
      url='https://github.com/sferanchuk/d3b_charts/',
      install_requires=[ "django<2", "numpy", "scipy", "scikit-learn", "scikit-bio<0.5" ]
     )
	  
	  
