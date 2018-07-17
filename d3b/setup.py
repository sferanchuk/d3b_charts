
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
	  
	  
def download_file( url, path ):
	if os.path.isfile( path ):
		return
	r = requests.get(settings.STATICMAP_URL.format(**data), stream=True)
	if r.status_code == 200:
		with open(path, 'wb') as f:
			for chunk in r.iter_content(1024):
				f.write(chunk)
				
download_file( "https://d3js.org/d3.v3.min.js", "d3b/static/d3.v3.min.js" )
download_file( "https://d3js.org/d3.v4.min.js", "d3b/static/d3.v4.min.js" )
download_file( "https://github.com/benfred/venn.js/blob/master/venn.js", "d3b/static/venn.js" )
				
				
def download_phantomjs:
	if platform.system() == "Windows":
		download_file( "https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-windows.zip", "phantomjs-2.1.1.zip" )
	elif platform.system() == "Linux":
		if platform.machine() == "x86_64":
			download_file( "https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-linux-x86_64.tar.bz2", "phantomjs-2.1.1.tar.bz2" )
		elif platform.machine() == "i686":
			download_file( "https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-linux-i686.tar.bz2", "phantomjs-2.1.1.tar.bz2" )
	elif platform.system() == "Darwin":
		download_file( "https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-macosx.zip", "phantomjs-2.1.1.zip" )


			
