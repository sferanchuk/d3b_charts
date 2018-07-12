source d3b_py2/bin/activate
#python2 manage.py runserver 0:8000
uwsgi --http :8000 --wsgi-file d3b/wsgi.py --async 10 --ugreen