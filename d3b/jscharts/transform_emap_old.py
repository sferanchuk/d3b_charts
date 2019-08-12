#!/usr/bin/python2

import sys
import math
from scipy.special import erf

iname = sys.argv[1]
oname = sys.argv[2]
ttype = sys.argv[3]

with open( iname ) as f:
	data = f.readlines()
	
ls0 = data[0].split()
for nfcol in range( 1, len( ls0 ) ):
	if ls0[ nfcol ] != str( nfcol + 1 ):
		break

print(nfcol)

tstype = ttype.split()
if True:	
	maxdata = 0.
	mindata = 0.
	sumdata = 0.
	sumsq = 0.
	cnt = 0
	for lc in range( 1, len( data ) ):
		ls = data[lc].strip().split( "\t" )
		fls = list(map( float, ls[ nfcol : ] )) 
		maxdata = max( maxdata, max( fls ) )
		mindata = min( mindata, min( fls ) )
		sumdata += sum( fls )
		sumsq += sum( [ v * v for v in fls ] )
		cnt += len( fls )
	scale = float( tstype[1] ) 
	mean = sumdata / cnt
	mdisp = math.sqrt( sumsq / cnt - mean * mean )
	with open( oname, "wt" ) as f:
		f.write( data[0] )
		for lc in range( 1, len( data ) ):
			ls = data[lc].strip().split( "\t" )
			fls = list(map( float, ls[ nfcol : ] )) 
			f.write( "\t".join( ls[ : nfcol ] ) + "\t" )
			if tstype[0] == "scale":
				tls = [ str( int( ( v - mindata ) * scale / ( maxdata - mindata ) ) ) for v in fls ]
			elif tstype[0] == "erf":
				tls = [ str( int( scale * ( 1. + erf( ( v - mean ) / mdisp ) ) ) ) for v in fls ]
			f.write( "\t".join( tls ) + "\t\n" )
	



