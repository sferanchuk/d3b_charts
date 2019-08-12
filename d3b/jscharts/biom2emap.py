#!/usr/bin/python2
import biom
import sys
if len( sys.argv ) == 1:
	print("arguments: biom_in")
	sys.exit(1)
	
try:
	t = biom.load_table( sys.argv[1] )
	samples = t.ids()
	spec = t.ids( axis = 'observation' )
	maxlevel = 0
	for i in range( len( spec ) ):
		tax = t.metadata( spec[i], 'observation' )[ 'taxonomy' ]
		maxlevel = max( maxlevel, len( tax ) )
	s = "\t".join( str( i + 1 ) for i in range( maxlevel + 1 ) ) + "\t" + "".join( v + "\t" for v in samples ) 
	print(s)
	for i in range( len( spec ) ):
		tax = t.metadata( spec[i], 'observation' )[ 'taxonomy' ]
		cs = ""
		for k in range( maxlevel ):
			cst = ""
			if k < len( tax ) and len( tax[k] ) > 3:
				if tax[k][2] == "_":
					cst = tax[k][3:] 
				else:
					cst = tax[k]
			cs += cst + "\t"
		cs += spec[i] + "\t"
		for kk in range( len( samples ) ):
			cs += str( int( t.get_value_by_ids( spec[i], samples[kk] ) ) ) + "\t"
		print(cs)
except TypeError as err:
	t = biom.load_table( sys.argv[1] )
	samples = t.ids()
	spec = t.ids( axis = 'observation' )
	s = "1\t2\t3\t" + "".join( v + "\t" for v in samples ) 
	print(s)
	for j in range( len( spec ) ):
		#cs = "%s\t" % spec[j]
		tax = t.metadata( spec[j], 'observation' )[ 'KEGG_Pathways' ]
		cs = ""
		for k in range( 3 ):
			if k < len( tax ):
				cs += "%s\t" % tax[k]
			else:
				cs += "*\t"
		for i in range( len( samples ) ):
			cs += str( int( t.get_value_by_ids( spec[j], samples[i] ) ) ) + "\t"
		print(cs)
	
