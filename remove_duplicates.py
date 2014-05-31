#!/opt/local/bin/ipython

# BLAST UCSC
#http://genome.ucsc.edu/cgi-bin/hgBlat

from pysequence import read_sequences, write_sequences
from collections import OrderedDict

code_description = """
-------------------------------------
Playing with Genome
-------------------------------------

Remove identical sequences.

"""

remove_exact_duplicates = True

if remove_exact_duplicates:

	infile = raw_input("type the name of the file to be processed:\n")
	outfile = raw_input("type the name of the file where to save the sequences:\n")
	keep_record = raw_input("do you want to keep track of the identical sequences? Much slower (T/F)\n")

	if keep_record.upper() == 'T':
		keep_record = True
	else:
		keep_record = False

	print ' Removing exact duplicates:---------------------------------'
	print " infile = %s" %infile
	print " outfile = %s" %outfile
	print " keep_record = %r" %keep_record

	[ seq_names, seq_values] = read_sequences( infile )

# The fastest smarter way of doing this is what is done at the end, using dictionaries
# However, this way the correspondence is lost, dict assigns the last entry to the seq
# I want to keep a record.

	check_dict = OrderedDict( zip(seq_names, seq_values) )
	#nmax = 1000
	values = seq_values#[0:nmax]
	names = seq_names#[0:nmax]
	Nseq = len( seq_values )
	
# ------ Removing exact matches using dictionary
	if keep_record:
		seq_dict = OrderedDict()
		cnt = 0
		for i,entry in enumerate(values):
			#print i, entry
			if entry not in seq_dict.keys():
				cnt += 1
				#print cnt
				seq_dict[entry] = [ names[i] ]
			else:
				str = seq_dict[entry]
				str.append( names[i] )
				seq_dict[entry] = str

		print len( seq_dict )

		print 'Single entry fraction: %f' %( float(cnt)/len(values) )

		subseq_values = []
		subseq_names = []

	# To save look-up table of the matches
		discard_names = []

		for key, value in seq_dict.iteritems():
			#print key, value, len(value)
			subseq_values.append( key )
			subseq_names.append( value[0] )
			if  len( value ) > 1:
				discard_names.append( value[0]+':'+','.join(value[1:]) )

		print '# exact duplicates = %i' %( len(check_dict)-len(seq_dict) )

		write_sequences(outfile, subseq_names, subseq_values)

# ------ Saving sequences matches
		outfile = outfile+'.single2matches'
		body = '\n'.join( discard_names )
		output = open( outfile, 'w')
		output.write(body)
		output.close()
		
	else:
# Smart fastest way of selecting single entries
		single_dict = OrderedDict( zip(seq_values, seq_names) )
		print len(single_dict)

		subseq_names = single_dict.values()
		subseq_values = single_dict.keys()
		write_sequences(outfile, subseq_names, subseq_values)

	print 'Done.'
		

				