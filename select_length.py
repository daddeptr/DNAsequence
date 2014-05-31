#!/opt/local/bin/ipython

import numpy as np
from collections import OrderedDict
from pysequence import read_sequences, write_sequences

# BLAST UCSC
#http://genome.ucsc.edu/cgi-bin/hgBlat

select_length = True
# ------ Select sequences based on length
if select_length:
#if select_40bases_plus:
	
	infile = raw_input("type the name of the file to be processed:\n")
	outfile = raw_input("type the name of the file where to save the sequences:\n\
	Note that only single entries are kept.\n")
	length = int( raw_input("type the desired length (sequences atr trimmed to L):\n") )

	print ' Trimming :------------------------------------------------'
	print " infile = %s" %infile
	print " outfile = %s" %outfile
	print " length = %4i" %length

	[ seq_names, seq_values ] = read_sequences( infile )

	print len(seq_names)
	print len(seq_values)

	seq_dict = OrderedDict( zip(seq_names, seq_values) )

	print len(seq_dict)
	
	selection = seq_values
	Nseq = len( selection )
	seq_length = np.zeros( Nseq, dtype=int )
	for i,seq in enumerate( selection ):
		seq_length[i] = len( seq )

	print 'Minimum length = %i' %min( seq_length )
	print 'Maximum length = %i' %max( seq_length )

	sort_seq_length = np.sort( seq_length )
	isort = np.argsort( seq_length )	
	sort_seq_names = []
	sort_seq_values = []

	ok = 0
	err = 0
	for i in range(Nseq):
		sort_seq_names.append( seq_names[isort[i]] )
		sort_seq_values.append( seq_values[isort[i]] )
		if seq_length[isort[i]] == sort_seq_length[i]:
			ok += 1
		else:
			err += 1
			print 'Sorting error'
	print ok, err
		#print ''.join( ['-']*75 )
	
	N_seq = np.sum( sort_seq_length >= (np.ones( Nseq, dtype=int ) * length ) )
	Nshort_seq = np.sum( sort_seq_length < (np.ones( Nseq, dtype=int ) * length ) )
#	N20short_seq = np.sum( sort_seq_length < (np.ones( Nseq, dtype=int ) * 20 ) )
	print 'Total sequences number: %6i' %( Nseq )
	print 'Kept sequences number: %6i' %( N_seq )
	print 'Short sequences number: %6i' %( Nshort_seq )
#	print '20-Short sequences number: %6i' %( N20short_seq )
	print 'Short sequences fraction: %f' %( float( Nshort_seq ) / Nseq )
	print 'Retained sequences fraction: %f' %( 1. - float( Nshort_seq ) / Nseq )

	subseq_values = []
	subseq_names  = []

	bkp_subseq_values = []
	bkp_subseq_names  = []

#	short_values = []
#	short_names  = []

## ------ 0 to 20
#	for i, seq in enumerate( sort_seq_values[ 0:N20short_seq ] ):
#		short_values.append( seq )
#		short_names.append( sort_seq_names[i] )
#		if Feedback and False:
#			print 'Long seqs: %s %6i, %4i, %s' %( short_names[-1], i,len(seq),seq )

# ------ 20 to 40
	for i,seq in enumerate( sort_seq_values[ 0:Nshort_seq ] ):
		bkp_subseq_values.append( seq )
		bkp_subseq_names.append( sort_seq_names[i] )
#		if Feedback and False:
#			print 'Long seqs: %s %6i, %4i, %s' %( bkp_subseq_names[-1], i,len(seq),seq )

# ------ Longer than 40 trimmed to 40
	for i,seq in enumerate( sort_seq_values[ Nshort_seq: ] ):
		subseq_values.append( seq[0:length] )
		#subseq_values.append( seq ) # whole sequence kept to help similar seqs search
		subseq_names.append( sort_seq_names[Nshort_seq+i] )
		
#	print 'Num 0-20 sequences = %i' %len(short_values)
	print 'Num Short sequences = %i' %len(bkp_subseq_values)
	print 'Num long+ sequences = %i' %len(subseq_values)
	print 'Total: %i' %( len(bkp_subseq_values) + len(subseq_values) )
	
	write_sequences( outfile, subseq_names, subseq_values )

	if len( bkp_subseq_names )>0:
		outfile = outfile+".discarded"
		write_sequences( outfile, bkp_subseq_names, bkp_subseq_values )

	print 'Done.'