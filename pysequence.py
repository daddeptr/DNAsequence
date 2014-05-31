#!/opt/local/bin/ipython

#import re # regular expression
#import glob
import os
#import sqlite3 as lite
import sys
#import time
#import datetime
#import csv
import numpy as np
#import scipy
#import matplotlib.pyplot as plt
#import matplotlib.ticker as tckr
#import pylab
#
import Levenshtein
import difflib
from collections import OrderedDict
#import subprocess
#import faast as fst
import kdifference as kd

# BLAST UCSC
#http://genome.ucsc.edu/cgi-bin/hgBlat

# ---------------------------------------------------------------------
def find_sequence_errors( infile, outfile, Verbose=False):
	code = ' >>> FIND_SEQUENCE_ERRORS: '
	if Verbose:
		print '\n'+code+'Checking %s VS %s' %(infile, outfile)

	if Verbose:
		print code+'Reading %s' %(infile)
	[ iseq_names, iseq_values ] = read_sequences( infile, Verbose=True )
			
	if Verbose:
		print code+'Reading %s' %(outfile)
	[ oseq_names, oseq_values ] = read_sequences( outfile, Verbose=True )

#	iseq_dict = OrderedDict( zip(iseq_names, iseq_values) )		
	iseq_dict = OrderedDict( zip(iseq_values, iseq_names) )		
		
#	oseq_dict = OrderedDict( zip(oseq_names, oseq_values) )		
	oseq_dict = OrderedDict( zip(oseq_values, oseq_names) )		
	
	if Verbose:
		print code+" len seq = %i, %i" %( len(iseq_names), len(oseq_names))
		print code+" len dic = %i, %i" %( len(iseq_dict), len(oseq_dict))
	err = 0
	ok = 0

	ikeys = iseq_dict.keys()
	okeys = oseq_dict.keys()

	if Verbose:
		print '# ikeys', len(ikeys), len(iseq_dict), len(iseq_names)
		print '# okeys', len(okeys), len(oseq_dict), len(oseq_names)

	error_list = []

	for i,key in enumerate( okeys ):
		if key in ikeys:
#			full = iseq_dict[key]
#			short = oseq_dict[key]
#			pos = full.find( short )
#			ln = len( short )	
#			if pos >= 0:
			ok += 1
#			else:
#				err += 1
#				print '%i - Error: %s' %(err,key)
#				print ' %s: %s' %( infile, iseq_dict[key] )
#				print ' %s: %s' %( outfile, oseq_dict[key] )
#				print ''.join(['-']*75)
#				error_list.append( key )	

		else:
			err += 1	
			print '%i - Error: %s' %(err,key)
			print '       %s' %( oseq_dict[key] )
			print ''.join(['-']*75)
			error_list.append( key )	
				
	if Verbose:
		print code+'Sequences OK: %i' %ok
		print code+'Sequences ERRORS: %i' %err			
		print code+'Sequences ERRORS: %i' %(err+ok)			

	return error_list

# ---------------------------------------------------------------------
def find_keys_errors( infile, outfile, Verbose=False):
	code = ' >>> FIND_KEYS_ERRORS: '
	if Verbose:
		print '\n'+code+'Checking %s VS %s' %(infile, outfile)

	if Verbose:
		print code+'Reading %s' %(infile)
	[ iseq_names, iseq_values ] = read_sequences( infile, Verbose=True )
			
	if Verbose:
		print code+'Reading %s' %(outfile)
	[ oseq_names, oseq_values ] = read_sequences( outfile, Verbose=True )

	iseq_dict = OrderedDict( zip(iseq_names, iseq_values) )		
		
	oseq_dict = OrderedDict( zip(oseq_names, oseq_values) )		
	
	if Verbose:
		print code+" len seq = %i, %i" %( len(iseq_names), len(oseq_names))
		print code+" len dic = %i, %i" %( len(iseq_dict), len(oseq_dict))
	err = 0
	ok = 0

	ikeys = iseq_dict.keys()
	okeys = oseq_dict.keys()

	if Verbose:
		print '# ikeys', len(ikeys), len(iseq_dict), len(iseq_names)
		print '# okeys', len(okeys), len(oseq_dict), len(oseq_names)

	error_list = []

	for i,key in enumerate( okeys ):
		if key in ikeys:
			full = iseq_dict[key]
			short = oseq_dict[key]
			pos = full.find( short )

			ln = len( short )	
			if pos >= 0:
				ok += 1
			else:
				err += 1
				print '%i - Error: %s' %(err,key)
				print ' %s: %s' %( infile, iseq_dict[key] )
				print ' %s: %s' %( outfile, oseq_dict[key] )
				print ''.join(['-']*75)
				error_list.append( key )	

		else:
			err += 1	
			print '%i - Error: %s' %(err,key)
			print '       %s' %( oseq_dict[key] )
			print ''.join(['-']*75)
			error_list.append( key )	
				
	if Verbose:
		print code+'Sequences OK: %i' %ok
		print code+'Sequences ERRORS: %i' %err			
		print code+'Sequences ERRORS: %i' %(err+ok)			

	return error_list

# ---------------------------------------------------------------------
def get_match_quality( string1, string2, Verbose=False, marker='v', is_linker=False, tol=0.2 ):
	code = ' >>> GET_MATCH_QUALITY: '
# Matching first string onto the second one. Quality established based on the minimum
# between the length of the string 1 and string2
	if Verbose:
		print code+'Comparing: \n 1)- %s \n 2)- %s' %(string1, string2)
		#print '                        using difflib (could be sped up by using Levenshtein)'

	len1 = len( string1 )
	len2 = len( string2 )
	# create matching object: operations to be applied to string1 to become string2
	# The last bit is the sequence we are interested in
	match = difflib.SequenceMatcher( None, string1, string2 )
	# get matching blocks to select last piece to be removed
	common_blocks = match.get_matching_blocks()
	if Verbose:
		print 'Common blocks: ', common_blocks
	if not is_linker:
		last_match = common_blocks[-2]
	else:
		last_match = common_blocks[0]
	
	if Verbose:
		print 'Last match: ', last_match

	str_match = []
	virus = []
	opcode = match.get_opcodes()
	#print opcode
	delta = 0
	Nequal = 0
	Ninsert = 0
	Ndelete = 0
	Nreplace = 0
	comparison_len = 0
	#Reconstructing string and quantifying the mismatch
	string1_pos = len( string2 )
	for tag, i1, i2, j1, j2 in opcode:
		if Verbose:
			print ("%7s a[%d:%d] (%s) b[%d:%d] (%s)" %(tag, i1, i2, string1[i1:i2], j1, j2, string2[j1:j2]))
		if tag == 'equal':
			# Identifying matching (virus) characters
			patch = [ marker ] * (i2-i1)
			mat1 = string1[i1:i2] #''.join(patch)
			mat2 = ''.join(patch) #mat1
			Nequal += (i2-i1) # Quality
			string1_pos = min( [string1_pos,j1] )				
		if tag == 'insert':
			mat1 = ''.join( ['+'] * (j2-j1) ) #string2[j1:j2]
			mat2 = string2[j1:j2] #''#.join(patch)
			#if ( (i1 != len(a)) and (i2 != 0) ): # We do not penalize the last insertion of the longest string
# If i1=i2 then nothing happens to that string
			if ( (i1 != len1) and (i2 != 0) ): # We do not penalize the last insertion of the longest string
				delta += (j2-j1)
				Ninsert += (j2-j1)
		if tag == 'delete':
			mat1 = string1[i1:i2]
			mat2 = ''.join( ['-'] * (i2-i1) )
			if not is_linker:
				delta += ( i2-i1 )
				Ndelete += (i2-i1)
			else:
				if j2 != len2:
					delta += ( i2-i1 )
					Ndelete += (i2-i1)
#				if i1 != len1:
#					delta += ( i2-i1 )
#					Ndelete += (i2-i1)
		if tag == 'replace':
			if (j2-j1) >= (i2-i1):
				dst = abs((j2-j1)-(i2-i1))
				mat1 = ''.join( ['>'] * (i2-i1) ) + ''.join( ['_'] * dst ) #string1[j1:j2]
				mat2 = ''.join( ['<'] * (j2-j1) ) #+ ''.join( ['_'] * abs((j2-j1)-(i2-i1)) ) #string2[i1:i2]
			if (j2-j1) < (i2-i1):
				dst = abs((j2-j1)-(i2-i1))
				mat1 = ''.join( ['>'] * (i2-i1) )
				mat2 = ''.join( ['<'] * (j2-j1) ) + ''.join( ['_'] * dst )
			delta += abs( (i2-i1) + (j2-j1) ) / 2
#			Nreplace += abs( (i2-i1) + (j2-j1) )
			Nreplace += max( [(j2-j1), (i2-i1)])

		if Verbose:
			print 'Delta: %i; Matches: %i, Insert: %i, Delete: %i, Replace: %i' %( delta, Nequal, Ninsert, Ndelete, Nreplace )

		virus.append( mat1 )
		str_match.append( mat2 )

#	if trim_end:
#		comparison_len = len(string2)-string1_pos
#	else:
#		if is_linker:
#			min_len = min( [len(string1), len(string2)] )
#			comparison_len = min_len
#		else:
#			comparison_len = len(string1)
			
#	match_quality = 1. - float(delta) / comparison_len
	comparison_len = Ninsert + Ndelete + Nreplace + Nequal
	if Nequal == 0:
		match_quality = 0.
	else:
		#match_quality = 1. - float(delta) / Nequal
		#match_quality = 1. - float(Ninsert+Ndelete+Nreplace) / Nequal
		match_quality = float( Nequal ) / comparison_len
	str_match = ''.join( str_match )
	virus = ''.join( virus )

	if Verbose:
		#print seq_names[i]
		print code+"String1: %s" %string1
		print code+'String2: %s' %string2
		print code+'String:  %s' %str_match
		print code+'Match :  %s' %virus

	trim_pos = -1
	inpos = -1
	finpos = -1
	if match_quality > (1.-tol):
		if not is_linker:
			trim_pos = last_match[1]+last_match[2]
		else:
			trim_pos = last_match[1]

		inpos = last_match[1]
		finpos = last_match[1]+last_match[2]

	if Verbose:
#		print code+'Quality: %f, trimming position: %i \n%sDone.' %(match_quality, trim_pos, code)
		print code+'Quality: %f, trimming position: %i %i \n%sDone.' %(match_quality, inpos, finpos, code)
		
	return match_quality, inpos, finpos
	
# ---------------------------------------------------------------------
def get_match_quality_lev( string1, string2, Verbose=False, marker='v', is_linker=False, tol=0.2 ):
	code = ' >>> GET_MATCH_QUALITY_LEV: '
# Matching first string onto the second one. Quality established based on the minimum
# between the length of the string 1 and string2
	if Verbose:
		print code+'Comparing: \n 1)- %s \n 2)- %s' %(string1, string2)
		print code+'Tolerance = %f' %tol
		#print '                        using difflib (could be sped up by using Levenshtein)'

	len1 = len( string1 )
	len2 = len( string2 )
	# create matching object: operations to be applied to string1 to become string2
	# The last bit is the sequence we are interested in
#	match = difflib.SequenceMatcher( None, string1, string2 )
	# get matching blocks to select last piece to be removed
#	common_blocks = match.get_matching_blocks()
#	if Verbose:
#		print 'Common blocks: ', common_blocks
#	if not is_linker:
#		last_match = common_blocks[-2]
#	else:
#		last_match = common_blocks[0]
#	
#	if Verbose:
#		print 'Last match: ', last_match

	str_match = []
	virus = []
#	opcode = match.get_opcodes()
	opcode = Levenshtein.opcodes( string1, string2 )
	#print opcode
	delta = 0
	Nequal = 0
	Ninsert = 0
	Ndelete = 0
	Nreplace = 0
	comparison_len = 0
	inpos = len2
	finpos = 0
#	string1_pos = len2
	#Reconstructing string and quantifying the mismatch
	for tag, i1, i2, j1, j2 in opcode:
		if Verbose:
			print ("%7s a[%d:%d] (%s) b[%d:%d] (%s)" %(tag, i1, i2, string1[i1:i2], j1, j2, string2[j1:j2]))
		if tag == 'equal':
			# Identifying matching (virus) characters
			patch = [ marker ] * (i2-i1)
			mat1 = string1[i1:i2] #''.join(patch)
			mat2 = ''.join(patch) #mat1
			Nequal += (i2-i1) # Quality
#			string1_pos = min( [string1_pos,j1] )				
			inpos = min( [inpos,j1] )				
			finpos = max( [finpos,j2] )				
		if tag == 'insert':
			mat1 = ''.join( ['+'] * (j2-j1) ) #string2[j1:j2]
			mat2 = string2[j1:j2] #''#.join(patch)
			#if ( (i1 != len(a)) and (i2 != 0) ): # We do not penalize the last insertion of the longest string
# If i1=i2 then nothing happens to that string
			if ( (i1 != len1) and (i2 != 0) ): # We do not penalize the last insertion of the longest string
				delta += (j2-j1)
				Ninsert += (j2-j1)
		if tag == 'delete':
			mat1 = string1[i1:i2]
			mat2 = ''.join( ['-'] * (i2-i1) )
			if not is_linker:
				delta += ( i2-i1 )
				Ndelete += (i2-i1)
			else:
				if j2 != len2:
					delta += ( i2-i1 )
					Ndelete += (i2-i1)
#				if i1 != len1:
#					delta += ( i2-i1 )
#					Ndelete += (i2-i1)
		if tag == 'replace':
			if (j2-j1) >= (i2-i1):
				dst = abs((j2-j1)-(i2-i1))
				mat1 = ''.join( ['>'] * (i2-i1) ) + ''.join( ['_'] * dst ) #string1[j1:j2]
				mat2 = ''.join( ['<'] * (j2-j1) ) #+ ''.join( ['_'] * abs((j2-j1)-(i2-i1)) ) #string2[i1:i2]
			if (j2-j1) < (i2-i1):
				dst = abs((j2-j1)-(i2-i1))
				mat1 = ''.join( ['>'] * (i2-i1) )
				mat2 = ''.join( ['<'] * (j2-j1) ) + ''.join( ['_'] * dst )
			delta += abs( (i2-i1) + (j2-j1) ) / 2
#			Nreplace += abs( (i2-i1) + (j2-j1) )
			Nreplace += max( [(j2-j1), (i2-i1)])

		if Verbose:
			print 'Delta: %i; Matches: %i, Insert: %i, Delete: %i, Replace: %i' %( delta, Nequal, Ninsert, Ndelete, Nreplace )

		virus.append( mat1 )
		str_match.append( mat2 )

#	if trim_end:
#		comparison_len = len(string2)-string1_pos
#	else:
#		if is_linker:
#			min_len = min( [len(string1), len(string2)] )
#			comparison_len = min_len
#		else:
#			comparison_len = len(string1)
			
#	match_quality = 1. - float(delta) / comparison_len
	comparison_len = Ninsert + Ndelete + Nreplace + Nequal
	if Nequal == 0:
		match_quality = 0.
	else:
		#match_quality = 1. - float(delta) / Nequal
		#match_quality = 1. - float(Ninsert+Ndelete+Nreplace) / Nequal
		match_quality = float( Nequal ) / comparison_len
	str_match = ''.join( str_match )
	virus = ''.join( virus )

	if Verbose:
		#print seq_names[i]
		print code+"String1: %s" %string1
		print code+'String2: %s' %string2
		print code+'String:  %s' %str_match
		print code+'Match :  %s' %virus

#	trim_pos = -1
#	inpos = -1
#	finpos = -1
#	if match_quality > (1.-tol):
#		if not is_linker:
#			trim_pos = last_match[1]+last_match[2]
#		else:
#			trim_pos = last_match[1]
#
#		inpos = last_match[1]
#		finpos = last_match[1]+last_match[2]

	if Verbose:
#		print code+'Quality: %f, trimming position: %i \n%sDone.' %(match_quality, trim_pos, code)
		print code+'Quality: %f, trimming position: %i %i \n%sDone.' %(match_quality, inpos, finpos, code)
		
	return match_quality, inpos, finpos
	
# ---------------------------------------------------------------------
class Printer():
    """
    Print things to stdout on one line dynamically
    """
 
    def __init__(self,data):
 
        sys.stdout.write("\r\x1b[K"+data.__str__())
        sys.stdout.flush()
        
# ---------------------------------------------------------------------
def write_sequences( outfile, seq_names, seq_values):

		print 'Saving in %s' %outfile
		nlen = len( seq_names )
		vlen = len( seq_values )
		if nlen != vlen:
			print 'Number of sequences and labels differ. STOP.'
			return
		output = open( outfile,'w' )
		body = []
		for i in range( nlen ):
			body.append( seq_names[i]+'\n' )
			if i<len( seq_names )-1:
				body.append( seq_values[i]+'\n' )
			else:
				body.append( seq_values[i] )

		output.write( ''.join(body) )
		output.close()
		
# ---------------------------------------------------------------------
def read_sequences( infile, Verbose=False ):
	code = ' >>> READ_SEQUENCES: '
	if Verbose:
		print code+'Reading %s' %( infile )
# ------ Reading sequences
	sequence_file = open( infile,'r' )
	sequences = sequence_file.read()
	sequence_file.close()
# ------ Removing new-line character
	if Verbose:
		print code+'Removing empty entries'
	sequences = sequences.strip()
	sequences = sequences.split('\n')
	# Odd lines are names
	seq_names = sequences[0::2]

	# Even lines are values (i.e. sequences)
	seq_values = sequences[1::2]
	# Upper case
	for i in range(len(seq_values)):
		seq_values[i] = seq_values[i].upper()

	# Check length of the names vs sequences to avoid i/o errors
	nseq = len(seq_names)
	nseq_chk = len(seq_values)
	if nseq != nseq_chk:
		print code+'ERROR: Different number of values and names! %6i%6i' %(nseq,nseq_chk)
	else:
		if Verbose:
			print code+'Nseq = %6i' %(nseq)

	return seq_names, seq_values

# ---------------------------------------------------------------------
def merge_files( files_array, outfile='merged.txt', verbose=False):
	code =  '>>> merge_files: '
	body = []
	for file in files_array:
		if verbose:
			print file
		content = open( file, 'r')
		for line in content:	
			body.append( line.strip('\n')+'\n' )
		content.close()
		if verbose:
			print len( body )

	output = open( outfile,'w')
	output.write( ''.join(body[0:-1]) )

	if verbose:
		print code+'Done.'

# ---------------------------------------------------------------------
def get_fasta_format( input_file, verbose=False, out_root='' ):
	code =  '>>> get_fasta_format: '

	sequences = open( input_file,'r' )

	all = []	
	for i,entry in enumerate(sequences):
		if entry[0] != ">":
			entry = entry.strip('\n')
		else:
			if i != 0:
				entry = "\n"+entry
	#print entry
		all.append(entry)

	input = ''.join(all)
	if verbose:
		print input
	if len(out_root) == 0:
		outfile = file.strip('.fna')+'.fasta'
	else:
		outfile = out_root
	if verbose:
		print outfile
	output = open( outfile,'w' )
	output.write( input )
	
	if verbose:
		print code+'Done.'
