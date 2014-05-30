#!/opt/local/bin/ipython

import numpy as np


# ---------------------------------------------------------------------
# def kd1_create_C( pattern, text, k ):
#  
#  	n = len( text )
#  	m = len( pattern )
#  
#  	C = np.zeros( (k-(-1)+1, n-m+2*(k+1)+1), dtype='int')
#  	min_d = 0 #-(k+1)
#  	max_d = n-m+k+1-(-(k+1)) #n-m+k+1
#  	min_e = 0 #-1
#  	max_e = k-(-1) #k
# 
# 	N = 1024
# #	for d in range(n-m+k+1 + 1)+(k+1):
# 	for d in range(k+1,k+1+n-m+k+1 + 1):
# 		print d
# 		C[min_e,d] = d-(k+1)-1
# 
# #	for d in range(k+1+1):
# 	for d in range(-(k+1),1):
# 		print d
# 		C[abs(d)-1,d] = -1
# 		C[abs(d)-2,d] = -N
# 
# 
# 	print C
#	for c in 0 to n-m+k do
#	for e in0 to k do
#	d:=c-e
#	col:=max(C(e-1,d-1)+1, C(e-1,d)+1, C(e-1,d+1))
#	while col<n and col<m+d and x(col+1)=y(col+1-d) do
#	col:=col+1
#	end while
#	C(e,d):=min(col, m+d, n)
#	end for
#	end for

# ---------------------------------------------------------------------
def kd0_create_T( pattern, text, k ):

	n = len( text )
	m = len( pattern )

	T = np.zeros( (m+1,n+1), dtype='int')
	T[0,:] = 0
	T[:,0] = range(m+1)

	for j in range(1,n+1):
		for i in range(1,m+1):
			if pattern[i-1] == text[j-1]:
				p = 0
#				patch = ''.join([patch,pattern[i-1]])
			else:
				p = 1
			diag = T[i-1,j-1]+p
			row = T[i,j-1]+1
			col = T[i-1,j]+1
			paths = [ diag, row, col ]
			#if Verbose:
			#	print paths
			T[i,j] = min( paths )
	
	return T

# ---------------------------------------------------------------------
def kd0_find_end( pattern, text, k, Verbose=False, output_result=False, ending_match=0 ):

	code = ' >>> KD0_FIND_END: '

	m = len( pattern )
	n = len( text )
	if ending_match == 0:
		ending_match = m

	if output_result and Verbose:
		print '\n'+code+'\n'
		print ' - Pattern: %s' %pattern
		print ' - Text:    %s' %text
		print ' - len(P) = %i' %m 
		print ' - len(T) = %i' %n 
		print ' - mismatches = %i' %k 
		print ' - ending_match = %i' %ending_match 

# 	T = np.zeros( (m+1,n+1), dtype='int')
# 	T[0,:] = 0
# 	T[:,0] = range(m+1)
	
# 	if Verbose and False:
# 		print ' - Starting point'
# 		print ''.join(['-']*72)
# 		print T

	matches = []
	epos = -1
	err = m
	lmatch = m
	q = 0
#	patch = ''
	
# 	for j in range(1,n+1):
# 		for i in range(1,m+1):
# 			if pattern[i-1] == text[j-1]:
# 				p = 0
# #				patch = ''.join([patch,pattern[i-1]])
# 			else:
# 				p = 1
# 			diag = T[i-1,j-1]+p
# 			row = T[i,j-1]+1
# 			col = T[i-1,j]+1
# 			paths = [ diag, row, col ]
# 			#if Verbose:
# 			#	print paths
# 			T[i,j] = min( paths )
	T = kd0_create_T( pattern, text, k )

	for j in range(1,n+1):
		if T[m,j] <= k:
			quality = 0
			a = pattern
			b = text[j-m:j]
				#print a
				#print b
			for indx in range(m):
				if a[indx] == b[indx]:
					quality += 1
			matches.append( (j,T[m,j],quality) )
#			matches.append( (j,T[m,j]) )
#			if err >= T[m,j]:
			if (err >= T[m,j]) and (q <= quality):
				err = T[m,j]
				ipos = j-1-(m-T[m,j])
				epos = j-1
				
			if Verbose:
				print ''.join(['-']*72)
				print T
				print ' - Pattern found ending at %i with %i mismatches' %(epos, T[m,j])
				print text
				print text[0:epos+1-(m-T[m,j])]+''.join(['v']*(m-T[m,j]))+text[epos+1:]
				print ''.join(['_']*(epos+1-(m-T[m,j])))+pattern
	
	if Verbose:
		print ''.join(['-']*72)
		t_ar = []
		for i in range(n):
			t_ar.append( text[i] )
		print T

		print matches

	check_partial = ( T[1:,n] == 0 )
	#print check_partial
	#print len(check_partial)
	check_partial[m-1] = False # either one character or the full path: neither have
								# to be accounted for here
	if any( check_partial ):
		for i in range(m):
			if check_partial[i] and ( (i+1) >= ending_match):
				lmatch = (i+1)
				epos = n-1
				err = 0
				partial_match = True
				if Verbose:
					print ' - Possible partial match at the end of the string'
					print '   of length = %i' %lmatch
				break

	if output_result:
		if epos>=0:
			print code
			print ' - output_result'
			print ' - Best match found ending at %i with %i mismatches' %(epos, err)
			print text
			rec = ''.join(['_']*(epos+1) )+text[epos+1:]				
			pat = ''.join(['_']*(epos+1-lmatch) )+pattern
			print rec
			print pat
		else:
			if Verbose:
				print ' - Pattern not found'

	return epos, err
				
# ---------------------------------------------------------------------
def kd0_find_start( pattern, text, k, Verbose=False, output_result=False, ending_match=0 ):

	code = ' >>> KD0_FIND_START: '

	if output_result and Verbose:
		print '\n'+code+'\n'
		print ' - Pattern: %s' %pattern
		print ' - Text:    %s' %text
		print ' - len(P) = %i' %len(pattern) 
		print ' - len(T) = %i' %len(text) 
		print ' - mismatches = %i' %k 
		print ' - ending_match = %i' %ending_match 

	text = text[::-1]

	m = len( pattern )
	n = len( text )

	P = pattern
	mP = len(pattern)
	found = False
	for shift in range(ending_match+1):
		
		pattern = P[0:mP-shift]		
		m = len( pattern )
	
		pattern = pattern[::-1]

		if Verbose:
			print ' - Pattern: %s' %pattern
			print ' - Text: %s' %text
			print ' - len(P) = %i' %m 
			print ' - len(T) = %i' %n 
			print ' - mismatches = %i' %k 
			print ' - ending_match = %i' %ending_match 
			print ' - shift = %i' %shift 

#		T = np.zeros( (m+1,n+1), dtype='int')
#		T[0,:] = 0
#		T[:,0] = range(m+1)
	
# 		if Verbose:
# 			print ' - Starting point'
# 			print ''.join(['-']*72)
# 			print T

		matches = []
		epos = -1
		err = m
		q = 0
#		patch = ''
	
# 		for j in range(1,n+1):
# 			for i in range(1,m+1):
# 				if pattern[i-1] == text[j-1]:
# 					p = 0
# #					patch = ''.join([patch,pattern[i-1]])
# 				else:
# 					p = 1
# 				diag = T[i-1,j-1]+p
# 				row = T[i,j-1]+1
# 				col = T[i-1,j]+1
# 				paths = [ diag, row, col ]
# 				#if Verbose:
# 				#	print paths
# 				T[i,j] = min( paths )
		T = kd0_create_T(pattern, text, k)

		for j in range(1,n+1):
			if T[m,j] <= k:
				quality = 0
				a = pattern
				b = text[j-m:j]
				#print a
				#print b
#				for indx in range(m):
#					if a[indx] == b[indx]:
#						quality += 1
				quality = np.sum( np.array(list(a)) == np.array(list(b)) )
#				print ' quality = ', quality
				matches.append( (j,T[m,j],quality) )
#				if err >= T[m,j]:
# If shift is larger than 0, meaning that we are looking for substrings we want them 
# to be at the end (beginning) of the text. Bug found when running the primer. Linker to be checked
				if shift == 0:
					if (err >= T[m,j]) and (q < quality):
						err = T[m,j]
						q = quality
						ipos = j-1-(m-T[m,j])
						epos = j-1
				else:
					if (err >= T[m,j]) and (q < quality) and (j==quality):
						err = T[m,j]
						q = quality
						ipos = j-1-(m-T[m,j])
						epos = j-1
					
				
				if Verbose:
					print ''.join(['-']*72)
					print T
					print ' - Pattern found ending at %i with %i mismatches' %(epos, err)
					print text
					print text[0:epos+1-(m-T[m,j])]+''.join(['v']*(m-T[m,j]))+text[epos+1:]
					print ''.join(['_']*(epos+1-(m-T[m,j])))+pattern
	
		if Verbose:
			print ''.join(['-']*72)
			t_ar = []
			for i in range(n):
				t_ar.append( text[i] )
			print T

			print matches

		k = max( [0, k-1]) # After full search if not found entirely we reduce k
	# reconstructing left-to-right strings
		if output_result:
			if epos>=0:
				found = True
				if found:
	#			print ''.join(['-']*72)
					pattern = pattern[::-1]
					text = text[::-1]

					epos = n - (epos+1)

					print code
					print ' - output_result for shift = %i' %shift
					print ' - Best match found ending at %i with %i mismatches. Matches %i' %(epos, err, q)
					print text
					rec = text[0:epos]+''.join(['_']*(n-(epos)))
					pat = ''.join(['_']*(epos))+pattern
					print rec
					print pat
					break
			else:
				if Verbose:
					print ' - Pattern not found. Exit.'

	return epos, err
				
# ---------------------------------------------------------------------

if False:
# Web example:
#	P = 'GATAA'
#	#T = 'CGAAGTCGTAACGCGACTAG'
#	T = 'CAGATAAGAGAA'

	P =   'AAGTCGTAAC'
	T = 'CGAAGTCGTAACGCGACTAG'
	T = 'TTAAGTCGTAACTTTTTTTTTTT'
	T = 'TTAAGTiCGTAACTTTTTTTTTTT'
#	T = 'AAGTCGTAACTTTTTTTTTTTTT'
#	T = 'TTTTTTTTTTTAAGTCiGTAACT'
#	T = 'TTTTTTTTTTTTAAGTCGTAACT'
#	T = 'TTTTTTTTTTTTTTTTTTAAGTCGT'
#	T = 'TTTTTTTTTTTTTTTTTTAAGTCGTTTTTT'

#	T = 'AAGTCGTACd'
#	P = 'AGATCGGAAGAGCTCGTATG'
#	T = 'TAGGTCAATTTTAGTATGTTGGTCAAAATCGGAAGAGCTCGTATGCCGTCT'
#	T = 'GGTCAAAATCGGAAGAGCTCGTATGCCGTCT'
	
	[p,e] = kd0_find_end(P,T,1, Verbose=False, output_result=True, ending_match=5 )	
	print p,e
 	print T[p+1:]

 	[p,e] = kd0_find_start(P,T,1, Verbose=False, output_result=True, ending_match=0 )	
 	print p,e
 	print T[0:p]

# 	P = 'adbbc'
# 	T = 'abbdadcbc'
# 	c = kd1_create_C(P,T,2)
# from multiprocessing import Pool
# 
# def f(x):
#     return x*x
# 
# if __name__ == '__main__':
#     pool = Pool(processes=4)              # start 4 worker processes
# #    result = pool.apply_async(f, [10])    # evaluate "f(10)" asynchronously
# #    print result.get(timeout=1)           # prints "100" unless your computer is *very* slow
#     print pool.map(f, range(10))          # prints "[0, 1, 4,..., 81]"



	
