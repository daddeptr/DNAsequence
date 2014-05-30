// BLAST UCSC
// http://genome.ucsc.edu/cgi-bin/hgBlat
// http://www.yuanlei.com/studies/presentations/lit-stringmatch.pdf

#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
//#include <set>
//#include <algorithm>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
	
using namespace std;

string code_description = "-------------------------------------\n\
        Playing with Genome\n\
-------------------------------------\n\n\
 Pattern search using the first algorithm described here:\n\
 http://www.yuanlei.com/studies/presentations/lit-stringmatch.pdf\n\n\
 The code accepts the following input parameters:\n\
  1) -e/-b pattern; specify the string to be searched and whether the \n\
     end or the beginning has to be returned. One of the two must be present;\n\
  2) -i input sequence file;\n\
  3) -o output sequence file. Sequences are trimmed downstream if -e if given\n\
     upstream if -b if provided;\n\
  4) -k number of allowed mismatches;\n\
  5) -r minimum number of bases at the end of the sequence;\n\
  6) -v to display comments.\n";
/*
1a- Look for IR (GTAAACTTCCGACTTCAACTG) (2 mismatches)\n\
1b- Remove everything upstream. These are the good ones\n\
1c- Collect discarded sequences;\n\
2a- Look for linker (AGATCGGAAGAGCTCGTATG) (2 mismatches)\n\
2b- remove linker and everything downstream\n\
3a- Look for primer (GTCCCTTAAGCGGAGCCCT) (2 mismatches)\n\
3b- remove primer and everything downstream\n";
// Working folder
//string dir = "/Users/dpietrob/Desktop/claudia-bio/";	
// IR
string IR = "GTAAACTTCCGACTTCAACTG";
// linker
string linker = "AGATCGGAAGAGCTCGTATG";
// Nested Primer
string primer = "GTCCCTTAAGCGGAGCCCT";
// Tolerance
int tol = 2;
// ------
bool Feedback                           = true;
//add_TA                             = False
bool search_ir                          = false;
bool search_linker                      = false;
bool search_primer                      = true;
bool select_40bases_plus_and_trim_to_40 = false;
bool remove_exact_duplicates            = false;
bool remove_similar                     = false;
bool remove_similar_kd0                 = false;
bool output_files                       = true;
bool check_output                       = false;
*/

// ---------------------------------------------------------------------
int iminval( int iarray[], int n ){
	// Determine length of the array
//	int n;
//	n = sizeof( iarray ) / sizeof(n);
//	cout << "\n - iminval: size of array = " << n;

	int result = iarray[0];

	for (int i=1; i<n; i+=1){
		if (iarray[i]<result){
			result = iarray[i];
		}
	}
	return result;
}

// ---------------------------------------------------------------------
void printMatrix(int M[], int m, int n){
		cout << "\n";
		for (int r=0; r<(m); r+=1){
			for (int c=0; c<(n); c+=1){
				cout << M[r*n+c] << "  ";
			}
			cout << "\n";
		}
}

// ---------------------------------------------------------------------
void kd0_populate_T( string pattern, string text, int k, int T[], bool Verbose=false ){

// Rewriting using 1D array in order to be able to leave the dimension free

	int n = text.length();
	int m = pattern.length();

	if (Verbose){
		cout << "\n - text length = " << n << "\n";
		cout << "\n - pattern length = " << m << "\n";

// Initializing first row of the matrix T equal to 0
		cout << "\n - Initializing first matrix to 0.";
	}
	for (int e=0; e<(n+1)*(m+1); e+=1){
		T[e] = 0;
	}

// Initializing the first column of the matrix to i, skipping (n+1) elements in the vectorized array
	for (int r=0; r<(m+1); r+=1){
		T[r*(n+1)] = r;
	}

	int diag;
	int row;
	int col;
	int p;
	int paths[3];
	
	for (int j=1; j<(n+1); j+=1){
		for (int i=1; i<(m+1); i+=1){
			if (pattern[i-1] == text[j-1]){
				p = 0;
				}
			else{
				p = 1;
			}
			diag = T[(i-1)*(n+1)+(j-1)] + p;
			row = T[i*(n+1)+(j-1)] + 1;
			col = T[(i-1)*(n+1)+j] + 1;

			paths[0] = diag;
			paths[1] = row;
			paths[2] = col;
//			cout << "\n paths:";// << 	sizeof(paths) / sizeof(1);
//			for (int e=0; e<3; e+=1){
//				cout << paths[e];
//			}

			T[i*(n+1)+j] = iminval( paths, 3 );
			
//			cout << "\n p - " << pattern[i-1] << " VS " << text[j-1] << "\n"; 
//			cout << p << "\n";
//			cout << diag << "\t" << row << "\t" << col << "\t" << T[i*(n+1)+j] << "\n";

//			cout << "\n p - " << pattern[i-1] << " VS " << text[j-1] << "\t" << T[i*(n+1)+j];
		}
/*		cout << "\n";
		for (int r=0; r<(m+1); r+=1){
			for (int c=0; c<(n+1); c+=1){
				cout << T[r*(n+1)+c];
			}
			cout << "\n";
		}
*/
//		printMatrix(T,m+1,n+1);
	}

//		printMatrix(T,m+1,n+1);

//	return T;
}

// ---------------------------------------------------------------------
void kd0_find_end( string pattern, string text, int k, int result[2], bool Verbose=false, bool output_result=false, int ending_match=0 ){

	string code = " >>> KD0_FIND_END: ";

	int m = pattern.length();
	int n = text.length();
	if (ending_match == 0){
		ending_match = m;
	}

	if (output_result && Verbose){
		cout << "\n" << code << "\n";
		cout << "\n" << " - Pattern: " << pattern;
		cout << "\n" << " - Text:    " << text;
		cout << "\n" << " - len(P) = " << m;
		cout << "\n" << " - len(T) = " << n; 
		cout << "\n" << " - mismatches = " << k; 
		cout << "\n" << " - ending_match = " << ending_match; 
	}

	int T[(m+1)*(n+1)]; 

//	int matches[];
	int epos = -1;
	int err = m;
	int lmatch = m;
	int q = 0;

	kd0_populate_T( pattern, text, k, T );

	if (Verbose){
		cout << "\n T matrix:\n";
		printMatrix(T,m+1,n+1);
		cout << "\nStarting search:\n";
		cout << "\n Last M row: ";
		for (int j=1; j<(n+1); j+=1){
			cout << T[m*(n+1)+j] << " ";
		}
		cout << "\n";
	}

	for (int j=1; j<(n+1); j+=1){
		if (T[m*(n+1)+j] <= k){
			//cout << "\n" << j << " " << T[m*(n+1)+j] << " " << text[j] << "\n";
			int quality = 0;
			string a = pattern;
			int istr = max(j-m,0);
			string b = text.substr(istr,m);
			//cout << "\n - a = " << a;
			//cout << "\n - b = " << b;
			for (int indx=0; indx<m; indx+=1){
				if (a[indx] == b[indx]){
					quality += 1;
				}
			}
//			matches.append( (j,T[m,j],quality) )
//			matches.append( (j,T[m,j]) )
//			if err >= T[m,j]:
			if ( (err >= T[m*(n+1)+j]) && (q <= quality) ){
				err = T[m*(n+1)+j];
//				ipos = j-1-(m-T[m*(n+1)+j]);
				epos = j-1;
				if (Verbose){
					cout << "\n------0------0------\n";
					//printMatrix(T,m+1,n+1);
					cout << " - Pattern found ending at " << epos << " with " << err <<" mismatches.\n";
					cout << text << "\n";
					cout << "\n";
					for (int indx=0; indx<epos; indx+=1){
						cout << text[indx];
					}
					for (int indx=epos; indx<n; indx+=1){
						cout << text[indx];
					}
					cout << "\n";
//				print ''.join(['_']*(epos+1-(m-T[m,j])))+pattern
				}	
			}
		}
	}		
/*	if Verbose:
		print ''.join(['-']*72)
		t_ar = []
		for i in range(n):
			t_ar.append( text[i] )
		print T

		print matches
*/
	bool check_partial[m];
	for (int i=1;i<=m; i+=1){
//		check_partial[i-1] = false;
//		cout << T[i*(n+1)];
//		if (T[i*(n+1)-1] == 0){
//			check_partial[i-1] = true;
//		}
		check_partial[i-1] = (T[i*(n+1)+n] == 0);
//		cout << check_partial[i-1];
	 }
	//cout << "\n" << sizeof(check_partial)/sizeof(check_partial[0]);
	
	check_partial[m-1] = false; // either one character or the full path: neither have
							   // to be accounted for here
	bool check_sum = false;
	for (int i=0; i<m; i+=1){
		check_sum = ( check_sum || check_partial[i]);
	}
//	cout << "\n";
//	cout << check_sum;
//	cout << "\n";
	if ( check_sum ){
		for (int i=0; i<m; i+=1){
			if (check_partial[i] && ( (i+1) >= ending_match)){
				lmatch = (i+1);
				epos = n-1;
				err = 0;
				//bool partial_match = true;
				if (Verbose){
					cout << "\n - Possible partial match at the end of the string.";
					cout << "\n   of length = " << lmatch;
				}
				break;
			}
		}	
	}
	if (output_result){
		if (epos>=0){
			cout << "\n" << code << "\n";
			cout << "\n - output_result:";
			cout << "\n - Best match found ending at " << epos << " with " << err << " mismatches.";
			cout << "\n" << text;
			string rec = "";
			string pat = "";
			for (int i=0; i<n; i+=1){
				if (i<=epos){
//					rec.append("_");
					rec += "_";
					if (i<=(epos-lmatch)){
//						pat.append("_");
						pat += "_";
					}
					else{
						if ((i-(epos-lmatch)) <= m){
							//cout << i-(epos-lmatch)-1;
//							pat.append(pattern.substr(i-(epos-lmatch)-1,1));
							pat += pattern.substr(i-(epos-lmatch)-1,1);
						}
						else{
//							pat.append(" ");
							pat += " ";
						}
					}
				}
				else{
//					rec.append(text.substr(i,1));
					rec += text.substr(i,1);
				}
			}
			//cout << "\n Out of the loop.";
			cout << "\n" << rec;
			cout << "\n" << pat;
			cout << "\n";
		}
		else{
			if (Verbose){
				cout << "\n - Pattern not found.\n";
			}
		}
	}
	
	result[0] = epos;
	result[1] = err;

	if (Verbose){
		cout << "\n" << code << "End.\n";
	}
}				

// ---------------------------------------------------------------------
void reverseStr( string& x, bool Verbose = false){
	if (Verbose){
		cout << "\n reverseSTR: Beginning...\n";
		cout << x;
	}
	string temp = "";
	int len = x.length();
	if (Verbose){
		cout << "\n string length = " << len;
	}
	for ( int c=(len-1); c>=0; c-=1){
		temp += x[c];
	}
	x = temp;
	if (Verbose){
		cout << "\n" << x << "\n";
		cout << " reverseSTR: Done.";
	}
}

// ---------------------------------------------------------------------
void kd0_find_start( string pattern, string text, int k, int result[2], bool Verbose=false, bool output_result=false, int ending_match=0 ){

	string code = " >>> KD0_FIND_START: ";

	int m = pattern.length();
	int n = text.length();

	if (output_result && Verbose){
		cout << "\n" << code << "\n";
		cout << "\n" << " - Pattern: " << pattern;
		cout << "\n" << " - Text:    " << text;
		cout << "\n" << " - len(P) = " << m;
		cout << "\n" << " - len(T) = " << n; 
		cout << "\n" << " - mismatches = " << k; 
		cout << "\n" << " - ending_match = " << ending_match; 
	}

//	cout << "\n - reverse test: " << text << "\n";
	reverseStr( text );
//	cout << text;
//	cout << "\n";

	string P = pattern;
	int mP = m;
	bool found = false;

	int epos = -1;
	int err = m;
	int lmatch = m;
	int q = 0;

	for (int shift=0; shift<=ending_match; shift+=1){
		
		k = floor( k * float(m-shift)/m );

		string pattern = P.substr(0,mP-shift);
		if (Verbose){
			cout << pattern;
		}		
		m = pattern.length();
	
		reverseStr( pattern );

		if (Verbose){
			cout << "\n" << code << "\n";
			cout << "\n" << " - Pattern: " << pattern;
			cout << "\n" << " - Text:    " << text;
			cout << "\n" << " - len(P) = " << m;
			cout << "\n" << " - len(T) = " << n; 
			cout << "\n" << " - mismatches = " << k; 
			cout << "\n" << " - ending_match = " << ending_match; 
			cout << "\n" << " - shift = " << shift; 
		}
//		matches = []
		int T[(m+1)*(n+1)]; 
	
		epos = -1;
		err = m;
		lmatch = m;
		q = 0;
	
		kd0_populate_T( pattern, text, k, T );
	
		if (Verbose){
			cout << "\n T matrix:\n";
			printMatrix(T,m+1,n+1);
			cout << "\nStarting search:\n";
			cout << "\n Last M row: ";
			for (int j=1; j<(n+1); j+=1){
				cout << T[m*(n+1)+j] << " ";
			}
			cout << "\n";
		}

		for (int j=1; j<=n; j+=1){
			if (T[m*(n+1)+j] <= k){
				int quality = 0;
				string a = pattern;
				int istr = max(j-m,0);
				string b = text.substr(istr,m);
				//cout << "\n - a = " << a;
				//cout << "\n - b = " << b;
				for (int indx=0; indx<m; indx+=1){
					if (a[indx] == b[indx]){
						quality += 1;
					}
				}
//				quality = np.sum( np.array(list(a)) == np.array(list(b)) )
//#				print ' quality = ', quality
//				matches.append( (j,T[m,j],quality) )
//#				if err >= T[m,j]:
// If shift is larger than 0, meaning that we are looking for substrings we want them 
// to be at the end (beginning) of the text. Bug found when running the primer. Linker to be checked
				if (shift == 0){
					if ( (err >= T[m*(n+1)+j]) && (q <= quality) ){
						err = T[m*(n+1)+j];
						epos = j-1;
						}	
				}
				else{
					if ( (err >= T[m*(n+1)+j]) && (q <= quality) && (j==quality)){
						err = T[m*(n+1)+j];
						epos = j-1;
					}
				}	
				
				if (Verbose){
					cout << "\n------0------0------\n";
					//printMatrix(T,m+1,n+1);
					cout << " - Pattern found ending at " << epos << " with " << err <<" mismatches.\n";
					cout << text << "\n";
					cout << "\n";
					for (int indx=0; indx<epos; indx+=1){
						cout << text[indx];
					}
					for (int indx=epos; indx<n; indx+=1){
						cout << text[indx];
					}
					cout << "\n";
//				print ''.join(['_']*(epos+1-(m-T[m,j])))+pattern
				}	
			}
		}
/*		if Verbose:
			print ''.join(['-']*72)
			t_ar = []
			for i in range(n):
				t_ar.append( text[i] )
			print T

			print matches */

//		k = max( 0, k-1); // After full search if not found entirely we reduce k
                          // reconstructing left-to-right strings
/*		if (output_result):
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
					print ' - Pattern not found. Exit.' */
		if (output_result){
			if (epos>=0){
				found = true;
				if (found){
					reverseStr( pattern );
					reverseStr( text );

					epos = n - (epos+1); 

					cout << "\n" << code << "\n";
					cout << "\n - output_result after shift " << shift;
					cout << "\n - Best match found starting at " << epos << " with " << err << " mismatches.";
					cout << "\n" << text;
					string rec = "";
					string pat = "";
					for (int i=0; i<n; i+=1){
						if (i<epos){
							rec += text.substr(i,1);
//							if (i<=(epos-lmatch)){
							pat += ("_");
//							}
//							else{
//								if (i>epos){
//									//cout << i-(epos-lmatch)-1;
//									pat.append(pattern.substr(i-(epos)-1,1));
//								}
//								else{
//									pat.append(" ");
//								}
//							}
						}
						else{
							rec += ("_");
							if ((i-epos)<m){
								pat += pattern.substr(i-epos,1);
							}
						}
					}
					//cout << "\n Out of the loop.";
					cout << "\n" << rec;
					cout << "\n" << pat;
					cout << "\n";
					break;
				}
			}
			else{
				if (Verbose){
					cout << "\n - Pattern not found.\n";
				}
			}
		}

	}
	result[0] = epos;
	result[1] = err;

}				

// --------------------------------------------------------------------
int main(int argc, char **argv){

	cout << code_description;
//	cout << "\nOption selection:\n";
//	cout << " - search_ir =" << search_ir << "\n";

// parsing inputs
       bool eflag = false;
       bool bflag = false;
       bool iflag = false;
       bool kflag = false;
       bool oflag = false;
       bool vflag = false;
       bool rflag = false;
       char *bvalue = NULL;
       char *evalue = NULL;
       char *kvalue = NULL;
       char *ivalue = NULL;
       char *ovalue = NULL;
       char *rvalue = NULL;
       int index;
       int c;
     
       opterr = 0;
     
       while ((c = getopt (argc, argv, "e:b:k:i:o:r:hv")) != -1)
         switch (c){
           case 'b':
             bflag = true;
             bvalue = optarg;
             break;
           case 'e':
             eflag = true;
             evalue = optarg;
             break;
           case 'i':
             iflag = true;
             ivalue = optarg;
             break;
           case 'k':
             kflag = true;
             kvalue = optarg;
             break;
           case 'o':
             oflag = true;
             ovalue = optarg;
             break;
           case 'r':
             rflag = true;
             rvalue = optarg;
             break;
           case 'v':
             vflag = true;
             break;
           case 'h':
             printf( "Four (4) inputs are allowed:\n\
             1) -i input_filename\n\
             2) -o output_filename\n\
             3) -e/b patter_to_searched_at_end/beginning\n\
             4) -k number_of_mismatches\n" );
             break;
           case '?':
             if (optopt == 'b')
               fprintf (stderr, "Option -%c requires an argument.\n", optopt);
             if (optopt == 'e')
               fprintf (stderr, "Option -%c requires an argument.\n", optopt);
             else if (optopt == 'i')
               fprintf (stderr, "Option -%c requires an argument.\n", optopt);
             else if (optopt == 'k')
               fprintf (stderr, "Option -%c requires an argument.\n", optopt);
             else if (optopt == 'o')
               fprintf (stderr, "Option -%c requires an argument.\n", optopt);
             else if (optopt == 'r')
               fprintf (stderr, "Option -%c requires an argument.\n", optopt);
             else if (isprint (optopt))
               fprintf (stderr, "Unknown option `-%c'.\n", optopt);
             else
               fprintf (stderr,
                        "Unknown option character `\\x%x'.\n",
                        optopt);
             return 1;
           default:
             abort ();
           }
     
       printf ("bflag = %d, bvalue = %s\neflag = %d, evalue = %s\niflag = %d,\
        ivalue = %s\noflag = %d, ovalue = %s\nkflag = %d, kvalue = %s\n\
        vflag = %d\n", bflag, bvalue, eflag, evalue, iflag, ivalue, oflag, ovalue, kflag, kvalue, vflag);

		if (bflag && eflag){
			printf(" Warning: both -e and -b option provided. Only one accepted.");
            abort ();
		}     

	for (index = optind; index < argc; index++)
	printf ("Non-option argument %s\n", argv[index]);

	if(!oflag){
		ovalue = "kdiff_strmatch_output.fasta";
	}
	if(!iflag){
		printf(" Warning: no input file provided. Stop.");
		abort();
	}
	int k = intptr_t(kvalue);
	string infile = ivalue;
	string outfile = ovalue;

	int res[2];
	int pos;
	int err;
	int r;
	string text;
	string T; // Added this line  

	if (eflag){

		string P = evalue;
		int m = P.length();
		if(!rflag){
			r = m/2;
		}
		else{
			r = m-intptr_t(rvalue);
		}


		ifstream file( infile.c_str() );

		int cnt = 0;

// Reading sequences in fasta format line by line and performing the search on the fly
		while (!file.eof()){
			getline (file, T);
			if (T.substr(0,1) == ">"){
				text.append (T+"\n"); // Added this line
			}
			else{
				cnt += 1;
				int n = T.length();
				int Table[(m+1)*(n+1)]; 
				kd0_populate_T( P, T, k, Table );
				kd0_find_end( P, T, k, res, false, vflag, r );
				pos = res[0];
				err = res[1];
				if (err <= k && pos>=0){
					text.append (T.substr(pos+1)+"\n"); // Added this piece
				}
				else{
					text.append (T+"\n"); // Added entire line			
				}
			}
			if (!(cnt % 10000)){
				cout << cnt << "\n";
			}
		}

		file.close();

		cout << "\nSaving into "+outfile+"\n";
		ofstream ofile;
		ofile.open ( outfile.c_str() );
		ofile << text;
		ofile.close();
	}

	if(bflag){
		string P = bvalue;
		int m = P.length();
		if(!rflag){
			r = m/2;
		}
		else{
			r = m-intptr_t(rvalue);
		}

		ifstream file( infile.c_str() );

		int cnt = 0;

		while (!file.eof()){
			getline (file, T);
			if (T.substr(0,1) == ">"){
				text.append (T+"\n"); // Added this line
			}
			else{
				cnt += 1;
				int n = T.length();
				int Table[(m+1)*(n+1)]; 
				kd0_populate_T( P, T, k, Table );
				kd0_find_start( P, T, k, res, false, vflag, r ); // Aggressive
				pos = res[0];
				err = res[1];
				if (err <= k && pos>=0){
					text.append (T.substr(0,pos+1)+"\n"); // Added this piece
				}
				else{
					text.append (T+"\n"); // Added entire line			
				}
			}
			if (!(cnt % 10000)){
				cout << cnt << "\n";
			}
		}

		file.close();

		cout << "\nSaving into "+outfile+"\n";
		ofstream ofile;
		ofile.open ( outfile.c_str() );
		ofile << text;
		ofile.close();

	}

	cout << "\nDone. Good bye.";
	return 0;
}
