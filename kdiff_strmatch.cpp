// BLAST UCSC
// http://genome.ucsc.edu/cgi-bin/hgBlat
// http://www.yuanlei.com/studies/presentations/lit-stringmatch.pdf

#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
//#include <set>
#include <algorithm>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

//#include <boost/algorithm/string.hpp>
//#include <vector>

#include "kdiff_lib.hpp"
	
using namespace std;
//using namespace boost;

string code_description = "-------------------------------------\n\
        Playing with Genome\n\
-------------------------------------\n\n\
 Pattern search using the first algorithm described here:\n\
 http://www.yuanlei.com/studies/presentations/lit-stringmatch.pdf\n\n\
 The code accepts the following input parameters:\n\
  1) -e/-b pattern; specify the string to be searched and whether the \n\
     downstream or the upstream has to be returned. Either one can be present;\n\
  2) -i input sequence file;\n\
  3) -o output sequence file. Sequences are trimmed downstream if -e if given\n\
     upstream if -b if provided;\n\
  4) -k number of allowed mismatches;\n\
  5) -r minimum number of bases at the end of the sequence;\n\
  6) -t length of the output sequence (trimming);\n\
  7) -v to display comments.\n";
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

// --------------------------------------------------------------------
int main(int argc, char** argv){

	cout << code_description;
//	vector <string> fields;
//	string cmds = argv.substr(0, argv.find()' ');
//	split( fields, s, is_any_of( " " ) );
//	cout << "\nArguments: " << argv[1];

//	int argc = 10;

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
       bool tflag = false;
       
       char *bvalue = NULL;
       char *evalue = NULL;
       char *kvalue = NULL;
       char *ivalue = NULL;
       char *ovalue = NULL;
       char *rvalue = NULL;
       char *tvalue = NULL;
/*
       string bvalue;// = NULL;
       string evalue;// = NULL;
       string kvalue;// = NULL;
       string ivalue;// = NULL;
       string ovalue;// = NULL;
       string rvalue;// = NULL;
*/
//       int index;
       char *c = NULL;
     
       opterr = 0;
     
//       while ((c = getopt (argc, argv, "e:b:k:i:o:r:hv")) != -1)
		if(argc > 1){
			printf("\n - reading parameters...\n");
			for (int i=1; i<argc; ++i){
				c = argv[i];
				if( c[0] == '-' ){
					printf(" Valid option: %s\n", c);
					printf(" command value: %s\n", argv[i+1]);
					if(c[1] == 'b'){
						 bflag = true;
						 bvalue = argv[i+1];
					}
					else if(c[1] == 'e'){
						 eflag = true;
						 evalue =  argv[i+1];
					}
					else if(c[1] == 'i'){
						 iflag = true;
						 ivalue =  argv[i+1];
					}
					else if(c[1] == 'k'){
						 kflag = true;
						 kvalue =  argv[i+1];
					}
					else if(c[1] == 'o'){
						 oflag = true;
						 ovalue =  argv[i+1];
					}
					else if(c[1] == 'r'){
						 rflag = true;
						 rvalue =  argv[i+1];
					}
					else if(c[1] == 't'){
						 tflag = true;
						 tvalue =  argv[i+1];
					}
					else if(c[1] == 'v'){
						 vflag = true;
					}
					else if(c[1] == 'h'){
//						 printf( "Four (4) inputs are allowed:\n\
//						 1) -i input_filename\n\
//						 2) -o output_filename\n\
//						 3) -e/b patter_to_searched_at_end/beginning\n\
//						 4) -k number_of_mismatches\n" );
						 cout << code_description;
					}
					else {
//						 printf( "At least (4) inputs are allowed:\n\
//						 1) -i input_filename\n\
//						 2) -o output_filename\n\
//						 3) -e/b patter_to_searched_at_end/beginning\n\
//						 4) -k number_of_mismatches\n" );
						 cout << code_description;
						 return 1;
						 abort ();
					}
/*
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
*/
	 			}
			}
			if (bflag && eflag){
				printf(" Warning: both -e and -b option provided. Only one accepted. Stop.");
				abort ();
			}
		}
		else{
			printf("No input provided. Terminated.");
			return 1;
		}

//		printf ("bflag = %d, bvalue = %s\neflag = %d, evalue = %s\niflag = %d,\
		ivalue = %s\noflag = %d, ovalue = %s\nkflag = %d, kvalue = %s\n\
		vflag = %d\n", bflag, bvalue, eflag, evalue, iflag, ivalue, oflag, ovalue, kflag, kvalue, vflag);

//	for (index = optind; index < argc; index++)
//	printf ("Non-option argument %s\n", argv[index]);

	if(!oflag){
		ovalue = "kdiff_strmatch_output.fasta";
	}
	if(!iflag){
		printf(" Warning: no input file provided. Stop.");
		abort();
	}

	int k = 0;
	if(eflag || bflag){
		int k = atoi(kvalue);
		printf(" - number of mismatches allowed: %d\n", k);
	}

	int len = -1;
	if(tflag){
		len = atoi(tvalue);
		printf(" - length of the returned sequence: %d\n", len);
	}
//	else{
//		int len = -1;
//	}
	
	string infile = ivalue;
	string outfile = ovalue;

	int res[2];
	int pos;
	int err;
	int r = 0;
	string text;
	string T; // Added this line  

//	printf(" - ending match = %d\n", r);
	if (rflag) {r = atoi(rvalue); }
//	printf(" - ending match = %d\n", r);

	if(eflag){

		string P = evalue;
		int m = P.length();
		if(!rflag){
			r = m;
		}
		else{
			r = m - r;
		}

		printf(" - ending match = %d\n", r);
		ifstream file( infile.c_str() );

		int cnt = 0;

//		string tmp = "";
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
					if(!tflag){
						text.append (T.substr(pos+1)+"\n"); // Added this piece
					}
					else{
						string tmp = T.substr(pos+1);
						int tmpl = tmp.length();
						tmpl = min(len,tmpl);
						text.append (tmp.substr(0,tmpl)+"\n"); // Added this piece					
					}
				}
				else{
					if(!tflag){
						text.append (T+"\n"); // Added entire line
					}
					else{
						text.append (T.substr(0,min(len,n))+"\n"); // Added this piece										
					}
//				delete[] Table;
				}
			}
			if (!(cnt % 10000) && (cnt>0)){
				cout << cnt << "\n";
			}
		}

		file.close();

		cout << "\nSaving into "+outfile+"\n";
		ofstream ofile;
		ofile.open ( outfile.c_str() );
		ofile << text;
		ofile.close();

		printf("Processed %i sequences.", cnt);

	}
	else if(bflag){
		string P = bvalue;
		int m = P.length();
		if(!rflag){
			r = 0;
		}
		else{
			r = m - r;
		}

		printf(" - ending match = %d\n", r);
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
					if(!tflag){
						text.append (T.substr(0,pos)+"\n"); // Added this piece
					}
					else
					{
						text.append (T.substr(0,min(len,n))+"\n"); // Added this piece					
					}
//				delete[] Table;
				}
				else{
					if(!tflag){
						text.append (T+"\n"); // Added entire line
					}
					else
					{
						text.append (T.substr(0,min(len,n))+"\n"); // Added this piece										
					}
				}
			}
			if (!(cnt % 10000) && (cnt>0)){
				cout << cnt << "\n";
			}
		}

		file.close();

		cout << "\nSaving into "+outfile+"\n";
		ofstream ofile;
		ofile.open ( outfile.c_str() );
		ofile << text;
		ofile.close();
		
		printf("Processed %i sequences.", cnt);

	}
	else if(tflag){

//		string P = evalue;
//		int m = P.length();
//		if(!rflag){
//			r = m;
//		}
//		else{
//			r = m - r;
//		}
//
//		printf(" - ending match = %d\n", r);
		ifstream file( infile.c_str() );

		int cnt = 0;

//		string tmp = "";
// Reading sequences in fasta format line by line and performing the search on the fly
		while (!file.eof()){
			getline (file, T);
			if (T.substr(0,1) == ">"){
				text.append (T+"\n"); // Added this line
			}
			else{
				cnt += 1;
				int n = T.length();
//				int Table[(m+1)*(n+1)]; 
//				kd0_populate_T( P, T, k, Table );
//				kd0_find_end( P, T, k, res, false, vflag, r );
//				pos = res[0];
//				err = res[1];
//				if (err <= k && pos>=0){
//					if(!tflag){
//						text.append (T.substr(pos+1)+"\n"); // Added this piece
//					}
//					else{
//						string tmp = T.substr(pos+1);
//				int tmpl = tmp.length();
				int tmpl = min(len,n);
				text.append (T.substr(0,tmpl)+"\n"); // Added this piece					
//					}
//				}
//				else{
//					if(!tflag){
//						text.append (T+"\n"); // Added entire line
//					}
//					else{
//						text.append (T.substr(0,min(len,n))+"\n"); // Added this piece										
//					}
//				delete[] Table;
//				}
			}
			if (!(cnt % 10000) && (cnt>0)){
				cout << cnt << "\n";
			}
		}

		file.close();

		cout << "\nSaving into "+outfile+"\n";
		ofstream ofile;
		ofile.open ( outfile.c_str() );
		ofile << text;
		ofile.close();

		printf("Processed %i sequences.", cnt);

	}


	cout << "\nDone. Good bye.";
	return 0;
}
