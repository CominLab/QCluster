#include <cstdlib>
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "seqio.hh"
#include "freqs.hh"
#include "em.hh"
using namespace std;


enum {PC=1};  // default for pseudocount
enum {TRIALS=100};  // number of bootstrap datasets for consensus clustering





void print_usage(char* progname)
{
	cout << "Centroid based (k-means-like) clustering of sequences "
		"in n-mer frequency space\n\n"
	"usage: "<<progname<<" [-h] [-c num_clusters] [-d dist_type] "
	"[-e num_method] [-k nmer_length] [-m max_iterations] [-N num_method] "
	"[-n] [-P num_method] [-p pseudocount] [-R] [-r] [-S seed] [-t num_trials]"
	" [-v] [-w] fastq_file\n\n"
	"\t-c num_clusters (5 clusters by default)\n"
	"\t-d dist_type:\n"
		"\t\ta: d2* distance\n"
		"\t\tc: chi square statistic\n"
		"\t\td: d2 distance\n"
		"\t\te: regular euclidean (L2) distance; default\n"
		"\t\tk: Kullback-Leibler divergence\n"
		"\t\ts: symmetrized Kullback-Leibler divergence\n"
	"\t-e num_method: method for computation of quality expected \n"
	"\t   value of words:\n"
		"\t\t1: average quality of the word over the full dataset (default)\n"
		"\t\t2: average quality of each base in the same word over the full \n"
		"\t\t   dataset and calculate the expected value with Markovian model\n"
		"\t\t3: average quality of each baseover the full dataset and \n"
		"\t\t   calculate the expected value with Markovian model\n"
	"\t-h print this message and exit\n"
	"\t-k nmer_length: length of word (2-mers by default)\n"
	"\t-m max_iterations: number of maximum iterations of the algorithm \n"
	"\t   without an improvement (min. 2) \n"
	"\t-N num_method: Divide each quality vector by the taxicab norm so \n"
	"\t   that the sum of the elements is (about) 1. Possible values:\n"
		"\t\t0: no normalization\n"
		"\t\t1: divide by the norm of the vector itself\n"
		"\t\t2: divide by the norm of the frequency vector (default)\n"
		"\t\t3: as method 2 and project the vector on the hyperplane 1*X-1=0\n"
	"\t-n normalize frequency matrix to make each column univariant; \n"
	"\t   implies L2 distance\n"
	"\t-P num_method: method for computation of expected frequancy of words:\n"
		"\t\t1: uses a markovian model to compute the expected frequency of\n"
		"\t\t   words basing only on the words of the cluster\n"
		"\t\t   (P2Local methond: default)\n"
		"\t\t2: Average frequency of every word over the entire dataset\n"
		"\t\t   (P1Global method)\n"
		"\t\t3: Markovian model based on average frequency of single bases\n"
		"\t\t   on the full dataset (P2Global method)\n"
	"\t-p pseudocount: (default: "<<PC<<"); can be fractional. \n"
	"\t   Must be nonzero with KL and simmetrized KL distance.\n"
	"\t-R do not redistribute missing quality among other bases. \n"
	"\t-r reverse complement and stack together\n"
	"\t-S seed: initial seed for random number generator\n"
	"\t-t num_trials: repeat clustering num_trials, choosing the best \n"
	"\t   partitioning. Default: 1\n"
	"\t-v output progress messages (repeat for increased verbosity)\n"
	"\t-w write sequences from each cluster to a file"<<endl;
	return;
}


int main(int argc, char **argv)
{
	int num_clusters = 5;  // target number of clusters
	char dist_type = 'e'; // default to L2
	int K = 2;  // k-mer size
	bool normalize_matrix_flag = false; // make columns univariant ("whiten");
		// off by default
	int normalize = 2; // normalize each row. Deafult: method 2
	bool redistribute = true; //Redistribute missing quality by default
	int e_method = 1;
	int p_method = 1;
	double pseudocount = PC; // pseudocount for k-mers
	bool rc_flag = false; // append reverse complement to each sequence;
		// off by default
	unsigned int random_seed = time(NULL); // initial seed
	int num_trials = 0; // later we check if it was specified using "-t" arg
		// and set the correct default if it was not; note that the default 
		// depends on whether consensus clustering is performed
	int verbose_level = 0;
	bool write_flag = false; // write sequences by cluster; off by default
	int default_max_iterations = 1000;

	int opt;
	while((opt = getopt(argc, argv, "c:d:e:k:m:N:nP:p:q:rRsS:t:vwh")) != -1 ){
		switch(opt){
            case 'c':
                num_clusters = atoi(optarg);
                break;
			case 'd':
                dist_type = *optarg;
                break;
			case 'e':
				e_method = atoi(optarg);
				break;
            case 'k':
                K = atoi(optarg);
                break;
			case 'm':
				max_iterations = atoi(optarg);
				//If max_iterations < 2 the default value will be restored
				if (max_iterations<2) max_iterations = default_max_iterations;
				break;
			case 'N':
				normalize = atoi(optarg);
				break;
            case 'n':
                normalize_matrix_flag = true;
                break;
			case 'P':
				p_method = atoi(optarg);
				if (p_method<1 || p_method>3) p_method = 2;
				break;				
			case 'p':
				pseudocount = atof(optarg);
				break;
			case 'R':
				//Do not redistribute quality
				redistribute = false;
				break;
            case 'r':
                rc_flag = true;
                break;
			case 'S':
				random_seed = atoi(optarg);
				break;
            case 't':
                num_trials = atoi(optarg);
                break;
            case 'v':
                ++verbose_level;
                break;
            case 'w':
                write_flag = true;
                break;
            case 'h':
				print_usage(argv[0]);
				return EXIT_SUCCESS;
            case '?':
                print_usage(argv[0]);
                return EXIT_FAILURE;
            default:
                // You won't actually get here
                break;
        }
	}

	if (dist_type != 'a'){
		p_method = 0;//No calculation of expected freq and expected quality
		e_method = 0;
	}

	// if "-t" argument was absent, use the appropriate default
	if (!num_trials) {
		num_trials = 1;
	}

	if(argc != optind + 1){
		cout<<"Missing/extra input file name"<<endl;
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}
	char* fastq_file_name = *(argv + optind );

	// initialize random number generator
	srand(random_seed);

	// initialize constants
	int L = 1; // length of the freq. vector L = NUM_NT**K
	for(int k=0; k<K; ++k){
		L *= NUM_NT;
	}

	// count sequences
	int N = CountReads(fastq_file_name);  // number of sequences
	if (verbose_level>0){
		cerr<<"Found "<<N<<" sequences\n";
	}

	// read sequences and build count matrix
	RecordGenerator rec_gen(fastq_file_name);
	double **freq = new double*[N];
	freq[0] = new double[N * L];
	double **quality = new double*[N];
	for (int i=0; i<N; i++) quality[i] = new double[L];
	
	double **freq_1  = new double*[N];
	freq_1[0] = new double[N * NUM_NT];
	double **quality_1 = new double*[N];
	for (int i=0; i<N; i++) quality_1[i] = new double[NUM_NT];
	
	//Vector for calculus of E(Pxi)2
	double *avg_quality_1 = new double[NUM_NT];
	for(int i=0; i<NUM_NT; i++) avg_quality_1[i] = 0;
	
	//Vector for calculus of E(Pxi)3
	double **avg_quality = new double*[L];
	for (int i=0; i<L; i++) {
		avg_quality[i] = new double[NUM_NT];
		for (int j=0; j<NUM_NT; j++) avg_quality[i][j] = 0;
	}
	
	//For each read we count the number of kmers
	for(int i=0; i<N; ++i){
		SeqRecord rec = rec_gen.next();
		string seq = rc_flag ? rec.seq() + rec.rc().seq() : rec.seq();
		freq[i] = freq[0] + i*L;

		//Calculate frequencies of individual bases
		freq_1[i] = freq_1[0] + i*NUM_NT;
		fill_overlap_count_vector(seq, rec.qual(), 1, freq_1[i], 
								  quality_1[i], false, pseudocount, 
								  avg_quality_1, NULL, NULL, false);
		
		//Calculate frequencies of kmers
		fill_overlap_count_vector(seq, rec.qual(), K, freq[i], quality[i],
								normalize, pseudocount, NULL, avg_quality, 
								freq_1[i], redistribute);
	}
	
	if (verbose_level>0){
		cerr<<"Read counts calculated\n";
	}

	
	//Compute the expected quality of each kmer
	double *expected_qual = NULL;
	if (e_method != 0) {
		expected_qual = new double[L];
		calculate_quality_expected_value(e_method, N, K, L, freq, quality,
						freq_1, avg_quality_1, avg_quality, expected_qual);
	}
	
	//Expected frequancy of each kmer
	double *expected_freq = NULL;
	//Instantiate the vector only if p_method is global (P1G or P2G)
	if (p_method==2 || p_method==3) expected_freq = new double[L];
	
	//P1G Average frequency of every word over the entire dataset
	if (p_method==2)
		expected_frequency_p1global(N, K, L, expected_freq, freq);
	
	//P2G Markovian model based on average frequency of single bases
	if (p_method==3)
		expected_frequency_p2global(N, K, L, expected_freq, freq_1);
	

	
	
	// whiten if requested; in this case reset the distance to euclidean L2
	if (normalize_matrix_flag){
		dist_type = 'e';
		normalize_freq_matrix(freq, quality, N, L);
	}
	if (verbose_level > 0){
		cerr<<"Word frequencies calculated\n";
	}


	// hard EM clustering
	int* assignment = new int[N];
	double* Z = new double[num_clusters * N];
		if (verbose_level > 0) {
			cerr<<"Running regular EM clustering "<<num_trials<<" times "
			"to chose the best partitioning"<<endl;
		}
		num_clusters = hard_em(num_clusters, N, L, freq, quality, expected_qual, 
			quality_1, expected_freq, NUM_NT, freq_1, assignment, Z, num_trials, 
			dist_type, verbose_level);
		
	if (verbose_level > 0){
		cerr<<"Clustering done\n";
	}

	// release memory previously allocated for frequency data
	delete[] freq[0];
	delete[] freq;
	delete[] freq_1[0];
	delete[] freq_1;
	
	for(int i=0; i<N; i++) {
		delete[] quality[i];
		delete[] quality_1[i];
	}
	for(int i=0; i<L; i++)
		delete[] avg_quality[i];
	delete[] quality;
	delete[] quality_1;
	delete[] avg_quality_1;
	delete[] avg_quality;
	if (p_method==2 || p_method==3) delete [] expected_freq;
	if (e_method != 0) delete[] expected_qual;

	// reassign cluster names so that they are sorted 
	// in the descending order of the number of members
	int* num_members = new int[num_clusters];
	int* ix = new int[num_clusters];
	int* rix = new int[num_clusters];
	for (int i=0; i<num_clusters; ++i){
		*(num_members + i) = 0;
		*(ix + i) = i;
	}
	for (int i=0; i<N; ++i){
		num_members[assignment[i]]++;
	}
	for (int i=0; i<num_clusters; ++i){
		int max_pos = i;
		for (int j=i+1; j<num_clusters; ++j){
			if (num_members[ix[j]] > num_members[ix[max_pos]]) {
				max_pos = j;
			}
		}
		int tmp = ix[i];
		ix[i] = ix[max_pos];
		ix[max_pos] = tmp;
	}
	for (int i=0; i<num_clusters; ++i){
		rix[ix[i]] = i;
	}
	for (int i=0; i<N; ++i) {
		*(assignment + i) = *(rix + *(assignment + i));
	}
	delete[] num_members;
	//delete[] ix and rix later

	cout.setf(ios_base::fixed);
	cout.precision(4);

	// output assignment to stdout
	RecordGenerator rec_gen_2(fastq_file_name);
	for(int n=0; n<N; ++n){
		cout<<rec_gen_2.next().id()<<"\t"<<*(assignment+n);
		cout<<endl;
	}
	delete[] ix;
	delete[] rix;

	// write sequences to individual files if requested
	if (write_flag) {
		vector<SeqRecord> rec_vec;
		FastqRead(fastq_file_name, rec_vec);
		int N = rec_vec.size();
		char fname[256];
		for (int i=0; i<num_clusters; ++i){
			sprintf(fname, "part_%i.fastq", i+1);
			ofstream ofs(fname);
			for (int j=0; j<N; ++j) {
				if (*(assignment + j) == i) {
					ofs<<"@"<<rec_vec[j].desc()<<endl;
					ofs<<rec_vec[j].seq()<<endl;
					ofs<<"+"<<endl;
					ofs<<rec_vec[j].qual()<<endl;
				}
			}
		}
	}

	// free allocated memory
	delete[] Z;
	delete[] assignment;
	return EXIT_SUCCESS;
}

