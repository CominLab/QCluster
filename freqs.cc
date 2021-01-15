#include <math.h>
#include "freqs.hh"
#include <iostream>
#include <cstdlib>
#include <string.h>
#include <unistd.h>
using namespace std;

static bool is_valid_nt(char c)
{
	c = toupper(c);
	if (c == 'A') {
		return true;
	} 
	if (c == 'C') {
		return true;
	} 
	if (c == 'G') {
		return true;
	} 
	if (c == 'T') {
		return true;
	}
	return false;
}


static int nt2int(char c)
{
	switch (c) {
		case 'A':
		case 'a': 
			return 0;
		case 'C':
		case 'c': 
			return 1;
		case 'G':
		case 'g': 
			return 2;
		case 'T':
		case 't': 
			return 3;
	}
	// should never get here
	return -1;
}

//Divide each element of the vector by the norm1 of the vector itself
void normalize_row(double* fv, int L)
{
	double total = 0;
	for (int l=0; l<L; l++){
		total += fv[l];	// ||fv||1
	}
	for (int l=0; l<L; l++){
		fv[l] /= total;
	}
	return;
}

//Divide each element of "qual" vector by the norm1 of "freq" vector
//Freq vector MUST NOT be normalized
void normalize_row_against(double* freq, double* qual, int L)
{
	double total = 0;
	for (int l=0; l<L; l++){
		total += freq[l]; // ||freq||1
	}
	for (int l=0; l<L; l++){
		qual[l] /= total;
	}
	return;
}

//Project Freq vector on the hyperplane x1 + x2 + x3... -1 = 0
//Freq and qual vector MUST be normalized (with normalize_row_against)
void normalize_row_projection(double* freq, double* qual, int L)
{
	double total = 0, total_qual = 0;
	for (int l=0; l<L; l++){
		total += freq[l];
		total_qual += qual[l];
	}
	double parziale = (total - total_qual) / L;
	for (int l=0; l<L; l++){
		qual[l] += parziale;
	}
	return;
}

//Note: freq_1 MUST NOT be normalized
void fill_overlap_count_vector(string seq, string seq_qual, int K,
						double* freq, double *quality_vector, int normalize,
						double pc, double *avg_quality_1, double **avg_quality, 
						double *freq_1, bool redistribute)
{
	int L = seq.length();
	int *kmer = new int[K];
	int *qual = new int[K];
	double *prob = new double[K];
	double *freq_bases = new double[NUM_NT];
	
	//If K==1, freq_1 is NULL and freq_bases is initialized to 1/4 for each base
	if (K==1){
		for (int i=0; i<NUM_NT; i++) freq_bases[i] = 1.0/4.0;
	}//Otherwise, freq_bases is the normalization of freq_1
	else{
		double somma =0;
		for (int i=0; i<NUM_NT; i++) somma += freq_1[i];
		for (int i=0; i<NUM_NT; i++) freq_bases[i] = freq_1[i]/somma;
	}
	
	double readqual=1;
	for(int i=0; i<L; ++i){
		seq[i] = toupper(seq[i]);
	}
	int N = 1;  // number of K-mers
	for(int k=0; k<K; ++k){
		N *= NUM_NT;
	}
	for(int i=0; i<N; ++i){
		freq[i] = pc;  // initialize with the pseudocount
		quality_vector[i] = pc;
	}
	bool valid_kmer = true;
	int index_kmer = 0;
	
	for (int i=0; i<L-K+1; i++){
		valid_kmer = true;
		index_kmer = 0;
		readqual = 1;
		//Fill the kmer, qual and prob vectors
		for (int j=0; j<K; j++){
			if (!is_valid_nt(seq[i+j])){
					valid_kmer = false;
					break;
				}
			kmer[j] = nt2int(seq[i+j]);
			qual[j] = int(seq_qual[i+j])-base;
			prob[j] = (1.0 - pow(10.0,-(qual[j])/10.0));
			readqual *= prob[j];
			//Calculate vector index of the kmer
			index_kmer = NUM_NT*index_kmer + kmer[j];
		}

		if (!valid_kmer) continue;

		freq[index_kmer] += 1;
		quality_vector[index_kmer] += readqual;
		
		if (avg_quality != NULL)
			for (int j=0; j<K; j++) avg_quality[index_kmer][kmer[j]] += qual[j];
			
		if (avg_quality_1 != NULL) 
			avg_quality_1[kmer[0]] += qual[0];
		
		if (!redistribute) continue;
		
		//Redistributing quality: for each letter in the kmer...
		for (int j=0; j<K; j++){
			int original_letter = kmer[j];
			double original_prob = prob[j];
			//...we try to replace it with another
			for (int letter=0; letter<NUM_NT; letter++){
				//we must redistribute also to the same letter
				//if (letter == original_letter) continue;
				kmer[j] = letter;
				//Equally subdivided probability
				//prob[j] = (1.0-original_prob) / double(NUM_NT-1.0);
				//Probability subidivided basing on the letter frequency
				prob[j] = (1.0-original_prob) * double(freq_bases[letter]);
				
				//Recalculation of the new kmer probability
				readqual = 1;
				index_kmer = 0;
				for (int f=0; f<K; f++){
					readqual *= prob[j];
					index_kmer = NUM_NT*index_kmer + kmer[j];
				}
				//Add a little probability to that kmer
				quality_vector[index_kmer] += readqual;
			}
			kmer[j] = original_letter;
			prob[j] = original_prob;
		}
	}
	
	
	switch(normalize){
		case 1:	
			normalize_row(freq, N);
			normalize_row(quality_vector, N);
			break;
		case 2:
			normalize_row_against(freq, quality_vector, N);
			normalize_row(freq, N);
			break;
		case 3:
			normalize_row_against(freq, quality_vector, N);
			normalize_row(freq, N);
			normalize_row_projection(freq, quality_vector, N);
			break;
		default:
			break;
	}
	
	delete[] kmer;
	delete[] qual;
	delete[] prob;
	if (K==1) delete[] freq_bases;
	return;
}




/**Function responsible for the calculation of the quality expected value for 
 * each k-word: the expected value is computed basing on the full dataset.
 * Refer to the documentation for the description of the possible methods. 
 * Default method: E1 */
void calculate_quality_expected_value(int method, int N, int K, int L, 
			double **freq, double **quality, double **freq_1,
			 double* avg_quality_1, double **avg_quality, double *expected_qual){
	switch (method){
	case 2:
		{
		int *times_to_use = new int[NUM_NT];
		for (int i=0; i<L; i++){
			for(int j=0; j<NUM_NT; j++) times_to_use[j] = 0;
			int divisore = 1;
			for(int k=0; k<K; k++){
				times_to_use[(i/divisore)%(NUM_NT)] += 1;
				divisore *= NUM_NT;
			}
			for(int j=0; j<NUM_NT; j++)
				if (times_to_use[j]!=0) avg_quality[i][j] /= times_to_use[j];
		}

		delete[] times_to_use;

		for (int i=0; i<L; i++){
			int divisore = 1;
			expected_qual[i] = 1;
			for(int k=0; k<K; k++){
				double num_occorrenze = 0;
				for (int j=0; j<N; j++) {num_occorrenze += freq[j][i];}
				expected_qual[i] *= 1 - 
					pow(10.0, -(avg_quality[i][(i/divisore)%(NUM_NT)])
									/(10.0*num_occorrenze));
				divisore *= NUM_NT;
			}
		}
		break;
		}
	case 3:
		{
		for(int i=0; i<NUM_NT; i++) {
			int tot = 0;
			for (int j=0; j<N; j++) tot += freq_1[j][i];
			avg_quality_1[i] = 1-pow(10.0,-avg_quality_1[i]/(10.0*tot));
		}
		int divisore = 1;
		
		for (int i=0; i<L; i++) expected_qual[i] = 1;
		for(int k=0; k<K; k++){
			for (int i=0; i<L; i++){
				expected_qual[i] *= avg_quality_1[(i/divisore)%(NUM_NT)];
			}
			divisore *= NUM_NT;
		}
		break;
		}
	case 1: //fall
		
	default:
		{
		for (int i=0; i<L; i++){ 
			expected_qual[i] = 0; 
			double denominatore = 0; 
			for (int j=0; j<N; j++) { 
				expected_qual[i] += quality[j][i]; 
				denominatore += freq[j][i]; 
			} 
			expected_qual[i] /= denominatore; 
		} 
		break;	
		}
	}//END SWITCH
	
}//END FUNCTION


void expected_frequency_p2global(int N, int K, int L, double *expected_freq,
								double **freq_1)
{
	double nt_freq[NUM_NT];
	memset(nt_freq, 0, NUM_NT * sizeof(nt_freq[0]));
	double S=0;
	// compute frequencies of individual nucleotides
	for (int n=0; n<N; n++){
		for (int l=0; l<NUM_NT; l++){
			nt_freq[l] += freq_1[n][l]; //= count of letter l on the dataset
		}
	}
	for (int l=0; l<NUM_NT; l++)
		S += nt_freq[l];	//S = total number of letters in the dataset

	for (int l=0; l<NUM_NT; l++){
		nt_freq[l] /= S;	//Average frequency of each letter in the dataset
	}
	// compute frequencies of words using Markov model
	for(int l=0; l<L; l++){
		double p = 1;
		for(int M=1; M<L; M*=NUM_NT){
			int digit = l/M % NUM_NT;
			p *= nt_freq[digit];
		}
		expected_freq[l] = p;
	}
	return;
}


void expected_frequency_p1global(int N, int K, int L, double *expected_freq,
								double **freq)
{
	int tot_parole = 0;
	for (int l=0; l<L; l++) {
		expected_freq[l] = 0;
		for (int n=0; n<N; n++)
			expected_freq[l] += freq[n][l];
		tot_parole += expected_freq[l];
	}	
	for (int l=0; l<L; l++) {
		expected_freq[l] /= tot_parole;
	}
	return;
}

// normalize frequency matrix to make its columns univariant
void normalize_freq_matrix(double** freq, double** qual, int N, int row_length)
{
	double tmp, sum_freq, sum_freq_square, V_freq;
	double sum_qual, sum_qual_square, V_qual;
	for(int l=0; l<row_length; ++l){
		sum_freq_square=0;
		sum_freq=0;
		sum_qual_square = 0;
		sum_qual = 0;
		for(int n=0; n<N; ++n){
			tmp = *(*(freq+n)+l);
			sum_freq_square += tmp*tmp;
			sum_freq += tmp;
			tmp = *(*(qual+n)+l);
			sum_qual_square += tmp*tmp;
			sum_qual += tmp;
		}
		//Variance:
		V_freq = sqrt( sum_freq_square/N - (sum_freq*sum_freq)/(N*N) );  
		V_qual = sqrt( sum_qual_square/N - (sum_qual*sum_qual)/(N*N) ); 
		if (V_freq>0 && V_qual>0){
			for(int n=0; n<N; ++n){
				*(*(freq+n)+l) /= V_freq;
				*(*(qual+n)+l) /= V_qual;				
			}
		}
	}
	return;
}
