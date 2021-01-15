#include <string.h> // memset
#include <math.h>
#include<iostream>
using namespace std;


// mean centroid
void mean_centroid(int N, int row_length, double** freq, int num_nt,
double** freq_1, double* centroid, double* centroid_tilde, 
double **quality, double *expected_qual, double **quality_1,
double *expected_freq)
{
	memset(centroid, 0, row_length * sizeof(*centroid));
	for(int n=0; n<N; n++){
		for(int l=0; l<row_length; l++){
			centroid[l] += quality[n][l];
		}
	}

	for(int l=0; l<row_length; l++){
		centroid[l] /= N;
	}
	return;
}


// d2 centroid
void d2_centroid(int N, int row_length, double** freq, int num_nt,
double** freq_1, double* centroid, double* centroid_tilde, 
double **quality, double *expected_qual, double **quality_1,
double *expected_freq)
{
	// add the total counts
	memset(centroid, 0, row_length * sizeof(*centroid));
	for(int n=0; n<N; ++n){
		for(int l=0; l<row_length; ++l){
			centroid[l] += quality[n][l];
		}
	}
	// now normalize the entries
	double total_sq_count = 0;
	for(int l=0; l<row_length; ++l){
		total_sq_count += centroid[l] * centroid[l];
	}
	double norm = sqrt(total_sq_count);
	for(int l=0; l<row_length; ++l){
		centroid[l] /= norm;
	}
	return;
}


// KL centroid --- operates on raw counts which do not have to be normalized
void kl_centroid(int N, int row_length, double** freq, int num_nt,
double** freq_1, double* centroid, double* centroid_tilde, 
double **quality, double *expected_qual, double **quality_1,
double *expected_freq)
{
	// add the total counts
	memset(centroid, 0, row_length * sizeof(*centroid));
	for(int n=0; n<N; ++n){
		for(int l=0; l<row_length; ++l){
			centroid[l] += quality[n][l];
		}
	}
	// now normalize the entries
	double total_count = 0;
	for(int l=0; l<row_length; ++l){
		total_count += centroid[l];
	}
	for(int l=0; l<row_length; ++l){
		centroid[l] /= total_count;
	}
	return;
}


// MM centroid --- evaluates frequencies of single nucleotides
// and coumputes frequencies of words using zero order Markov model;
// i.e., as product of single nucleotide frequencies:
// f_{w1...w_k} = f_w1 f_w2 ... f_wk
void mm_centroid(int N, int row_length, double** freq, const int num_nt, 
		double** freq_1, double* centroid, double* centroid_tilde, 
		double **quality, double *expected_qual, double **quality_1,
		double *expected_freq)
{
	//If p_method is global (P1G or P2G) we don't need to calculate centroids
	if (expected_freq != NULL) return;
		
	//P2L: It uses a markovian model to compute the expected frequency of words
	double nt_freq[num_nt];
	memset(nt_freq, 0, num_nt * sizeof(nt_freq[0]));
	double S=0;
	// compute frequencies of individual nucleotides
	for (int n=0; n<N; n++){
		for (int l=0; l<num_nt; l++){
			//Number of occurrencies of LETTERS in the cluster:
			nt_freq[l] += freq_1[n][l]; 
		}
	}
	for (int l=0; l<num_nt; l++)
		S += nt_freq[l];	//S = number of letter in the cluster

	for (int l=0; l<num_nt; l++){
		nt_freq[l] /= S;	//Average frequency of letters in the cluster
	}
	// compute frequencies of words using Markov model
	for(int l=0; l<row_length; l++){
		double p = 1;
		for(int M=1; M<row_length; M*=num_nt){
			int digit = l/M % num_nt;
			p *= nt_freq[digit];
		}
		centroid[l] = p;
	}
	return;
	
	/* P1L: Do not use: produces extremely unbalanced clustering*/
	/*int tot_word = 0;
	for (int l=0; l<row_length; l++) {
		centroid[l] = 0;
		for (int n=0; n<N; n++)
			centroid[l] += freq[n][l];
		tot_word += centroid[l];
	}
	
	//Minimum value: 0.01, to avoid division per zero later
	for (int l=0; l<row_length; l++) {
		if (centroid[l] == 0){
			centroid[l] = 0.01;			
		}
		centroid[l] /= tot_word;
	}*/
}


// d2* centroid
void d2ast_centroid(int N, int row_length, double** freq, int num_nt,
double** freq_1, double* centroid, double* centroid_tilde, 
double **quality, double *expected_qual, double **quality_1,
double *expected_freq)
{
	// compute frequencies from MM
	mm_centroid(N, row_length, freq, num_nt, freq_1, centroid, NULL, 
				quality, expected_qual, quality_1, expected_freq);
	// compute X_tilde
	memset(centroid_tilde, 0, row_length * sizeof(*centroid_tilde));
	
	double total_count = 0;
	for(int n=0; n<N; ++n){
		for(int l=0; l<row_length; ++l){
			centroid_tilde[l] += quality[n][l];
			total_count += freq[n][l];
		}
	}
	
	//Determines what data to use: global(expected_freq) or local(centroid)
	double * exp_freq;
	if (expected_freq==NULL) exp_freq = centroid;
	else exp_freq = expected_freq;
	
	double S = 0;
	for(int l=0; l<row_length; ++l){
		centroid_tilde[l] -= total_count * exp_freq[l] * expected_qual[l];
		centroid_tilde[l] /= sqrt(exp_freq[l] * expected_qual[l]);
		S += centroid_tilde[l] * centroid_tilde[l];
	}
	S = sqrt(S);
	for(int l=0; l<row_length; ++l){
		centroid_tilde[l] /= S;
	}
	return;
}
