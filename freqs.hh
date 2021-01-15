#pragma once

#include <string>
#include <vector>

using namespace std;

// number of characters
char const base ='!'; // Illumina 1.8+ Phred+33

enum {NUM_NT = 4};

// normalize row dividing by the total counts
void normalize_row(double* fv, int L);

void fill_overlap_count_vector(string seq, string qual, int K, 
		double* freq_vector, double *quality, int normalize=2, 
		double pc=1, double *avg_quality_1=NULL, double **avg_quality = NULL, 
		double *freq_1 = NULL, bool redistribute = true);

/**Function responsible for the calculation of the quality expected value for 
 * each k-word: the expected value is computed basing on the full dataset.
 * Refer to the documentation for the description of the possible methods. 
 * Default method: E1 */
void calculate_quality_expected_value(int method, int N, int K, int L, 
			double **freq, double **quality, double **freq_1,
			 double* avg_quality_1, double **avg_quality, double *expected_qual);
	
void expected_frequency_p2global(int N, int K, int L, double *expected_freq,
								double **freq_1);

void expected_frequency_p1global(int N, int K, int L, double *expected_freq,
								double **freq);

// normalize frequency matrix to make its columns univariant
void normalize_freq_matrix(double** data, double** qual, int N, int row_length);

