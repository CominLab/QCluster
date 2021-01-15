#pragma once

extern int max_iterations;

/**
 * Hard EM clustering
 * K = number of clusters
 * N = number of samples
 * row_length = row length of the freq matrix
 * freq = frequencies  matrix
 * quality = quality matrix
 * expected_qual = vector of expected quality of each word
 * quality_1 = quality matrix for single bases
 * expected_freq = vector of expected frequencies of each word
 * num_nt = number of different bases (4)
 * freq_1 = frequencies matrix for single bases
 * assignment = vector of cluster assignment; must be pre-allocated
 * Z = flat array for storing assignment confidence; must be pre-allocated
 * num_trials = number of times to repeat clustering before chosing the result
 * 	with the minimal distortion
 * dist = distance: e = Euclidean, k = KL, s = symmetrized KL
 * verbose = verbosity level
 * Returns the number of clusters formed
 */
int hard_em(int K, int N, int row_length, double** freq, 
double **quality, double *expected_qual, double **quality_1, double *expected_freq, int num_nt,
double** freq_1, int* assignment, double* Z, int num_trials=1,
char dist_type='e', int verbose=0);


