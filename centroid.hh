/**
 * Centroid evaluation functions to be used with EM clustering routine
 */

// mean centroid
void mean_centroid(int N, int row_length, double** freq, int num_nt,
double** freq_1, double* centroid, double* centroid_tilde, 
double **quality, double *expected_qual, double **quality_1,
double *expected_freq);

// d2 centroid
void d2_centroid(int N, int row_length, double** freq, int num_nt,
double** freq_1, double* centroid, double* centroid_tilde, 
double **quality, double *expected_qual, double **quality_1,
double *expected_freq);

// KL centroid --- operates on raw counts which do not need to be normalized
void kl_centroid(int N, int row_length, double** freq, int num_nt,
double** freq_1, double* centroid, double* centroid_tilde, 
double **quality, double *expected_qual, double **quality_1,
double *expected_freq);

// MM centroid --- evaluates frequencies of single nucleotides
// and coumputes frequencies of words using zero order Markov model
void mm_centroid(int N, int row_length, double** freq, const int num_nt, 
double** freq_1, double* centroid, double* centroid_tilde, 
double **quality, double *expected_qual, double **quality_1,
double *expected_freq);

// d2* centroid
void d2ast_centroid(int N, int row_length, double** freq, int num_nt,
double** freq_1, double* centroid, double* centroid_tilde, 
double **quality, double *expected_qual, double **quality_1,
double *expected_freq);
