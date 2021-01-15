/**
 * Distance functions to be used with EM clustering routine
 * Arguments are as follows: 
 * q = centroid coordinates, appropriately normalized
 * p = word counts, appropriately normalized
 * L = length of the two vectors p and q
 * centroid_tilde = vector of X_tilde values for centroid; only used in the calculation 
 * 		of d2* distance, otherwise this argument is irrelevant and kept to
 * 		preserve the call signature
 * Appropriate normalization (e. g., raw counts, frequencies, norm equal to 1)
 * 		depends on the distance used. It should be done elsewhere
 */

// squared euclidean distance
// word count vector needs to be normalized so that
// its components add to 1
double euclidean_distance(double *x, double *y, int L, double* xt, 
						  double *quality, double *expected_qual, double *expected_freq);

// KL divergence; expects raw counts for p and frequencies for q
double kl_distance(double *q, double *p, int L, double* centroid_tilde, 
				   double *quality, double *expected_qual, double *expected_freq);

// symmetrized KL divergence; expects raw counts for p and frequencies for q
double symkl_distance(double *q, double *p, int L, double* centroid_tilde, 
					  double *quality, double *expected_qual, double *expected_freq);

// d2 distance; expects any normalization for p and unit norm for q
double d2_distance(double *q, double *p, int L, double* centroid_tilde, 
				   double *quality, double *expected_qual, double *expected_freq);

// chi2 distance, inspired by d2* distance; expects raw counts for p and 
//		frequencies for q
double chi2_distance(double *q, double *p, int L, double* centroid_tilde, 
					 double *quality, double *expected_qual, double *expected_freq);

// d2* distance; expects any normalization for p, frequencies for q and unit
// 		norm for centroid_tilde
double d2ast_distance(double *q, double *p, int L, double* centroid_tilde, 
					  double *quality, double *expected_qual, double *expected_freq);
