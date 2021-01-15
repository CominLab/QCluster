#include "em.hh"
#include "freqs.hh"
#include "dists.hh"
#include "centroid.hh"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

using namespace std;


int *best_assignment;
double best_distortion;
int iterations;
int max_iterations = 1000;


/**Initialize the vector that mantains the best assignment found so far. */
void initialize_cache(int N){
	iterations = 0;
	best_assignment = new int[N];
	best_distortion = 1000000000000000;//infinite	
	for (int i=0; i<N; i++) best_assignment[i] = -1;
	return;
}

/**Function that controls if the number of iterations without improvement of 
 * the best solution exceed a bound. If so, it substitutes the current solution
 * with the best found so far and return true. */
bool exceed_max_iterations(int N, int K, int *assignment, double &distortion){
	iterations++;
	if (iterations>max_iterations){
		distortion = best_distortion;
		for(int n=0; n<N; n++) assignment[n] = best_assignment[n];	
		return true;
	}
	else{
		if (distortion < best_distortion){
			iterations = 0;
			for(int n=0; n<N; n++) best_assignment[n] = assignment[n];	
			best_distortion = distortion;
		}
	}	
	return false;
}


void deallocate_cache(){
	delete [] best_assignment;
	return;
}



/**
 * Evaluate confidence of assignment to each centroid
 */
static void eval_confidence(int K, int N, int row_length, double **data,
double **centroids, double *Z, double **quality, double *expected_qual, double *expected_freq)
{
	static const double KL_DIST_CUTOFF = 10; // set exponential to zero
		// if we exceed this negative exponent
	for (int n=0; n<N; n++){
		double* row = Z + n*K;
		for(int k=0; k<K; k++){
			row[k] = kl_distance(centroids[k], data[n], row_length, NULL, 
								 quality[n], expected_qual, expected_freq);
		}
		double min_dist = *row;
		for(int k=1; k<K; k++){
			if (min_dist > row[k]){
				min_dist = row[k];
			}
		}
		double S = 0;
		for(int k=0; k<K; k++){
			row[k] -= min_dist;
			S += row[k] = row[k] > KL_DIST_CUTOFF ? 0 : exp(-row[k]);
		}
		for(int k=0; k<K; k++){
			row[k] /= S;
		}
	}
	return;
}


// Auxiliary function to select procedures for evaluation of 
// distances and centroids
void select_dists_cent(char dist_type, 
						double (**distf)(double*, double*, int,
						double*, double*, double*, double *), 
						void (**eval_centroid)(int, int, double**,
						int, double**, double*, double*, double **, double *, 
						double**, double*))
{
	switch (dist_type){
		case 'a':
			*distf = &d2ast_distance;
			*eval_centroid = &d2ast_centroid;
			break;
		case 'c':
			*distf = &chi2_distance;
			*eval_centroid = &mm_centroid;
			break;
		case 'd':
			*distf = &d2_distance;
			*eval_centroid = &d2_centroid;
			break;
		case 'k':
			*distf = &kl_distance;
			*eval_centroid = &kl_centroid;
			break;
		case 's':
			*distf = &symkl_distance;
			*eval_centroid = &kl_centroid;
			break;			
		case 'e': // fall through
		default:
			*distf = &euclidean_distance;
			*eval_centroid = &mean_centroid;
	}

	return;
}


/** EM clustering routine
 * This is an auxiliary routine, and it assumes that the memory for the 
 * necessary computations is already pre-allocated. This routine is called
 * from other routines. Calling routines can do a single clustering run
 *  or multiple runs with consensus clustering.
 * K = number of clusters
 * N = number of samples
 * row_length = row length
 * freq = frequencies matrix
 * quality = quality matrix
 * expected_qual = vector of expected quality of each word
 * quality_1 = quality matrix for single bases
 * expected_freq = vector of expected frequencies of each word
 * num_nt = number of different bases (4)
 * freq_1 = frequencies matrix for single bases
 * assignment = assignment of samples to clusters
 * numMembers = number of members in each cluster
 * centroids = matrix for centroids
 * tmp_data = array of double* of length N for storing elements in each cluster;
 * 		used in centroid calculation
 * distf = function for evaluating distance
 * eval_centroid = function for evaluating centroids
 * 
 * Return value is the distortion (sum of the distances from each data 
 * point to its corresponding centroid; in the case of the squared 
 * Euclidean distance it gives intra-cluster variance).
 */
static double em_routine(int K, int N, int row_length, double** freq, 
		double **quality, double *expected_qual, double **quality_1, double *expected_freq, 
		int num_nt, double** freq_1, int* assignment, int* numMembers,
		double** centroids, double** centroids_tilde, double** tmp_data,
		double**tmp_data_1, char dist_type,
		int verbose=0, int allow_empty_clusters = 0)
{
	double **tmp_qfreq = new double*[N];
	double ** tmp_qfreq_1 = new double*[N];

	
	double distortion = 0;
	initialize_cache(N);
	// choose auxiliary functions
	double (*distf)(double*, double*, int, double*, double *, double*, double *);
	void (*eval_centroid)(int, int, double**, int, double**, double*, 
						  double*, double **, double *, double**, double*);
	select_dists_cent(dist_type, &distf, &eval_centroid);

	while (true){ // repeat clustering attempts till we get a result
restart:
		if (verbose > 0) {
			cerr<<"\tEM routine: starting clustering attempt!"<<endl;
		}
		
		memset(assignment, 0, N*sizeof(*assignment));
		memset(numMembers, -1, K*sizeof(*numMembers));
		// initial centroid assignment --- initialize with randomly chosen 
		// freq items
		for (int k=0; k<K; ++k) {
			int point_number = rand()% N;
			eval_centroid(1, row_length, freq + point_number, num_nt, freq_1 +
			point_number, centroids[k], centroids_tilde[k], 
			quality + point_number, expected_qual, 
			quality_1 + point_number, expected_freq);
		}
		
		
		// start iterations of EM clustering
		while (true){
			// compute the distances to the new centroids 
			// and compute the new assignments
			bool assignmentChanged = false;  // set to true if 
			// at least one element changes its cluster assignment
			for(int n=0; n<N; n++){ // update assignment of each element
				int new_assignment=-1;  // new assignment
				double min_dist=-1;  // distance to the closest centroid
				// initial value set here is irrelevant as it 
				// gets changed later
				for(int k=0; k<K; k++){ // loop through centroids to 
				// find the closest one
					if (!allow_empty_clusters || numMembers[k]){
						// evaluate the distance; note the order of 
						// arguments: centroid followed by freq --- it's 
						// important for non-symmetric distances
						double dist = (*distf)(centroids[k], freq[n], 
							row_length, centroids_tilde[k], quality[n], 
							expected_qual, expected_freq);
						if (new_assignment==-1 || dist < min_dist) {
							min_dist = dist;
							new_assignment = k;
						}
					}
				}
				// check for assignment change
				assignmentChanged = assignmentChanged || (assignment[n] != 
					new_assignment);
				// update assignment
				assignment[n] = new_assignment;
			}
			
			// evaluate the number of members in each cluster
			memset(numMembers, 0, K*sizeof(*numMembers));
			for(int n=0; n<N; n++){
				numMembers[assignment[n]] += 1;
			}

			// check for empty clusters
			if (!allow_empty_clusters){
				for(int k=0; k<K; k++){
					if(!numMembers[k]){
						if (verbose > 0) {
							cerr<<"\tEM routine: empty cluster"<<endl;
						}
						// got an empty cluster; continue from the beginning
						deallocate_cache();
						initialize_cache(N);
						goto restart;
					}
				}
			}
			
			distortion = 0;
			for (int n=0; n<N; n++) {
				distortion += distf(centroids[assignment[n]], freq[n],
				row_length, centroids_tilde[assignment[n]], quality[n], 
				expected_qual, expected_freq);
			}
			if (!assignmentChanged){
				// clustering succeded; return
				if (verbose > 0) {
					cerr<<"\tEM routine: attempt succeded"<<endl;
				}
				deallocate_cache();
				delete [] tmp_qfreq;
				delete [] tmp_qfreq_1;
				return distortion;
			}
			
			//If the number of iterations has excedeed a maximum, the 
			//algorithm stops and the best solution found is returned.
			if (exceed_max_iterations(N, K, assignment, distortion)){
				if (verbose > 0) {
					cerr<<"\tEM routine: Cycling. Returning"<<endl;
				}
				deallocate_cache();
				delete [] tmp_qfreq;
				delete [] tmp_qfreq_1;
				return distortion;
			}			

			// evaluate the new centroids
			for(int k=0; k<K; k++){
				if(!allow_empty_clusters || numMembers[k]){
					int elems_in_cluster = 0;
					for(int n=0; n<N; ++n){
						if(assignment[n] == k){
							tmp_data[elems_in_cluster] = freq[n];
							tmp_data_1[elems_in_cluster] = freq_1[n];
							
							tmp_qfreq[elems_in_cluster] = quality[n];
							tmp_qfreq_1[elems_in_cluster] = quality_1[n];
							
							elems_in_cluster++;
						}
					}
					eval_centroid(elems_in_cluster, row_length, tmp_data,
						num_nt, tmp_data_1, centroids[k], centroids_tilde[k],
						tmp_qfreq, expected_qual, tmp_qfreq_1,expected_freq);
				}
			}
		}
	}
}


// Count number of distinct clusters from assignment vector.
// Cluster numbers are non-negative. They can be non-contiguous:
// e. g., 0, 2, 4, 5, 7
static int count_num_clusters(int const N, const int * assignment)
{
	int max_id = 0;
	for (int n=0; n<N; n++){
		if (assignment[n] > max_id){
			max_id = assignment[n];
		}
	}
	max_id++; // increment to get the size of the index vector
	int *ix = new int[max_id];
	memset(ix, 0, max_id*sizeof(*ix));
	for (int n=0; n<N; n++) {
		ix[assignment[n]] = 1;
	}
	int  num_clusters=0;
	for (int k=0; k<max_id; k++) {
		num_clusters += ix[k];
	}
	delete[] ix;

	return num_clusters; 
}


// implementation of a publicly accessible function
int hard_em(int K, int N, int row_length, double** data, 
		double **quality, double *expected_qual, double **quality_1, double *expected_freq,
		int num_nt, double** data_1, int* assignment, double *Z, 
		int num_trials, char dist_type, int verbose)
{
	// Allocate matrices for centroids, distances, assignment and the 
	// vector for the number of members
	double**  centroids = new double*[K];
	centroids[0] = new double[row_length*K];
	for(int k=1; k<K; ++k) {
		centroids[k] = centroids[0] + k*row_length;
	} // K*row_length matrix for centroid locations

	double**  centroids_tilde = new double*[K];
	centroids_tilde[0] = new double[row_length*K];
	for(int k=1; k<K; ++k) {
		centroids_tilde[k] = centroids_tilde[0] + k*row_length;
	} // K*row_length matrix for X_tilde vars for centroids

	int* tmp_assignment = new int[N];
	int *numMembers = new int[K];  // number of members in each cluster
	double** tmp_data = new double*[N];
	double** tmp_data_1 = new double*[N];

	// call clustering routine
	double min_distortion = 0;
	for (int t=0; t<num_trials; t++) {
		if (verbose > 0) {
			cerr<<"Calling EM routine, attempt "<<t+1<<endl;
		}
		double distortion =	em_routine(K, N, row_length, data, quality, 
				expected_qual, quality_1, expected_freq, num_nt, data_1,
				tmp_assignment, numMembers, centroids, centroids_tilde,
				tmp_data, tmp_data_1, dist_type, verbose);
		if (verbose > 0) {
			cerr<<"Resulting distortion: "<<distortion<<endl;
		}
		if (!t || (distortion < min_distortion)) {
			if (verbose > 0) {
				cerr<<"Updating best partitioning and minimal distortion"<<endl;
			}
			min_distortion = distortion;
			memcpy(assignment, tmp_assignment, N*sizeof(*assignment));
		}
	}
	if (verbose > 0) {
		cerr<<"Resulting minimal distortion: "<<min_distortion<<endl;
	}
	//evaluate confidence if using KL with raw counts
	if (dist_type == 'k') {
		eval_confidence(K, N, row_length, data, centroids, Z, quality, expected_qual, expected_freq);
	}
	// delete allocated arrays
	delete[] centroids[0];
	delete[] centroids;
	delete[] centroids_tilde[0];
	delete[] centroids_tilde;
	delete[] numMembers;
	delete[] tmp_data;
	delete[] tmp_data_1;
	delete[] tmp_assignment;

	// return the number of clusters
	return count_num_clusters(N, assignment);
}


