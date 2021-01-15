#include <string.h> // memcpy
#include <math.h>
#include<iostream>
using namespace std;



double euclidean_distance(double *centroid, double *y, int L, double* xt, 
						  double *quality, double *expected_qual, double *expected_freq)
{
	double S = 0;
	double z;
	for(int i=0; i<L; i++){
		z = centroid[i] - quality[i];
		S += z*z;
	}
	return S;
}


double kl_distance(double *centroid, double *p, int L, double* centroid_tilde, 
				   double *quality, double *expected_qual, double *expected_freq)
{
	double total_count = 0;
	for(int i=0; i<L; ++i){
		total_count += quality[i];
	}
	double S = 0;
	for(int i=0; i<L; i++){
		S += quality[i] * log(quality[i]/(total_count * centroid[i]));
	}
	return S;
}


double symkl_distance(double *centroid, double *p, int L, double* centroid_tilde, 
					  double *quality, double *expected_qual, double *expected_freq)
{
	double total_count = 0;
	for(int i=0; i<L; ++i){
		total_count += quality[i];
	}
	double S = 0;
	for(int i=0; i<L; i++){
		S += (quality[i]-total_count*centroid[i]) * log(quality[i]/(total_count*centroid[i]));
	}
	return S;
}


double d2_distance(double *centroid, double *p, int L, double* centroid_tilde, 
				   double *quality, double *expected_qual, double *expected_freq)
{
	double norm = 0;
	double S = 0;
	for(int i=0; i<L; i++){
		norm += quality[i] * quality[i];
	}
	norm = sqrt(norm);
	for(int i=0; i<L; i++){
		S += centroid[i] * quality[i] / norm;
	}
	return 1 - S;
}


double chi2_distance(double *centroid, double *p, int L, double* centroid_tilde, 
					 double *quality, double *expected_qual, double *expected_freq)
{
	double chi2 = 0;
	double total_count = 0;
	for(int i=0; i<L; i++){
		total_count += quality[i];
	}
	for(int i=0; i<L; i++){
		double exp_count = centroid[i] * total_count;
		chi2 += (quality[i] - exp_count) * (quality[i] - exp_count) / exp_count;
	}
	//return (chi2 - (L-3) * log(chi2))/2;
	return chi2;
}

double d2ast_distance(double *centroid, double *p, int L, double* centroid_tilde,
					  double *quality, double *expected_qual, double *expected_freq)
{
	//Determines what data to use: global(expected_freq) or local(centroid)
	double * exp_freq;
	if (expected_freq==NULL) exp_freq = centroid;
	else exp_freq = expected_freq;
	
	double x[L];
	memcpy(x, quality, L*sizeof(*quality));
	double total_count = 0;
	for(int i=0; i<L; i++){
		total_count += p[i];
	}
	double S = 0;
	for(int l=0; l<L; ++l){
		x[l] -= total_count * exp_freq[l] * expected_qual[l];
		x[l] /= sqrt(exp_freq[l] * expected_qual[l]);
		S += x[l] * x[l];
	}
	S = sqrt(S);
	for(int l=0; l<L; ++l){
		x[l] /= S;
	}	
	S = 0;
	for(int l=0; l<L; ++l){
		S += x[l] * centroid_tilde[l];
	}
	return (1-S)/2;
}

