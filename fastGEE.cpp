// Biostats 615 - Statistical Computing //
// fastGEE.cpp - Tahmeed Tureen, Henghsi Yu, and Weiqi Zhou
// Supplement C++ Code for Final Project

// Load typical packages
#include<cmath>
#include<iostream>
#include<iomanip>

// Typiccal C libraries
#include<cmath>
#include<cstdlib>

// C++ Data Structures
#include<set>
#include<map>

// Include Eigen Package for Matrices
#include <Eigen/Core>

// Load Matrix615 header file to read in from file
#include "Matrix615.h"

using namespace std;
// we avoid using namespace Eigen to be able to clearly see where we are calling the Eigen Package

// Requires: A Eigen Matrix
// Modifies: Nothing
// Effects: Returns the number of unique observations (patients, study participants) in the data set

template <typename any>
int count_unique(const Eigen::EigenBase<any> &mat_in) {
	std::set<double> unique_set;
	for (int i = 0; i < int(mat_in.rows()); ++i) {
		unique_set.insert( mat_in(i, 0) );
	}

	return int(unique_set.size());
}


// Driver Code
int main(int argc, char* argv[]) {

	// Steps corresponding to the commented code is associated with the following documents:
	
	// Tentative Read Matrix using Matrix615.h
	
	Matrix615<double> read_Design; // throw-away Matrices to read in data
	Matrix615<double> read_Y;
	Matrix615<double> read_ID;

	// Read the data in
	read_Design.readFromFile(argv[3]); // t
	read_Y.readFromFile(argv[2]);
	read_ID.readFromFile(argv[1]);

	// Step (1): Store the design matrix, response, and Id into an Eigen Matrix
	Eigen::MatrixXd	design_X;
	Eigen::MatrixXd outcome_Y;

	read_Design.cloneToEigen(design_X);
	read_Y.cloneToEigen(outcome_Y);

	// cout << design_X.rows() << " " << design_X.cols() << endl << endl;
	// cout << design_X << endl << endl << endl;
	// cout << outcome_Y.rows() << " " << outcome_Y.cols() << endl << endl;
	// cout << outcome_Y << endl << endl << endl;

	// Step (2): Count number of unique ID's, patients
	Eigen::MatrixXd ID_Mat;
	read_ID.cloneToEigen(ID_Mat);

	std::set<int> unique_set; // to dynamically count n
	std::map<int, int> id_map; // to dynamically count m

	for (int i = 0; i < int(outcome_Y.rows()); ++i) {
		unique_set.insert(  ID_Mat(i,0) );
		id_map[ID_Mat(i, 0)] += 1; 
	}

	int n = int( unique_set.size() );
	int m = int( id_map[ID_Mat(1,0)] );

	cout << "n is " << n << endl << endl;
	cout << "m is " << m << endl << endl;

	// Step (3): Initialize betas *********** Should we incorporate the beta's that come from a fast linear regression or no?

	Eigen::MatrixXd Betas_new(int(design_X.cols()) , 1); // dimensions (p+1)  x 1
	Betas_new.setZero();

	cout << "Betas ***********************" << endl;
	cout << Betas_new << endl << endl;

	// Step (4): Initialize rho and phi (rho: off diagonals of correlation structure & phi: )
	double rho = 0.0;
	double phi = 0.0;

	// Step (5): Set up variables for iteration convergence check

	// Should we make this dynamic? 

	double diff_beta = 1; // updated difference in beta estimation {beta (new) - beta (old)}
	
	double diff_threshold = 0.00000010; // threshold for the update difference in betas
	
	int iteration_threshold = 1000; // threshold on how many iterations there will be

	// Step (6): Initialize the iteration
	int iteration_count = 0;

	// Step (7): Assign appropriate value to n*
	double n_star = (0.5) * double(n) * double(m) * double(m - 1);

	cout << " n* is " << n_star << endl;

	// Step (8): Start the GEE Estimation

	Eigen::MatrixXd Betas_updated = Betas_new;

	cout << Betas_updated << endl;

	while ( (diff_beta > diff_threshold) && (iteration_count < iteration_threshold) ) {
		// (1)
		Betas_updated = Betas_new;

		// (2) Vector EE ((p+1) x 1)
		Eigen::MatrixXd EE( int(design_X.cols()), 1);
		EE.setZero();

		// (3) Matrix GI ( (p+1) x (p+1) )
		// First derivative of EE; Model-based Variance
		Eigen::MatrixXd GI( int(design_X.cols()) , int(design_X.cols()));
		GI.setZero();

		// (4) Matrix G ( (p+1) x (p+1) )
		// Meat of the Sandwhich Estimator
		Eigen::MatrixXd G( int(design_X.cols()) , int(design_X.cols()));
		G.setZero();		

		// (5) Initialize phi(sum) and tau(sum)
		double phi_sum = 0.0;
		double tau_sum = 0.0;

		// inner loop in while loop Time-Complexity: O(n^2)??

		for (int i; i < n; ++i) {
			start = end + 1;
			end ;
		}

		break;
	}


	return 0;

}