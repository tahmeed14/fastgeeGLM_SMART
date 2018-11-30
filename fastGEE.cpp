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
#include<vector>

// Include Eigen Package for Matrices
#include <Eigen/Core>

// Load Matrix615 header file to read in from file
#include "Matrix615.h"

using namespace std;
// we avoid using namespace Eigen to be able to clearly see where we are calling the Eigen Package

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

	Eigen::MatrixXd Weights(outcome_Y.rows(), 1);
	Weights.setOnes();

	// cout << Weights.rows() << endl << endl;
	// cout << Weights << endl;


	if (int (argc) == 5) {
		Matrix615<double> read_Weights;
		read_Weights.readFromFile(argv[4]);
	} 

	// cout << design_X.rows() << " " << design_X.cols() << endl << endl;
	// cout << design_X << endl << endl << endl;
	// cout << outcome_Y.rows() << " " << outcome_Y.cols() << endl << endl;
	// cout << outcome_Y << endl << endl << endl;

	// Step (2): Count number of unique ID's, patients
	Eigen::MatrixXd ID_Mat;
	read_ID.cloneToEigen(ID_Mat);

	std::set<int> unique_ID; // to dynamically count n
	std::map<int, int> id_map; // to dynamically count m
	std::map<int, double> weight_map;

	for (int i = 0; i < int(ID_Mat.rows()); ++i) {
		unique_ID.insert(  ID_Mat(i,0) );
		id_map[ID_Mat(i, 0)] += 1;
		weight_map[ID_Mat(i, 0)] = Weights(i, 0); 
	}

	int n = int( unique_ID.size() );
	// int m = int( id_map[ID_Mat(1,0)] );

	// q = (p + 1) : # number of paramters to be estimated using GEE
	int q = int ( design_X.cols() );

	// Correlation Structure

	// Step : Weights

	// Step (3): Initialize betas *********** Should we incorporate the beta's that come from a fast linear regression or no?

	Eigen::MatrixXd Betas_new(q , 1); // dimensions (p+1)  x 1
	Betas_new.setZero();

	// cout << "Betas ***********************" << endl;
	// cout << Betas_new << endl << endl;

	// Step (4): Initialize rho and phi (rho: off diagonals of correlation structure & phi: )
	double rho = 0.0;
	// double phi = 0.0;

	// Step (5): Set up variables for iteration convergence check

	// Should we make this dynamic? 

	double diff_beta = 1; // updated difference in beta estimation {beta (new) - beta (old)}
	
	double diff_threshold = 0.00000010; // threshold for the update difference in betas
	
	int iteration_threshold = 1000; // threshold on how many iterations there will be

	// Step (6): Initialize the iteration
	int iteration_count = 0;

	// Step (7): Assign appropriate value to n*
	// double n_star = (0.5) * double(n) * double(m) * double(m - 1);

	// double n_star = (0.5) * double(n) * double(5) * double(5 - 1);

	// cout << " n* is " << n_star << endl;

	double n_star_sum = 0.0;

	// Step (8): Start the GEE Estimation

	Eigen::MatrixXd Betas_updated = Betas_new;

	Eigen::MatrixXd sandwhich_Mat;

	// cout << Betas_updated << endl << endl;

	while ( (diff_beta > diff_threshold) && (iteration_count < iteration_threshold) ) {
		// (1)
		Betas_updated = Betas_new;

		// (2) Vector EE ((p+1) x 1)
		Eigen::MatrixXd EE(q , 1);
		EE.setZero();

		// (3) Matrix GI ( (p+1) x (p+1) )
		// First derivative of EE; Model-based Variance
		Eigen::MatrixXd GI(q , q);
		GI.setZero();

		// (4) Matrix G ( (p+1) x (p+1) )
		// Meat of the Sandwhich Estimator
		Eigen::MatrixXd G(q , q);
		G.setZero();		

		// (5) Initialize phi(sum) and tau(sum)
		double phi_sum = 0.0;
		double tau_sum = 0.0;

		// (6) Initialize start and end
		int start = 0;
		int end = -1;

		// inner loop in while loop Time-Complexity: O(n^2)??

		for (int i = 0; i < n; ++i) {

			int m = id_map[i + 1];
			start = end + 1;
			end = start + m - 1;

			//update n*
			n_star_sum = n_star_sum + (0.5) * double(m) * double(double(m) - 1.0);

			// assign mu_i ( m x (p+1) ) for ith person
			// Eigen.block(starting_row = , starting_col = , dim_row = , dim_col = )
			Eigen::MatrixXd mu_i = design_X.block(start, 0, m, q) * Betas_updated;

			// assign r_i (m x 1)
			Eigen::MatrixXd r_i = outcome_Y.block(start, 0, m, 1) - mu_i;

			// Loop to update phi-sum
			for (int j = 0; j < int(m); ++j) phi_sum = phi_sum + (r_i(j,0) * r_i(j,0));

			// Loop to update tau-sum
			for (int j = 0; j < int(m - 1); ++j) {
				for (int k = (j + 1); k < m; ++k) tau_sum = tau_sum + r_i(j,0) * r_i(k, 0);
			}

			// create R matrix (m x m) for each observation/patient
			Eigen::MatrixXd R(m, m);
			R.setConstant(rho); // off-diagonals will be rho
			for (int d = 0; d < int(m); ++d) R(d, d) = 1; //diagonals will be 1

			// update EE ((p+1) x 1)
			EE = EE + ((design_X.block(start, 0, m, q).transpose()) * (R.inverse() * (weight_map[i+1] * r_i) )); // weight_map

			// update GI ((p+1) x (p+1)) or (q x q)
			GI = GI + ((design_X.block(start, 0, m, q).transpose()) * R.inverse() *  design_X.block(start, 0, m, q));

			// update G ((p+1) x (p+1)) or (q x q)
			G = G + ((design_X.block(start, 0, m, q).transpose()) *
				( ((R.inverse()) * r_i) * (r_i.transpose()*(R.inverse())) ) * design_X.block(start, 0, m, q));
		}

		// Update beta using either Newton Raphson or Gradient Boosting *** we can add an if statement
		Betas_new = Betas_updated + ( GI.inverse() * EE);

		// calculate the difference in the betas
		diff_beta = (Betas_new - Betas_updated).norm();

		// update rho
		rho = ( (( double(n) - double(q) ) * tau_sum) / ( (double(n_star_sum) - double(q) ) * phi_sum) );


		sandwhich_Mat = GI.inverse() * G * GI.inverse();
		// increment iteration count
		iteration_count += 1;
	}

	cout << endl << "Estimated Betas are: " << endl;
	cout << Betas_new.transpose() << endl << endl << endl;

	// cout << "The sandwhich matrix is: " << endl << endl;
	// cout << sandwhich_Mat << endl << endl;

	cout << "The beta S.E.'s are: " << endl << endl;
	cout << sandwhich_Mat.diagonal().array().sqrt() << endl << endl;
	
	return 0;

}