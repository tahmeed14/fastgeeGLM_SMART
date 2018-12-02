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

// To calculate Time
#include<ctime>

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

	// cout << Weights.rows() << endl << endl;
	// cout << Weights << endl;
	Eigen::MatrixXd Weights(outcome_Y.rows(), 1);
	// Weights.setOnes();

	// cout << "argc is " << argc << endl;

	if (int (argc) == 5) {
		Matrix615<double> read_Weights;
		read_Weights.readFromFile(argv[4]);

		Eigen::MatrixXd Weights;
		read_Weights.cloneToEigen(Weights);
	}

	else {
		Weights.setOnes();
	}


	// Step (2): Count number of unique ID's, patients
	Eigen::MatrixXd ID_Mat;
	read_ID.cloneToEigen(ID_Mat);

	std::set<int> unique_ID; // to dynamically count n
	std::map<int, int> m_map; // to dynamically count m
	std::map<int, double> weight_map;

	for (int i = 0; i < int(ID_Mat.rows()); ++i) {
		unique_ID.insert(  ID_Mat(i,0) );
		m_map[ID_Mat(i, 0)] += 1;
		weight_map[ID_Mat(i, 0)] = Weights(i, 0); 
	}

	// int n = int( unique_ID.size() );
	
	// cout << "n is << " << n << endl;

	// cout << "size of the id to m map is " << m_map.size() << endl;

	// cout << "m of 1st patient is " << m_map[1] << " " << m_map[100] << " " << m_map[2] << endl << endl;
	// cout << "size of weight map " << weight_map.size() << endl;

	// q = (p + 1) : # number of paramters to be estimated using GEE
	int q = int ( design_X.cols() );

	// cout << "q is " << q << endl;

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
	double n_star_sum = 0.0;
	double m_sum = 0.0;

	// calculate n_star_sum

	// for (std::set<int>::iterator iter = unique_ID.begin(); iter != unique_ID.end(); ++iter) {
	// 	n_star_sum = n_star_sum + (0.5) * double(m_map[*iter]) * double(double(m_map[*iter]) - double(1.0));
	// }

	// Step (8): Start the GEE Estimation

	Eigen::MatrixXd Betas_updated = Betas_new;

	Eigen::MatrixXd sandwhich_Mat;

	Eigen::MatrixXd tempGI(q, q);

	// // cout << Betas_updated << endl << endl;

	int before_EE = clock();

	// Step (10)
	while ( (diff_beta > diff_threshold) && (iteration_count < iteration_threshold) ) {
		// (1)
		Betas_updated = Betas_new;

		// (2) Vector EE ((px1) x 1) or (q x 1)
		Eigen::MatrixXd EE(q, 1);
		EE.setZero();

		// (3) Matrix GI (q x q)
		Eigen::MatrixXd GI(q, q);
		GI.setZero();

		// // (4) Matrix G ( (p+1) x (p+1) )
		// // Meat of the Sandwhich Estimator
		// Eigen::MatrixXd G(q , q);
		// G.setZero();

		// (5) Initialize phi(sum) and tau(sum)
		long double phi_sum = 0.0;
		long double tau_sum = 0.0;

		// (6) Initialize indexes for block matrix multiplications
		int start = 0;
		int end = -1;

		// inner loop, we will loop using the unique_ID set (size = n)
		// Therefore, it loops n times
		for (std::set<int>::iterator iter = unique_ID.begin(); iter != unique_ID.end(); ++iter) {

			int m = m_map[*iter];
			start = end + 1;
			end = start + m - 1;

			n_star_sum = n_star_sum + (0.5) * double(m) * double(double(m) - 1.0);
			m_sum = m_sum + double(m);

			// cout << "start is " << start << endl;
			// cout << "end is " << end << endl << endl;

			// assign mu for ith observation
			// Eigen.block(starting_row = , starting_col = , dim_row = , dim_col = )
			Eigen::MatrixXd mu_i = design_X.block(start, 0, m, q) * Betas_updated;

			// cout << "mu_i is here " << endl;
			// cout << mu_i << endl << endl;

			// assign r_i (mi x 1)
			Eigen::MatrixXd r_i = outcome_Y.block(start, 0, m, 1) - mu_i;

			// cout << "r_i is here " << endl;
			// cout << r_i << endl;

			// add r_i ^2 to phi_sum from the ith observation
			for (int j = 0; j < int(m); ++j) {
				phi_sum = phi_sum + (r_i(j,0) * r_i(j,0));
			}

			// cout << "phi_sum is " << phi_sum << endl; 

			// add to tau_sum
			for (int j = 0; j < ( int(m) - 1 ); ++j) {
				for (int k = (j + 1); k < int(m); ++k) {
					// cout << " j is " << j << " and k is " << k << endl;
					tau_sum = tau_sum + (r_i(j,0) * r_i(k,0));
				}
			}

			// cout << "tau_sum is " << tau_sum << endl;

			// create R matrix (mi x mi) for each observation
			Eigen::MatrixXd R(m, m);
			R.setConstant(rho);
			for (int d = 0; d < int(m); ++d) {
				R(d, d) = 1;
			}

			// cout << R << endl << endl << endl;
			// update EE ((p+1) x 1)
			EE = EE + ((design_X.block(start, 0, m, q).transpose()) * (R.inverse() * (r_i) )); // weight_map

			// update GI ((p+1) x (p+1)) or (q x q)
			GI = GI + ((design_X.block(start, 0, m, q).transpose()) * R.inverse() *  design_X.block(start, 0, m, q));

			// update G ((p+1) x (p+1)) or (q x q)
			// G = G + ((design_X.block(start, 0, m, q).transpose()) * ( ((R.inverse()) * r_i) * (r_i.transpose()*(R.inverse())) ) * design_X.block(start, 0, m, q));

		}

		// Update beta using Newton Raphson method
		Betas_new = Betas_updated + (GI.inverse() * EE);

		// calculate the difference in the betas
		diff_beta = (Betas_new - Betas_updated).norm();
		// diff_beta = (Betas_new - Betas_updated).array().abs().sum();

		// update rho
		// cout << "tau_sum is " << tau_sum << endl;
		// cout << "phi_sum is " << phi_sum << endl;

		rho = ( (( double(m_sum) - double(q) ) * tau_sum) / ( (double(n_star_sum) - double(q) ) * phi_sum) );

		// cout << "rho is " << rho << endl << endl;

		// sandwhich_Mat = GI.inverse() * G * GI.inverse();
		tempGI = GI;

		iteration_count += 1;

	}

	// (4) Matrix G ( (p+1) x (p+1) )
	// Meat of the Sandwhich Estimator
	Eigen::MatrixXd G(q , q);
	G.setZero();

	// To calculate the sandwhich estimator
	// for (std::set<int>::iterator iter = unique_ID.begin(); iter != unique_ID.end(); ++iter) {
	// 	n_star_sum = n_star_sum + (0.5) * double(m_map[*iter]) * double(double(m_map[*iter]) - double(1.0));
	// }

	int start = 0;
	int end = -1;

	for (std::set<int>::iterator iter = unique_ID.begin(); iter != unique_ID.end(); ++iter) {
		int m = m_map[*iter];
		start = end + 1;
		end = start + m - 1;

		Eigen::MatrixXd mu_i = design_X.block(start, 0, m, q) * Betas_updated;

		// cout << "mu_i is here " << endl;
		// cout << mu_i << endl << endl;

		// assign r_i (mi x 1)
		Eigen::MatrixXd r_i = outcome_Y.block(start, 0, m, 1) - mu_i;

		Eigen::MatrixXd R(m, m);

		R.setConstant(rho);

		for (int d = 0; d < int(m); ++d) {
			R(d, d) = 1;
		}

		// update G ((p+1) x (p+1)) or (q x q)
		G = G + ((design_X.block(start, 0, m, q).transpose()) * ( ((R.inverse()) * r_i) * (r_i.transpose()*(R.inverse())) ) * design_X.block(start, 0, m, q));
	}

	sandwhich_Mat = tempGI.inverse() * G * tempGI.inverse();

	int after_EE = clock();


	// cout << " n start sum ended up being " << n_star_sum << endl;
	// cout << "m sum ended up being " << m_sum << endl;
	cout << "**************** GEE Results ****************" << endl;
	cout << endl << endl << "Iteration count at: " << iteration_count << endl;

	cout << "The correlation parameter: " << rho << endl;

	cout << endl << "Estimated Betas are: " << endl;
	cout << Betas_new.transpose() << endl << endl << endl;

	// cout << "The sandwhich matrix is: " << endl << endl;
	// cout << sandwhich_Mat << endl << endl;

	cout << "The beta Robust S.E.'s are: " << endl << endl;
	cout << sandwhich_Mat.diagonal().array().sqrt() << endl << endl;

	cout << " Estimated convergence time for EE is: " << endl;
	cout << " Time (secs): " << (after_EE - before_EE) / double(CLOCKS_PER_SEC) << endl;
	
	return 0;

}