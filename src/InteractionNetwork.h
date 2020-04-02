/*
 * InteractionNetwork.h
 *
 *  Created on: Sep 8, 2010
 *      Author: mgernot
 */

#ifndef INTERACTIONNETWORK_H_
#define INTERACTIONNETWORK_H_

#include "Constants.h"
#include "Particle.h"
#include "Exceptions.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>
#include <iostream>

using namespace boost::numeric;

class InteractionNetwork{

	Assembly *F;

	ublas::vector<double> frc;
	ublas::vector<double> vel_free;
	ublas::vector<double> vel;

	std::vector<int> random_i; // contains a randomized list of all contact indices

	ublas::matrix<double> H; // fabric tensor, (2N x C+NC)
	ublas::vector<double> dist; // particle distance
	ublas::matrix<double> M_free;
	ublas::matrix<double> b; // b = H^T * M^(-1), mass matrix for conversion of external forces into accelerations (vector b), (C+NC x 2N)
	ublas::matrix<double> A; // A = H^T * M^(-1) * H, mass matrix for conversion of particle interaction forces into accelerations, (C+NC x C+NC)

	double epsilon; // precision of the iterative solver
	double mu; // coefficient of friction
	ublas::vector<double> x_periodicity;

	ublas::vector<double> p_i;
	ublas::vector<double> p_j;

public:
	InteractionNetwork(Assembly *fab, const ublas::vector<double> x_prd)
	: F(fab), epsilon(0.000001), mu(0.2), x_periodicity(x_prd) {   }

	InteractionNetwork(const InteractionNetwork &src) {
		// TODO: implement copy-constructor
		}

	~InteractionNetwork() {   }

	void updateNetwork();
	int setNetwork();
	void setFriction(const double m) {	mu = m; }
	void setPrecision(const double e) { epsilon = e; }

	void getForceWFriction(double dt) throw (slowConvergenceException);
	void GetMomentum(double dt) throw (slowConvergenceException);
};

#endif /* INTERACTIONNETWORK_H_ */
