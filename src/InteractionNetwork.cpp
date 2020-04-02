/*
 * InteractionNetwork.cpp
 *
 *  Created on: Sep 8, 2010
 *      Author: mgernot
 */

#include "InteractionNetwork.h"


inline double signum(double val) {
	if (val < 0.) return (-1.);
	else return (1.);
}


void InteractionNetwork::updateNetwork() {

	for (unsigned int r = 0; r < F->P_count; r++) {
	// flip particles to the left side of the domain, if they move out at the right side
			if( F->P[r]->position(0) - F->P[r]->radius > x_periodicity(0) )
				F->P[r]->position(0) = F->P[r]->position(0) - x_periodicity(0);

	// flip particles to the right side of the domain, if they move out at the left side
			if( F->P[r]->position(0) + F->P[r]->radius < 0.)
				F->P[r]->position(0) = F->P[r]->position(0) + x_periodicity(0);
	}

}


int InteractionNetwork::setNetwork() {

	// TODO:
	//      *) exceptions

	// TODO: replace by a global threshold which is eventually compared against
	//       flight distances later...

	// counter for the dist vector and later for the matrix H
	unsigned int c = 0;
	// threshold for considering contacts
	double threshold_ = 0.;

	// clean up the contact list
	F->cleanUpContacts();

	// iteration over all potential contacts
	for (unsigned int i = 0; i < F->P_count - 1; i++) {
		for (unsigned int j = i + 1;	j < F->P_count; j++) {

			// do not consider contacts between two controlled particles
			if ( F->P[i]->flag > 0 && F->P[j]->flag > 0)
				continue;

			// calculate squares of maximum particle distances
			threshold_ = (F->P[i]->radius + F->P[j]->radius) *
						(F->P[i]->radius + F->P[j]->radius) *
						_CONTACT_APPROX_TRSH_ * _CONTACT_APPROX_TRSH_;

			// compare the particle distance square to 'threshold_'
			if ((   (F->P[j]->position[0] - F->P[i]->position[0]) *
					(F->P[j]->position[0] - F->P[i]->position[0]) +
				    (F->P[j]->position[1] - F->P[i]->position[1]) *
				    (F->P[j]->position[1] - F->P[i]->position[1]))
				    < threshold_)
				F->addContact(i, j, ublas::zero_vector<double>(2), ublas::zero_vector<double>(2));

			// compare the backward periodic particle distance square to 'threshold_'
			if ((   (F->P[j]->position[0] - F->P[i]->position[0] + x_periodicity[0]) *
				    (F->P[j]->position[0] - F->P[i]->position[0] + x_periodicity[0]) +
				    (F->P[j]->position[1] - F->P[i]->position[1]) *
				    (F->P[j]->position[1] - F->P[i]->position[1]))
				    < threshold_)
				F->addContact(i, j, x_periodicity, ublas::zero_vector<double>(2));

			// compare the forward periodic particle distance square to 'threshold_'
			if ((   (F->P[j]->position[0] - x_periodicity[0] - F->P[i]->position[0]) *
					(F->P[j]->position[0] - x_periodicity[0] - F->P[i]->position[0]) +
					(F->P[j]->position[1] - F->P[i]->position[1]) *
					(F->P[j]->position[1] - F->P[i]->position[1]))
					< threshold_)
				F->addContact(i, j, ublas::zero_vector<double>(2), x_periodicity);

		}
	}

	// clean up 1) the random vector, 2) the 'dist' vector and 3) the fabric tensor 'H'
	random_i.clear();
	dist.resize( F->C_count );
	H.resize(3 * F->P_count, 2 * F->C_count );

#ifdef VERBOSE
		std::cout << std::endl << "*) Contacts" << std::endl;
#endif

	for (std::map<unsigned int, Contact *>::iterator iter  = F->C.begin(); iter != F->C.end(); iter++) {

#ifdef VERBOSE
		(*iter).second->printContact();
#endif
		// initialize the random_i vector with integers [1:contact_count]
		random_i.push_back((*iter).second->number);

		if ((*iter).second->gap > 0.)
			dist(c) = (*iter).second->gap;
		else
			dist(c) = 0.;

		for (unsigned int j = 0; j < F->P_count; j++){

			// line to represent 'i' to 'j' interaction
			if (j == (*iter).second->Particle_i->number) {
			H(3*j, 2*c) = (*iter).second->n_ij(0);
			H(3*j+1, 2*c) = (*iter).second->n_ij(1);
			H(3*j+2, 2*c) = 0.;

			H(3*j, 2*c+1) = (*iter).second->t_ij(0);
			H(3*j+1, 2*c+1) = (*iter).second->t_ij(1);
			H(3*j+2, 2*c+1) = 1/(*iter).second->Particle_i->radius;
			}

			// line to represent 'j' to 'i' interaction
			else if (j == (*iter).second->Particle_j->number) {
			H(3*j, 2*c) = (-1) * (*iter).second->n_ij(0);
			H(3*j+1, 2*c) = (-1) * (*iter).second->n_ij(1);
			H(3*j+2, 2*c) = 0.;

			H(3*j, 2*c+1) = (-1) * (*iter).second->t_ij(0);
			H(3*j+1, 2*c+1) = (-1) * (*iter).second->t_ij(1);
			H(3*j+2, 2*c+1) = 1/(*iter).second->Particle_j->radius;
			}

			// line to represent 'j' to 'i' interaction
			else {
			H(3*j, 2*c) = 0.;
			H(3*j+1, 2*c) = 0.;
			H(3*j+2, 2*c) = 0.;

			H(3*j, 2*c+1) = 0.;
			H(3*j+1, 2*c+1) = 0.;
			H(3*j+2, 2*c+1) = 0.;
			}


		}
		c++;
	}

#ifdef VERBOSE
	std::cout << std::endl << "*) Fabric tensor" << std::endl;
	std::cout.precision(3);
	for (ublas::matrix<double>::iterator1 p = H.begin1(); p < H.end1(); p++) {
		for (ublas::matrix<double>::iterator2 q = p.begin(); q < p.end(); q++)
			std::cout << std::fixed << *q << "\t";
		std::cout << std::endl;

	}

	std::cout << std::endl << "*) Gaps vector" << std::endl;
	for (ublas::vector<double>::iterator p = dist.begin(); p < dist.end(); p++) {
		std::cout << std::fixed << *p << std::endl;
	}
	std::cout.precision(6);
#endif

	// trim all vectors that will contain properties of the interaction network
	// to the right size
	frc.resize(2 * F->C_count);
	vel_free.resize(2 * F->C_count);
	vel.resize(2 * F->C_count);

	// compute the fabric matrix A and the vector b for the iterative solver
	b = ublas::prod(ublas::trans(H), F->invM);
	A = ublas::prod(b, H);

#ifdef VERBOSE
	std::cout << std::endl << "*) Mass matrix A" << std::endl;
	std::cout.precision(3);
	for (ublas::matrix<double>::iterator1 p = A.begin1(); p < A.end1(); p++) {
		for (ublas::matrix<double>::iterator2 q = p.begin(); q < p.end(); q++)
			std::cout << std::fixed << *q << "\t";
		std::cout << std::endl;
	}
#endif

return ( F->C_count );
}


void InteractionNetwork::getForceWFriction(double dt)
	throw (slowConvergenceException) {

	double gap = 0.; // prediction for the gap between two particles at time t+dt

	double restitution = 0.;

	double frc_r = 0.; // normal interaction force at the current contact
	double frc_t = 0.; //tangential interaction force at the current contact
	double vel_new = 0.; // normal contact velocity at the current contact
	double w_new = 0.; // tangential contact velocity at the current contact

	double norm_old = _HUGE_;
	double norm_new = 0.;

	double lambda = 1.0;

	int i, j;

	const unsigned int contact_count = F->C_count; // number of potential contacts at a current configuration
	const unsigned int particle_count = F->P_count; // number of particles


	// transform particle velocities into contact velocities using the fabric tensor H
	vel = ublas::prod(ublas::trans(H), F->getParticleVelocities());

	// calculate the acceleration due to external forces b
	vel_free = vel + (ublas::prod(b, F->getParticleForceExt()) * dt);

	// set the contact forces back to zero
	frc.clear();

	j = 0;
	while (j++ < (int) (contact_count * contact_count) && std::fabs(
			(norm_new - norm_old)/norm_new) > epsilon) {

		// TODO: implement the slowConvergencException
		// check the maximum number of iterations
//		if (j > ???) throw slowConvergenceException();

		// randomize the index vector
		std::random_shuffle(random_i.begin(), random_i.end());

		norm_old = norm_new;
		// iterate over all contacts and find the ones with g < 0.
		for (std::vector<int>::iterator iter = random_i.begin(); iter < random_i.end(); iter ++) {

			// assign the randomized index from random_i to i
			i = *(iter);

			// get the i-th column of A
			ublas::matrix_column<ublas::matrix<double> > A_col_n(A, 2 * i);
			ublas::matrix_column<ublas::matrix<double> > A_col_t(A, 2 * i + 1);

			// calculate the new gap between particles under free (only external forces) acceleration
			gap = dist(i) + vel_free(2 * i) * dt;

			// particles do not touch, no interaction
			if (gap > 0.) {
				// normal interaction
				vel_new = vel_free(2 * i);
				frc_r = 0.;
				// tangential interaction
				w_new = vel_free(2 * i + 1);
				frc_t = 0.;
			}

			// particles touch, interaction takes place
			else {
				// normal interaction
				vel_new = -1*(restitution) * dist(i) / dt;
				frc_r = -1 / dt * lambda / A_col_n(2 * i) * (vel_free(2 * i) - vel_new);
				// tangential interaction
				w_new = 0.;
				// tangential force is divided into a rotational part and a
				// translational part by division through the tangential interaction
				// entry in A
				frc_t = -1 / dt * lambda / A_col_t(2 * i + 1) * vel_free(2 * i + 1);
				// check if the Coulomb criterion is fulfilled
				if (std::fabs(frc_t) > frc_r * mu) {
					frc_t = mu * frc_r * signum(frc_t);
					// the minus sign cancels out since w_new points into opposite direction of frc_t
					w_new = dt * A_col_t(2 * i + 1) * frc_t + vel_free(2*i+1);
				}
			}

			// calculate the new vel_free due to external forces, internal forces and particle velocities
			frc(2 * i) = frc(2 * i) + frc_r;
			frc(2 * i + 1) = frc(2 * i + 1) + frc_t;
			vel(2 * i) = vel_new;
			vel(2 * i + 1) = w_new;
			vel_free = vel_free + (A_col_n * frc_r * dt) + (A_col_t * frc_t * dt);
		}

		norm_new = ublas::norm_2(frc);
	}

#ifdef VERBOSE
	std::cout.precision(4);

	std::cout << std::endl << "*) Contact forces" << std::endl;
	for (ublas::vector<double>::iterator p = frc.begin(); p < frc.end(); p++) {
		std::cout << std::fixed << *p << std::endl;
	}
	std::cout << std::endl << "*) Contact velocities" << std::endl;
	for (ublas::vector<double>::iterator p = vel.begin(); p < vel.end(); p++) {
		std::cout << std::fixed << *p << std::endl;
	}
//	std::cout << std::endl << "*) Contact status" << std::endl;
//	for (std::map< unsigned int, Contact * >::iterator p = F->C.begin(); p != F->C.end(); p++) {
//		std::cout << "Contact " << (*p).first << ": " << std::fixed << (*p).second->slip << std::endl;
//	}
	std::cout << std::endl << "*) Number of iterations" << std::endl << j << std::endl;
	std::cout.precision(6);
#endif

	F->setInternalForces(ublas::prod(H, frc));

	j = 0;
// update the contacts with new forces and velocities
	for (std::map<unsigned int, Contact*>::iterator iter  = F->C.begin(); iter != F->C.end(); iter++) {
		iter->second->updateContactForce(frc(2*j), frc(2*j+1), vel(2*j), vel(2*j+1));
		j++;
	}
}
