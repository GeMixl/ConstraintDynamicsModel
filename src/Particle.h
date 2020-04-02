/*
 * Particle.h
 *
 *  Created on: Oct 3, 2010
 *      Author: mixlmay
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "Constants.h"
#include "Exceptions.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <vector>
#include <map>

using namespace boost::numeric;


class Particle {

	static unsigned int counter;

public:
	ublas::c_vector<double, 2> position;
	double phi;
	ublas::c_vector<double, 2> velocity;
	double phi_dot;
	ublas::c_vector<double, 3> force_ext;
	ublas::c_vector<double, 3> force_int;
	ublas::c_vector<double, 3> force_old;

	double radius;
	double mass;

	double spring_constant;

	short int flag;

	unsigned int number;

	Particle(double p_x, double p_y, double r, double m, double f_x, double f_y, short int f);
	~Particle() {   }
	Particle(const Particle &src) {   }
};


class Contact {

public:
	Particle *Particle_i;
	Particle *Particle_j;

	ublas::c_vector<double, 2> n_ij;
	ublas::c_vector<double, 2> t_ij;

	ublas::c_vector<double, 2> f;
	ublas::c_vector<double, 2> v;

	ublas::c_matrix<double, 2, 2> DelassiusMatrix;

	double distance;
	double gap;

	static unsigned int counter;
	unsigned int number;
	unsigned int t_survivial;
	bool slip;
	bool in_touch;

	Contact(Particle *P_i, Particle *P_j, const ublas::c_vector<double, 2> shiftI, const ublas::c_vector<double, 2> shiftJ);
	~Contact();
	Contact(const Contact &src) {
	}

	bool updateContactPosition (Particle *P_i, Particle *P_j, const ublas::c_vector<double, 2> shiftI, const ublas::c_vector<double, 2> shiftJ);
	bool updateContactForce (double f_n, double f_t, double v_n, double v_t);
	void printContact ();

private:
	void setContactTime();
};


class Assembly {

	public:
	std::vector<Particle *> P;
	std::vector<Particle *> P_free;
	std::vector<Particle *> P_f_ctrl;
	std::vector<Particle *> P_v_ctrl;

//	std::vector<Contact *> C;
	std::map<unsigned int, Contact *> C;

	ublas::banded_matrix<double> M;
	ublas::banded_matrix<double> invM;

	unsigned int P_count, C_count;

	int addParticle(double p_x, double p_y, double r, double m, double f_x, double f_y, short int f);
	int addContact (const int i, const int j, const ublas::c_vector<double, 2> shiftI, const ublas::c_vector<double, 2> shiftJ);
	int cleanUpContacts();
	void generateMassMatrix();

	void setInternalForces(ublas::vector<double> f_i);
	void updateContacts();

	ublas::vector<double> getParticleVelocities();
	ublas::vector<double> getParticleForceExt();

	Assembly() : P_count(0), C_count(0) {   }
	~Assembly();
	Assembly(const Assembly &src) {   }
};


#endif /* PARTICLE_H_ */
