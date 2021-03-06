/*
 * Particle.cpp
 *
 *  Created on: Oct 5, 2010
 *      Author: mgernot
 */

#include "Particle.h"

unsigned int Particle::counter = 0;

Particle::Particle(double p_x, double p_y, double r, double m, double f_x, double f_y, short int f) {

	number = counter++;

	position[0]  = p_x; position[1]  = p_y; phi  = 0.0;
	velocity[0]  = 0.0; velocity[1]  = 0.0; phi_dot  = 0.0;

	force_ext[0] = f_x; force_ext[1] = f_y; force_ext[2] = 0.0;
	force_int[0] = 0.0; force_int[1] = 0.0; force_int[2] = 0.0;
	force_old[0] = 0.0; force_old[1] = 0.0; force_old[2] = 0.0;

	radius = r;
	mass   = m;
	flag   = f;

}


unsigned int Contact::counter = 0;


Contact::Contact( Particle *P_i, Particle *P_j, const ublas::c_vector<double, 2> shiftI, const ublas::c_vector<double, 2> shiftJ){

	number = counter++;

	Particle_i = P_i;
	Particle_j = P_j;

	distance = ublas::norm_2(Particle_j->position - shiftJ - Particle_i->position + shiftI);
	gap = distance - Particle_i->radius - Particle_j->radius;

	n_ij(0) = ( Particle_i->position(0) - shiftI(0) - Particle_j->position(0) + shiftJ(0) ) / distance;
	n_ij(1) = ( Particle_i->position(1) - shiftI(1) - Particle_j->position(1) + shiftJ(1) ) / distance;

	t_ij(0) = ( Particle_i->position(1) - shiftI(1) - Particle_j->position(1) + shiftJ(1) ) / distance;
	t_ij(1) = ( Particle_j->position(0) - shiftJ(0) - Particle_i->position(0) + shiftI(0) ) / distance;

	t_survivial = 0;
	slip = false;
	in_touch = false;

}

bool Contact::updateContactPosition (Particle *P_i, Particle *P_j, const ublas::c_vector<double, 2> shiftI, const ublas::c_vector<double, 2> shiftJ){

	distance = ublas::norm_2(Particle_j->position - shiftJ - Particle_i->position + shiftI);
	gap = distance - Particle_i->radius - Particle_j->radius;

	n_ij(0) = ( Particle_i->position(0) - shiftI(0) - Particle_j->position(0) + shiftJ(0) ) / distance;
	n_ij(1) = ( Particle_i->position(1) - shiftI(1) - Particle_j->position(1) + shiftJ(1) ) / distance;

	t_ij(0) = ( Particle_i->position(1) - shiftI(1) - Particle_j->position(1) + shiftJ(1) ) / distance;
	t_ij(1) = ( Particle_j->position(0) - shiftJ(0) - Particle_i->position(0) + shiftI(0) ) / distance;

	if (gap > _SMALL_) {in_touch = false; t_survivial = 0;}
	else {in_touch = true; t_survivial++;}

	return(in_touch);
}

bool Contact::updateContactForce (double f_n, double f_t, double v_n, double v_t){

	f(0) = f_n;
	f(1) = f_t;

	v(0) = v_n;
	v(1) = v_t;

	if (v_t < _SMALL_) slip = false;
	else slip = true;

	return(slip);
}

void Contact::printContact () {

	std::cout << "Contact " << number << " between particle " << Particle_i->number
			<< " and particle " << Particle_j->number << std::endl
			<< "\t\tn_ij = (" << n_ij(0) << "," << n_ij(1) << ")"<< std::endl
			<< "\t\tt_ij = (" << t_ij(0) << "," << t_ij(1) << ")" << std::endl;

}

Contact::~Contact() {
	counter--;
}


Assembly::~Assembly() {

	for (unsigned int i= 0; i<P.size();i++)
		delete P.at(i);
}


int Assembly::addParticle(double p_x, double p_y, double r, double m,
		double f_x, double f_y, short int f) {

	P.push_back(new Particle(p_x, p_y, r, m, f_x, f_y, f));
	return(++P_count);
}

int Assembly::addContact (const int i, const int j, const ublas::c_vector<double, 2> shiftI, const ublas::c_vector<double, 2> shiftJ) {

	unsigned int key;
	std::pair<std::map <unsigned int, Contact *>::iterator, bool> ret;

	unsigned int N = this->P_count;

	if(i < j)
		key = N * i - ( (i - 1) * i ) / 2 + j - i;
	else
		key = N * j - ( (j - 1) * j ) / 2 + i - j;

	ret = C.insert (std::pair< unsigned int, Contact * >( key, new Contact(this->P[i], this->P[j], shiftI, shiftJ) ));
	if (ret.second == false) {
		// TODO: find a better solution for that or manage the counter from outside throughout the program...
		ret.first->second->counter--;
		ret.first->second->updateContactPosition(this->P[i], this->P[j], shiftI, shiftJ); }
	else
		C_count++;

return(C_count);

}

int Assembly::cleanUpContacts() {

	for(std::map<unsigned int, Contact *>::iterator iter = C.begin(); iter != C.end(); iter++) {
		if ((*iter).second->distance > ((*iter).second->Particle_i->radius + (*iter).second->Particle_j->radius) *  _CONTACT_APPROX_TRSH_)
			delete (*iter).second;
	}
}


void Assembly::generateMassMatrix() {

	int size_ = int (P.size());

	M.resize(3 * size_, 3 * size_, 0, 0);
	invM.resize(3 * size_, 3 * size_, 0, 0);

	for (int i = 0; i < size_; i++) {

		switch (P[i]->flag) {
		case 0: {
			M(i * 3, i * 3) = P[i]->mass;
			M(i * 3 + 1, i * 3 + 1) = P[i]->mass;
			M(i * 3 + 2, i * 3 + 2) = 1.0 * P[i]->mass * P[i]->radius * P[i]->radius;
			invM(i * 3, i * 3) = 1 / P[i]->mass;
			invM(i * 3 + 1, i * 3 + 1) = 1 / P[i]->mass;
			invM(i * 3 + 2, i * 3 + 2) =  1 /  P[i]->mass / P[i]->radius / P[i]->radius;

			P_free.push_back(P[i]);

			break;
		}
		case 1: {
			M(i * 3, i * 3) = _HUGE_;
			M(i * 3 + 1, i * 3 + 1) = _HUGE_;
			M(i * 3 + 2, i * 3 + 2) = _HUGE_;
			invM(i * 3, i * 3) = 0.;
			invM(i * 3 + 1, i * 3 + 1) = 0.;
			invM(i * 3 + 2, i * 3 + 2) = 0.;

			P_v_ctrl.push_back(P[i]);

			break;
		}
		case 2: {
			M(i * 3, i * 3) = P[i]->mass;
			M(i * 3 + 1, i * 3 + 1) = P[i]->mass;
			M(i * 3 + 2, i * 3 + 2) =  P[i]->mass * P[i]->radius * P[i]->radius;
			invM(i * 3, i * 3) = 1 / P[i]->mass;
			invM(i * 3 + 1, i * 3 + 1) = 1 / P[i]->mass;
			invM(i * 3 + 2, i * 3 + 2) = 1 / P[i]->mass / P[i]->radius / P[i]->radius;

			P_f_ctrl.push_back(P[i]);

			break;
		}
		}
	}
}

void Assembly::setInternalForces(ublas::vector<double> f_i){

	for (unsigned int i = 0; i < P_count; i++){
		P[i]->force_int(0)  = f_i(3*i);
		P[i]->force_int(1)  = f_i(3*i+1);
		P[i]->force_int(2)  = f_i(3*i+2);
	}
}


ublas::vector<double> Assembly::getParticleVelocities() {

	ublas::vector<double> v_ (3 * P_count);

	for (unsigned int i = 0; i < P_count; i++){
		v_(3*i) = 	P[i]->velocity(0);
		v_(3*i+1) = P[i]->velocity(1);
		v_(3*i+2) = P[i]->phi_dot;
	}

	return(v_);
}

ublas::vector<double> Assembly::getParticleForceExt() {

	ublas::vector<double> f_ (3 * P_count);

	for (unsigned int i = 0; i < P_count; i++){
		f_(3*i) = 	P[i]->force_ext(0);
		f_(3*i+1) = P[i]->force_ext(1);
		f_(3*i+2) = P[i]->force_ext(2);
	}

	return(f_);

}
