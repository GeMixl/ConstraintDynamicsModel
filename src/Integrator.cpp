/*
 * Integrator.cpp
 *
 *  Created on: Oct 13, 2010
 *      Author: mgernot
 */

#include "Integrator.h"

void Integrator::setControlVelocity(const double &vx, const double &vy) {
	velocity_ctrl(0) = vx; velocity_ctrl(1) = vy;
}


void Integrator::setControlForce(const double &fx, const double &fy, const double &fr){
	force_ctrl(0) = fx; force_ctrl(1) = fy; force_ctrl(2) = fr;
}

void Integrator::setControlDomain(const double &x, const double &y, const double &w, const double &h) {
	domain(0) = x; domain(1) = y; domain(2) = w; domain(3) = h;
}

void Integrator::setControlledParticleStiffness(const double &s) {
	stiffness = s;
}


void Integrator::YlevelControlledParticles(std::vector<Particle *>::iterator b, std::vector<Particle *>::iterator e) {
	double _avrgYpos = 0.;

	// get the average y-position of all force controlled particles
	for (std::vector<Particle *>::iterator iter = b; iter < e; iter++) {
		_avrgYpos = _avrgYpos + (*iter)->position(1);
	}
	_avrgYpos = _avrgYpos / (double)F->P_f_ctrl.size();

	// set additional forces on the force controlled particles with respect to their relative position
	// f_add = k * delta pos_y     ...with     delta_pos_y = pos_y - averg_pos_y
	for (std::vector<Particle *>::iterator iter = b; iter < e; iter++) {
		(*iter)->force_ext(1) = force_ctrl(1) + stiffness * (_avrgYpos - (*iter)->position(1));
		(*iter)->force_ext(0) = (-1)*(*iter)->force_int(0);
	}

}


double Integrator::doIntegration(double dt) {

	double tmp = 0.;

	YlevelControlledParticles(F->P_f_ctrl.begin(), F->P_f_ctrl.end());
	//YlevelControlledParticles(F->P_v_ctrl.begin(), F->P_v_ctrl.end());

	for (std::vector<Particle *>::iterator iter = F->P_f_ctrl.begin(); iter < F->P_f_ctrl.end(); iter++){
	//	(*iter)->force_ext = (*iter)->force_ext + force_ctrl;
	}

	for (std::vector<Particle *>::iterator iter = F->P_v_ctrl.begin(); iter < F->P_v_ctrl.end(); iter++){
		(*iter)->velocity = velocity_ctrl;

		// ##############################
		// tmp = tmp + (*iter)->force_int(0);
		// ##############################

		(*iter)->force_int = -1*(*iter)->force_ext;
	}
	// ##############################
	// std::cout << tmp << "\t" << F->P[27]->phi *180 / 3.1415 << "\t" << F->P[28]->phi *180 / 3.1415 << std::endl;
	// ##############################

	// perform integration with all particles
	for (unsigned int i = 0; i < F->P_count; i++) {
		F->P[i]->velocity(0) = F->P[i]->velocity(0) + dt * F->invM(3*i, 3*i) * (F->P[i]->force_ext(0) + F->P[i]->force_int(0)) ;
		F->P[i]->velocity(1) = F->P[i]->velocity(1) + dt * F->invM(3*i+1, 3*i+1) * (F->P[i]->force_ext(1) + F->P[i]->force_int(1)) ;
		F->P[i]->phi_dot = F->P[i]->phi_dot + dt * F->invM(3*i+2, 3*i+2) * (F->P[i]->force_ext(2) + F->P[i]->force_int(2)) ;

//		F->P[i]->phi_dot = dt * F->invM(3*i+2, 3*i+2) * (F->P[i]->force_ext(2) + F->P[i]->force_int(2)) ;


//		if (std::fabs(F->P[i]->phi_dot + dt * F->invM(3*i+2, 3*i+2) * (F->P[i]->force_ext(2) + F->P[i]->force_int(2))) > 0.15)
//			F->P[i]->phi_dot = F->P[i]->phi_dot + dt * F->invM(3*i+2, 3*i+2) * (F->P[i]->force_ext(2) + F->P[i]->force_int(2)) ;
//		else
//			F->P[i]->phi_dot = 0.;

		F->P[i]->position(0) = F->P[i]->position(0) + dt * F->P[i]->velocity(0);
		F->P[i]->position(1) = F->P[i]->position(1) + dt * F->P[i]->velocity(1);
		F->P[i]->phi = F->P[i]->phi + dt * F->P[i]->phi_dot;
	}
	//std::cout << "# " << F_H << "\t" << F_V << "\t" << F_H/F_V << std::endl;


	return (dt);
}
