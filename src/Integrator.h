/*
 * Integrator.h
 *
 *  Created on: Oct 13, 2010
 *      Author: mgernot
 */

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include "Particle.h"

class Integrator {

protected:

	Assembly *F;

	ublas::c_vector<double, 3> force_ctrl;
	ublas::c_vector<double, 2> velocity_ctrl;
	ublas::c_vector<double, 4> domain;
	double stiffness;

public:
	Integrator(Assembly *A) : F(A), stiffness(0.) {   }
	~Integrator() {   };
	Integrator(const Integrator &src);

	double doIntegration(const double dt);
	void YlevelControlledParticles(std::vector<Particle *>::iterator b, std::vector<Particle *>::iterator e);
	void setControlForce(const double &fx, const double &fy, const double &fr);
	void setControlVelocity(const double &vx, const double &vy);
	void setControlDomain(const double &x, const double &y, const double &w, const double &h);
	void setControlledParticleStiffness(const double &s);


};

#endif /* INTEGRATOR_H_ */
