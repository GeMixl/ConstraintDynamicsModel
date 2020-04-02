/*
 * IOAdapter.h
 *
 *  Created on: Sep 8, 2010
 *      Author: mgernot
 */

#ifndef IOADAPTER_H_
#define IOADAPTER_H_

#include "Particle.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <string>
#include <iterator>
#include <sstream>

using namespace boost::numeric;
using boost::lexical_cast;
using boost::bad_lexical_cast;

class InAdapter{
public:
	InAdapter(std::string filename, bool dummy);
	InAdapter(const InAdapter &src);
	~InAdapter();
	int GetNumberOfTimeSteps();
	double GetTimeStepLength();
	double GetContactFrictionCoefficient();
	double GetFlagTwo_vx();
	double GetFlagTwo_vy();
	double GetFlagThree_Fx();
	double GetFlagThree_Fy();
	double GetFlagThree_Fr();
	double GetFlagTwo_kx();
	double GetFlagTwo_ky();
	double GetDomainOriginX();
	double GetDomainOriginY();
	double GetDomainWidth();
	double GetDomainHeight();
	int    GetBoundaryConditionX();
	int    GetBoundaryConditionY();
	double GetIntegrationTolerance();
	int PushConfigToArray(Assembly &A);

private:
	int getFieldByKey(const char *k, std::string &v);

	std::fstream inData, File;
	std::vector<std::string *> Content;
	int N;
};

class OutAdapter {
public:
	OutAdapter(const std::string n, const Assembly *a);
	~OutAdapter();
	OutAdapter(const OutAdapter &src);

	void getParticleData(double t);
	void getContactData(double t);

	void printData();

private:
	std::string NameBase;
	std::vector< std::vector<double> *> ParticleData, ContactData;
	std::fstream outParticleData ,outContactData;
	const Assembly *ProcData;
};

#endif /* IOADAPTER_H_ */
