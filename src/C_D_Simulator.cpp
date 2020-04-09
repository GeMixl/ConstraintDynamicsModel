/*
 * C_D_Simulator.cpp
 *
 *  Created on: Sep 1, 2010
 *      Author: mgernot
 */

// TODO:
/* exceptions:
 *   -) periodicity is smaller than 3-times the largest particle diameter
 */

//#include "ContactNetwork.h"
#include "IOAdapter.h"
#include "InteractionNetwork.h"
#include "Exceptions.h"
#include "Integrator.h"
#include "Particle.h"
#include "Visualization.h"

#include <fstream>
#include <string>

Fl_Text_Buffer *textbuf;
Fl_Text_Buffer *stylebuf;

char *checkremove(int *argc, char ***argv, char *opt)
{
  int i=0;
  int n=strlen(opt);
  char* ret=0;
  while (i<*argc)
    {
      if (strncmp((*argv)[i],opt,n)==0)
	{
	  ret=(*argv)[i]+n;
	  for (int j=i;j<*argc-1;j++) (*argv)[j]=(*argv)[j+1];
	  (*argc)--;
	} else i++;
    }
  return ret;
}


int main(int argc, char **argv) {

	int flag = 0;
	bool i_flag, o_flag;
	char *infile_name, *outfile_name;

	infile_name = checkremove(&argc,&argv,"--infile=");
	i_flag = (infile_name!=NULL);
	if (!i_flag) {
		infile_name = "setup.ini";
	}

	outfile_name = checkremove(&argc,&argv,"--outfile=");
	o_flag = (outfile_name!=NULL);
	if (!o_flag) {
		outfile_name = "out";
	}

	Assembly TheAssembly;
	InAdapter TheInAdapter(std::string(infile_name), 0);
	OutAdapter TheOutAdapter(std::string(outfile_name), &TheAssembly);

	double dt = TheInAdapter.GetTimeStepLength(); // integration time step
	int step = 0; //initial step
	int step_max = TheInAdapter.GetNumberOfTimeSteps();
	TheInAdapter.PushConfigToArray(TheAssembly);

	TheAssembly.generateMassMatrix();

	// TODO: find a better solution for that...
	ublas::vector<double> x_prd_tmp(2);
	x_prd_tmp(0) = 16.001;
	x_prd_tmp(1) = 0.;

	InteractionNetwork TheNetwork(&TheAssembly, x_prd_tmp);
	TheNetwork.setFriction(TheInAdapter.GetContactFrictionCoefficient());
	TheNetwork.setPrecision(TheInAdapter.GetIntegrationTolerance());

	Integrator TheIntegrator(&TheAssembly);
	TheIntegrator.setControlVelocity(
			TheInAdapter.GetFlagTwo_vx(),
			TheInAdapter.GetFlagTwo_vy() );
	TheIntegrator.setControlForce(
			TheInAdapter.GetFlagThree_Fx(),
			TheInAdapter.GetFlagThree_Fy(),
			TheInAdapter.GetFlagThree_Fr() );
	TheIntegrator.setControlDomain(
			TheInAdapter.GetDomainOriginX(),
			TheInAdapter.GetDomainOriginY(),
			TheInAdapter.GetDomainWidth(),
			TheInAdapter.GetDomainHeight() );
	TheIntegrator.setControlledParticleStiffness(TheInAdapter.GetFlagTwo_kx());

#ifdef VERBOSE
	std::cout << std::endl << "*) Particles" << std::endl;
	std::cout << "x pos.\ty pos.\tradius\tmass\tF_x\tF_y\tflag" << std::endl;
	for (unsigned int q = 0; q < TheAssembly.P_count; q++)
		std::cout << TheAssembly.P[q]->position[0] << "\t"
		          << TheAssembly.P[q]->position[1] << "\t"
		          << TheAssembly.P[q]->radius << "\t"
   		          << TheAssembly.P[q]->mass << "\t"
   		          << TheAssembly.P[q]->force_ext[0] << "\t"
		          << TheAssembly.P[q]->force_ext[1] << "\t"
		          << TheAssembly.P[q]->flag << "\t" << std::endl;

	std::cout << std::endl << "*) Mass matrix" << std::endl;
	for (ublas::banded_matrix<double>::iterator1 p = TheAssembly.M.begin1();
			p < TheAssembly.M.end1(); p++) {
		for (ublas::banded_matrix<double>::iterator2 q = p.begin();
				q < p.end(); q++)
			std::cout << "M(" << q.index1() << ";" << q.index2() << ") = "<< *q << "\t";
		std::cout << std::endl;
	}

	std::cout << std::endl << "*) Inverted mass matrix" << std::endl;
	for (ublas::banded_matrix<double>::iterator1 p = TheAssembly.invM.begin1();
			p < TheAssembly.invM.end1(); p++) {
		for (ublas::banded_matrix<double>::iterator2 q = p.begin();
				q < p.end(); q++)
			std::cout << "invM(" << q.index1() << ";" << q.index2() << ") = " << *q << "\t";
		std::cout << std::endl;
	}
	std::cout << std::endl;
#endif

	TheNetwork.setNetwork();

	// Moreau algorithm, by Unger&Kertesz using the iterative solver for force
	// redistribution.
	// =======================================================================

	// TODO: implement the periodicity zero case in *InteractionNetwork*
	// for the time being DO NOT SET periodicity to zero!!!

	while (step++ < step_max) {
		try {

#ifdef VERBOSE
		std::cout << std::endl << "### START STEP " << step << " ###" << std::endl;
#endif

//		Fabric_II->UpdateForces();
		TheNetwork.updateNetwork();
		TheNetwork.setNetwork();
		TheNetwork.getForceWFriction(dt);

#ifndef VERBOSE
		std::cout << "Computing step "<< step << ", t = " <<  step*dt << "\n";
#endif

		TheIntegrator.doIntegration(dt);

#ifdef VERBOSE
	std::cout.precision(3);

	std::cout << std::endl << ">>>PARTICLE STATUS REPORT<<<" << std::endl;
	std::cout << std::endl << "*) Particle positions" << std::endl;
	for (std::vector<Particle *>::iterator p = TheAssembly.P.begin(); p < TheAssembly.P.end(); p++) {
		std::cout << std::fixed << (*p)->position << ", phi:" << (*p)->phi << std::endl;
	}
	std::cout << std::endl << "*) Particle velocities" << std::endl;
	for (std::vector<Particle *>::iterator p = TheAssembly.P.begin(); p < TheAssembly.P.end(); p++) {
		std::cout << std::fixed << (*p)->velocity  << ", omega:" << (*p)->phi_dot << std::endl;
	}
	std::cout << std::endl << "*) Particle forces" << std::endl;
	for (std::vector<Particle *>::iterator p = TheAssembly.P.begin(); p < TheAssembly.P.end(); p++) {
		std::cout << std::fixed << (*p)->force_int + (*p)->force_ext << std::endl;
	}

	std::cout << std::endl << ">>>CONTACT STATUS REPORT<<<" << std::endl;
	for (std::map<unsigned int, Contact*>::iterator iter  = TheAssembly.C.begin(); iter != TheAssembly.C.end(); iter++) {
		std::cout << "Contact " << (*iter).first << " slip status = " << (*iter).second->slip << std::endl;
		std::cout << "          " << "in touch = " << (*iter).second->in_touch << std::endl;
		std::cout << "          " << "touch since " << (*iter).second->t_survivial << " timesteps" << std::endl;
		}

	std::cout << std::endl << "### END STEP " << step << " ###" << std::endl;
	std::cout.precision(6);
#endif

	TheOutAdapter.getContactData((double)step*dt);
	TheOutAdapter.getParticleData((double)step*dt);

		}

		catch (slowConvergenceException){
			std::cout <<            "*** Exception!"
					<< std::endl << "*** Slow convergence of iteration"
					<< std::endl << "*** found at timestep " << step << "! " << std::endl;
			break;
		}


	}

	TheOutAdapter.printData();

	textbuf = new Fl_Text_Buffer(0);
	style_init();
	MyWindow win(900, 600, "C_D_Simulator", std::string(outfile_name));
    win.resizable(win);
    win.show();

	return(Fl::run());
}
