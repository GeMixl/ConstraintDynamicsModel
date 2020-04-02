/*
 * IOAdapter.cpp
 *
 *  Created on: Sep 8, 2010
 *      Author: mgernot
 */
#include "IOAdapter.h"

inline std::string Double2String(const double v){
	std::ostringstream output;
	if(!(output << v))
		throw BadConversion("Double2String could not convert value!");
	return(output.str());
}

inline std::string Float2String(const float v){
	std::ostringstream output;
	if(!(output << v))
		throw BadConversion("Float2String could not convert value!");
	return(output.str());
}

inline std::string Int2String(const int v){
	std::ostringstream output;
	if(!(output << v))
		throw BadConversion("Int2String could not convert value!");
	return(output.str());
}

inline std::string UInt2String(const unsigned int v){
	std::ostringstream output;
	if(!(output << v))
		throw BadConversion("UInt2String could not convert value!");
	return(output.str());
}

inline double DoubleFromString(const std::string &in, const int nr, const char *delim = NULL){
	double out;
	if (!delim)
		delim = "\f\n\r\t\v ";
	unsigned int begin = in.find_first_not_of(delim);
	unsigned int end   = in.find_first_of(delim, begin);

	for(int i=0;i<nr;i++) {
		if (begin==std::string::npos)
			throw InvalidDataReqest("DoubleFromString tries to extract " + Int2String(nr) + " from \"" + in + "\"");
		begin = in.find_first_not_of(delim, end);
		end   = in.find_first_of(delim, begin);
	}
	std::istringstream input(in.substr(begin, end));
	if (!(input >> out))
		throw BadConversion("DoubleFromString could not convert string \"" + in + "\"");
	return (out);
}

inline float FloatFromString(const std::string &in, const int nr, const char *delim = NULL){
	float out;
	if (!delim)
		delim = "\f\n\r\t\v ";
	unsigned int begin = in.find_first_not_of(delim);
	unsigned int end   = in.find_first_of(delim, begin);

	for(int i=0;i<nr;i++) {
		if (begin==std::string::npos)
			throw InvalidDataReqest("FloatFromString tries to extract " + Int2String(nr) + " from \"" + in + "\"");
		begin = in.find_first_not_of(delim, end);
		end   = in.find_first_of(delim, begin);
	}
	std::istringstream input(in.substr(begin, end));
	if (!(input >> out))
		throw BadConversion("DoubleFromString could not convert string \"" + input.str() + "\"");
	return (out);
}

inline int IntFromString(const std::string &in, const int nr, const char *delim = NULL){
	int out;
	if (!delim)
		delim = "\f\n\r\t\v ";
	unsigned int begin = in.find_first_not_of(delim);
	unsigned int end   = in.find_first_of(delim, begin);

	for(int i=0;i<nr;i++) {
		if (begin==std::string::npos)
			throw InvalidDataReqest("IntFromString (\"" + in + "\")");
		begin = in.find_first_not_of(delim, end);
		end   = in.find_first_of(delim, begin);
	}
	std::istringstream input(in.substr(begin, end));
	if (!(input >> out))
		throw BadConversion("IntFromString could not convert string \"" + input.str() + "\"");
	return (out);
}

inline unsigned int UIntFromString(const std::string &in, const int nr, const char *delim = NULL){
	unsigned int out;
	if (!delim)
		delim = "\f\n\r\t\v ";
	unsigned int begin = in.find_first_not_of(delim);
	unsigned int end   = in.find_first_of(delim, begin);

	for(int i=0;i<nr;i++) {
		if (begin==std::string::npos)
			throw InvalidDataReqest("UIntFromString (\"" + in + "\")");
		begin = in.find_first_not_of(delim, end);
		end   = in.find_first_of(delim, begin);
	}
	std::istringstream input(in.substr(begin, end));
	if (!(input >> out))
		throw BadConversion("UIntFromString could not convert string \"" + input.str() + "\"");
	return (out);
}

inline unsigned int SIntFromString(const std::string &in, const int nr, const char *delim = NULL){
	short int out;
	if (!delim)
		delim = "\f\n\r\t\v ";
	unsigned int begin = in.find_first_not_of(delim);
	unsigned int end   = in.find_first_of(delim, begin);

	for(int i=0;i<nr;i++) {
		if (begin==std::string::npos || end==std::string::npos)
			throw InvalidDataReqest("SIntFromString tries to extract " + Int2String(nr) + " from \"" + in + "\"");
		begin = in.find_first_not_of(delim, end);
		end   = in.find_first_of(delim, begin);
	}
	std::istringstream input(in.substr(begin, end));
	if (!(input >> out))
		throw BadConversion("SIntFromString could not convert string \"" + input.str() + "\"");
	return (out);
}

inline void toString(const int value, char* output)
{ sprintf(output, "%i", value) ; }

inline void toString(const unsigned int value, char* output)
{ sprintf(output, "%u", value) ; }


inline double toDouble(std::string const& s)
{
  std::istringstream i(s);
  double x;
  if (!(i >> x))
	  throw BadConversion("convertToDouble(\"" + s + "\")");
  return x;
}

inline std::string splitString(std::string const& s, const int nr, const char* delim = NULL) {
	if(!delim)
		delim = "\f\n\r\t\v ";
	unsigned int begin = s.find_first_not_of(delim);
	unsigned int end   = s.find_first_of(delim, begin);

	for(int i = 1;i<nr;i++) {
		if (begin==std::string::npos)
			throw InvalidDataReqest("splitString (\"" + s + "\")");
		begin = s.find_first_not_of(delim, end);
		end   = s.find_first_of(delim, begin);
	}
	return(s.substr(begin,end-begin));
}

bool checkForComment(std::string &s, const char* flag = NULL) {

	if(!flag)
		flag = "%#";
	std::size_t idx = s.find_first_of(flag);
	if (idx==0) {
		return(true); }
	else if (idx==std::string::npos) {
		return(false);}
	else{
		s=s.substr(0,idx);
		return (false);}
}

InAdapter::InAdapter(std::string filename, bool dummy) {

	std::string tmp;
	inData.open(filename.c_str(), std::fstream::in);

	while (!inData.eof()) {
		getline(inData, tmp, '\n');
		if (!checkForComment(tmp, "#"))
			Content.push_back(new std::string(tmp));
		tmp.clear();
	}

	inData.close();
}

InAdapter::~InAdapter() {


}

int InAdapter::getFieldByKey(const char *k, std::string &v) {

	/* USAGE:
	 * text file contains command lines starting with '*'
	 * e.g.: "*COMMAND 1"
	 * The function finds the line from InAdapter::Content by passing the key "*COMMAND" in k
	 * It returns the line number and writes an argument, if present, into v.
	 */
	int i = -1;
	const char *delim = "\f\n\r\t\v ";
	std::string::size_type beg, end;

	for(std::vector<std::string *>::iterator iter = Content.begin(); iter != Content.end(); iter++){
		if((*iter)->find("*")==0){
			end = (*iter)->find_first_of(delim);
			if((*iter)->substr(0,end) == k) {
				i = std::distance(Content.begin(), iter);
				std::cout << "found requested key " << k << " at position " << i << std::endl;
				beg = (*iter)->find_first_not_of(delim, end);
				if (beg!=std::string::npos){
					end = (*iter)->find_first_of(delim, beg);
					v=(*iter)->substr(beg, end);
				}
			}
		}
	}
	return(i);
}

int InAdapter::PushConfigToArray(Assembly &A) {

	unsigned int DATA_FILE_PX_POS = 0;
	unsigned int DATA_FILE_PY_POS = 1;
	unsigned int DATA_FILE_RA_POS = 2;
	unsigned int DATA_FILE_MA_POS = 3;
	unsigned int DATA_FILE_FX_POS = 4;
	unsigned int DATA_FILE_FY_POS = 5;
	unsigned int DATA_FILE_FG_POS = 6;

	std::string dummy;
	int sta = getFieldByKey("*PARTICLE_DATA_START", dummy)+1;
	int sto = getFieldByKey("*PARTICLE_DATA_STOP", dummy)-1;
	unsigned int N = sto-sta;

	for (int i=sta;i<=sto;i++){
		A.addParticle(
				DoubleFromString(*(Content[i]), DATA_FILE_PX_POS),
				DoubleFromString(*(Content[i]), DATA_FILE_PY_POS),
				DoubleFromString(*(Content[i]), DATA_FILE_RA_POS),
				DoubleFromString(*(Content[i]), DATA_FILE_MA_POS),
				DoubleFromString(*(Content[i]), DATA_FILE_FX_POS),
				DoubleFromString(*(Content[i]), DATA_FILE_FY_POS),
				SIntFromString(*(Content[i]), DATA_FILE_FG_POS));
		}
	return(N);
}

int InAdapter::GetNumberOfTimeSteps() {

	std::string buf;
	int val;
	if (getFieldByKey("*NUMBER_OF_TIME_STEPS", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetNumberOfTimeSteps() found value " << buf << std::endl;
	val = boost::lexical_cast<int>(buf);
	if (val<1 || val>100000)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}

double InAdapter::GetTimeStepLength() {

	std::string buf;
	double val;
	if (getFieldByKey("*TIME_STEP_LENGTH", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetTimeStepLength() found value " << buf << std::endl;
	val = boost::lexical_cast<double>(buf);
	if (val<0.000000000001 || val>100)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}

double InAdapter::GetContactFrictionCoefficient(){

	std::string buf;
	double val;
	if (getFieldByKey("*CONTACT_FRICTION_COEFFICIENT", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetContactFrictionCoefficient() found value " << buf << std::endl;
	val = boost::lexical_cast<double>(buf);
	if (val<0. || val>1.)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}

double InAdapter::GetFlagTwo_vx() {

	std::string buf;
	double val;
	if (getFieldByKey("*FLAG_TWO_PARTICLE_V_X", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetFlagTwo_vx() found value " << buf << std::endl;
	val = boost::lexical_cast<double>(buf);
	if (val<-1000. || val>1000.)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}

double InAdapter::GetFlagTwo_vy() {

	std::string buf;
	double val;
	if (getFieldByKey("*FLAG_TWO_PARTICLE_V_X", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetFlagTwo_vy( found value " << buf << std::endl;
	val = boost::lexical_cast<double>(buf);
	if (val<-1000. || val>1000.)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}

double InAdapter::GetFlagTwo_kx() {

	std::string buf;
	double val;
	if (getFieldByKey("*FLAG_TWO_PARTICLE_K_X", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetFlagTwo_kx() found value " << buf << std::endl;
	val = boost::lexical_cast<double>(buf);
	if (val<-1000. || val>1000.)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}

double InAdapter::GetFlagTwo_ky() {

	std::string buf;
	double val;
	if (getFieldByKey("*FLAG_TWO_PARTICLE_K_Y", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetFlagTwo_ky() found value " << buf << std::endl;
	val = boost::lexical_cast<double>(buf);
	if (val<-1000. || val>1000.)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}

double InAdapter::GetFlagThree_Fx() {

	std::string buf;
	double val;
	if (getFieldByKey("*FLAG_THREE_PARTICLE_F_X", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetFlagThree_Fx() found value " << buf << std::endl;
	val = boost::lexical_cast<double>(buf);
	if (val<-1000. || val>1000.)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}

double InAdapter::GetFlagThree_Fy() {

	std::string buf;
	double val;
	if (getFieldByKey("*FLAG_THREE_PARTICLE_F_Y", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetFlagThree_Fy() found value " << buf << std::endl;
	val = boost::lexical_cast<double>(buf);
	if (val<-1000. || val>1000.)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}

double InAdapter::GetFlagThree_Fr() {

	std::string buf;
	double val;
	if (getFieldByKey("*FLAG_THREE_PARTICLE_F_R", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetFlagThree_Fr() found value " << buf << std::endl;
	val = boost::lexical_cast<double>(buf);
	if (val<-1000. || val>1000.)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}


double InAdapter::GetDomainOriginX() {

	std::string buf;
	double val;
	if (getFieldByKey("*DOMAIN_ORIGIN_X", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetDomainOriginX() found value " << buf << std::endl;
	val = boost::lexical_cast<double>(buf);
	if (val<-1000. || val>1000.)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}

double InAdapter::GetDomainOriginY() {

	std::string buf;
	double val;
	if (getFieldByKey("*DOMAIN_ORIGIN_Y", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetDomainOriginY() found value " << buf << std::endl;
	val = boost::lexical_cast<double>(buf);
	if (val<-1000. || val>1000.)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}

double InAdapter::GetDomainWidth() {

	std::string buf;
	double val;
	if (getFieldByKey("*DOMAIN_WIDTH", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetDomainWidth() found value " << buf << std::endl;
	val = boost::lexical_cast<double>(buf);
	if (val<-1000. || val>1000.)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}

double InAdapter::GetDomainHeight() {

	std::string buf;
	double val;
	if (getFieldByKey("*DOMAIN_HEIGHT", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetDomainHeight() found value " << buf << std::endl;
	val = boost::lexical_cast<double>(buf);
	if (val<-1000. || val>1000.)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}

int InAdapter::GetBoundaryConditionX() {

	std::string buf;
	int val;
	if (getFieldByKey("*BOUNDARY_CONDITION_X", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetBounaryConditionX() found value " << buf << std::endl;
	val = boost::lexical_cast<int>(buf);
	if (val<-1000 || val>1000)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}

int InAdapter::GetBoundaryConditionY() {

	std::string buf;
	int val;
	if (getFieldByKey("*BOUNDARY_CONDITION_Y", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetBounaryConditionX() found value " << buf << std::endl;
	val = boost::lexical_cast<int>(buf);
	if (val<-1000 || val>1000)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}

double InAdapter::GetIntegrationTolerance() {

	std::string buf;
	double val;
	if (getFieldByKey("*EPSILON_FORCE", buf)<0)
		;//throw an exception --> the string was not defined in the control file
	std::cout << "InAdapter::GetContactFrictionCoefficient() found value " << buf << std::endl;
	val = boost::lexical_cast<double>(buf);
	if (val<0. || val>1000.)
		;//throw an exception --> the number of time steps is too large/small
	return(val);
}

OutAdapter::OutAdapter(std::string n, const Assembly *a) {

	NameBase = n;
	ProcData = a;
}

OutAdapter::~OutAdapter() {

}

OutAdapter::OutAdapter(const OutAdapter &src) {

}

void OutAdapter::getParticleData(double t) {

	ParticleData.push_back(new std::vector<double>);

	for (unsigned int q = 0; q < ProcData->P_count; q++) {
		if (q==0) ParticleData.back()->push_back(t);
		ParticleData.back()->push_back(ProcData->P[q]->position[0]);
		ParticleData.back()->push_back(ProcData->P[q]->position[1]);
		ParticleData.back()->push_back(ProcData->P[q]->phi);
		ParticleData.back()->push_back(ProcData->P[q]->radius);
		ParticleData.back()->push_back(ProcData->P[q]->mass);
		ParticleData.back()->push_back(ProcData->P[q]->velocity[0]);
		ParticleData.back()->push_back(ProcData->P[q]->velocity[1]);
		ParticleData.back()->push_back(ProcData->P[q]->phi_dot); }
}

void OutAdapter::getContactData(double t) {

	ContactData.push_back(new std::vector<double>);
	for (std::map<unsigned int, Contact*>::const_iterator iter  = ProcData->C.begin(); iter != ProcData->C.end(); iter++) {
		if (iter==ProcData->C.begin()) ContactData.back()->push_back(t);
		ContactData.back()->push_back(iter->second->Particle_i->number);
		ContactData.back()->push_back(iter->second->Particle_j->number);
		ContactData.back()->push_back(iter->second->n_ij[0]);
		ContactData.back()->push_back(iter->second->n_ij[1]);
		ContactData.back()->push_back(iter->second->t_ij[0]);
		ContactData.back()->push_back(iter->second->t_ij[1]);
		ContactData.back()->push_back(iter->second->v[0]);
		ContactData.back()->push_back(iter->second->v[1]);
		ContactData.back()->push_back(iter->second->f[0]);
		ContactData.back()->push_back(iter->second->f[1]); }
}


void OutAdapter::printData() {

	std::string ptcl_filename = NameBase + std::string("_particles.dat");
	std::string ctct_filename = NameBase + std::string("_contacts.dat");

	outParticleData.open(ptcl_filename.c_str(), std::fstream::out);
	outContactData.open(ctct_filename.c_str(), std::fstream::out);
	outParticleData.precision(3);
	outContactData.precision(3);

	outParticleData << std::fixed << "t\t" << "x_pos\t" << "y_pos\t" << "phi\t" << "rad\t" << "mass\t"
					<< "x_vel\t" << "y_vel\t" << "phi_dot\t" << std::endl;

	outContactData << std::fixed << "t\t" << "P_i\t" << "P_j\t" << "x_n_ij\t" << "y_n_ij\t" << "x_t_ij\t" << "y_t_ij\t"
				   << "x_vel\t" << "y_vel\t" << "x_frc\t" << "y_frc\t" << std::endl;

	for (std::vector< std::vector<double> *>::iterator iter  = ParticleData.begin();
			iter != ParticleData.end(); iter++) {
		for(std::vector<double>::iterator jter = (*iter)->begin();
				jter != (*iter)->end(); jter++) {
			outParticleData << (*jter) << "\t";
			}
		outParticleData << std::endl;
		}
	for (std::vector< std::vector<double> *>::iterator iter  = ContactData.begin();
			iter != ContactData.end(); iter++) {
		for(std::vector<double>::iterator jter = (*iter)->begin(); jter != (*iter)->end(); jter++) {
			outContactData << (*jter) << "\t";
			}
		outContactData << std::endl;
		}

	outParticleData.close();
	outContactData.close();
}

