/*
 * Exceptions.h
 *
 *  Created on: Feb 7, 2011
 *      Author: mixlmay
 */

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#include <stdexcept>

class BadConversion : public std::runtime_error {
 public:
   BadConversion(std::string const& s)
     : std::runtime_error(s)
     { }
 };

class InvalidDataReqest: public std::runtime_error {
public:
	InvalidDataReqest(std::string const& s)
    : std::runtime_error(s)
    { }
};

class Exception {
protected:
	static unsigned int current_id;
	unsigned int id;

public:
	Exception() {   }
	Exception(const Exception &src){
		id = src.id;
	}
	virtual ~Exception(){}
};

class slowConvergenceException : virtual Exception {
protected:
	unsigned int NumberOfIterations;
public:
	slowConvergenceException() {   }
	slowConvergenceException(slowConvergenceException &src) { NumberOfIterations = src.NumberOfIterations; }
	virtual ~slowConvergenceException() {   }
};

#endif /* EXCEPTIONS_H_ */
