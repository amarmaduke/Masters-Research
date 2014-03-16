#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "../utility/def.h"
#include "parameter.h"
#include "../json/json.h"
#include <iostream>
#include <cstdio>

#ifndef JACOBI_H
#define JACOBI_H

typedef double value_type;
typedef boost::numeric::ublas::vector<value_type> vector_type;
typedef boost::numeric::ublas::matrix<value_type> matrix_type;

struct jacobi_functor
{

	parameter& state;

	jacobi_functor(parameter& p) : state(p) { }

	void operator() (	const vector_type&,
										matrix_type&,
										const value_type&,
										vector_type&);
};


#endif

