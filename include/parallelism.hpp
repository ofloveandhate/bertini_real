
#ifndef BR_PARALLELISM_H
#define BR_PARALLELISM_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpf2mpfr.h>

#include <time.h>
#include <float.h>
#include <limits.h>

#include <sstream>
#include <iostream>






#include "missing_bertini_headers.hpp"

#include "programConfiguration.hpp"


#include "fileops.hpp"
#include "checkSelfConjugate.hpp"
#include "curveNotSelfConj.hpp"
#include "curveSelfConj.hpp"
#include "data_type.hpp"
#include "isosingular.hpp"
#include "programConfiguration.hpp"
#include "witnessSet.hpp"
#include "output.hpp"
#include "surface.hpp"
#include "curve.hpp"
#include "solver_nullspace_left.hpp"






class process
{

public:
	BR_configuration program_options;
	solver_configuration solve_options;
	
	int MPType;
	virtual int main_loop();
};


class ubermaster_process : public process
{
	
public:
	
	ubermaster_process(BR_configuration & new_options, solver_configuration & new_solve_options){
		this->program_options = new_options;
		this->solve_options = new_solve_options;
	}
	
	int main_loop();
	};




class worker_process : public process
{
public:
	
	worker_process(BR_configuration & new_options, solver_configuration & new_solve_options){
		this->program_options = new_options;
		this->solve_options = new_solve_options;
	}
	
	
	
	
	
	int main_loop();

};

#endif





