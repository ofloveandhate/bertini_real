#ifndef BR_PARALLELISM_H
#define BR_PARALLELISM_H



/** \file parallelism.hpp */


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






#include "bertini_headers.hpp"

#include "programConfiguration.hpp"


#include "fileops.hpp"
#include "checkSelfConjugate.hpp"

#include "data_type.hpp"
#include "isosingular.hpp"
#include "programConfiguration.hpp"

#include "surface.hpp"
#include "curve.hpp"


/**
\defgroup mpienabled MPI-enabled classes

 */

/**
 \defgroup send MPI Sends
 */

/**
 \defgroup receive MPI Receives
 */

/**
 \defgroup bcast MPI Bcasts (sends and receives)
 */


/**
 \brief BR process base class -- holds current state of program and solver.
 
 In order to make the program_options and solve_options 'globally' accessibly, we place them into the containing process.
 */
class process
{
protected:
	BR_configuration program_options;///< holds the current state of Bertini_real
	solver_configuration solve_options; ///< holds the current state of the solver
	
	
public:

	
	int MPType; ///< operating MP type.
	
	
	virtual int main_loop() = 0;
	
	virtual ~process()
	{
		
	};
	
	
	
};



/**
 \brief Master process, level 0.
 */
class ubermaster_process : public process
{
	
public:
	
	ubermaster_process(BR_configuration & new_options, solver_configuration & new_solve_options){
		this->program_options = new_options;
		this->solve_options = new_solve_options;
	}
	
	
	/**
	 \brief Master Bertini_real procedure.
	 
	 Loads the witness_data, tracker config, and decomposes components of the user's choosing.
	 \return An integer flag indicating the success of the loop.
	 */
	int main_loop();

	void critreal(witness_set & W, vec_mp *pi, vertex_set & V);
	~ubermaster_process()
	{
		
	}
};



/**
 \brief worker process, level 1.
 */
class worker_process : public process
{
public:
	
	worker_process(BR_configuration & new_options, solver_configuration & new_solve_options){
		this->program_options = new_options;
		this->solve_options = new_solve_options;
	}
	
	
	
	~worker_process()
	{
		
	}
	
	/**
	 \brief the listen-work loop for workers in an MPI ring.
	 
	 \return SUCCESSFUL flag from process.
	 */
	int main_loop();

};










#endif





