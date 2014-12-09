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
class Process
{
protected:
	BertiniRealConfig program_options;///< holds the current state of Bertini_real
	SolverConfiguration solve_options; ///< holds the current state of the solver
	
	
public:

	
	int MPType; ///< operating MP type.
	
	
	virtual int main_loop() = 0;
	
	virtual ~Process()
	{
		
	};
	
	
	
};



/**
 \brief Master process, level 0.
 */
class UbermasterProcess : public Process
{
	
public:
	
	UbermasterProcess(BertiniRealConfig & new_options, SolverConfiguration & new_solve_options){
		this->program_options = new_options;
		this->solve_options = new_solve_options;
	}
	
	
	/**
	 \brief Master Bertini_real procedure.
	 
	 Loads the witness_data, tracker config, and decomposes components of the user's choosing.
	 \return An integer flag indicating the success of the loop.
	 */
	int main_loop();

	void bertini_real(WitnessSet & W, vec_mp *pi, VertexSet & V);
	
	
	void critreal(WitnessSet & W, vec_mp *pi, VertexSet & V);
	
	
	
	
	~UbermasterProcess()
	{
		
	}
};



/**
 \brief worker process, level 1.
 */
class WorkerProcess : public Process
{
public:
	
	WorkerProcess(BertiniRealConfig & new_options, SolverConfiguration & new_solve_options){
		this->program_options = new_options;
		this->solve_options = new_solve_options;
	}
	
	
	
	~WorkerProcess()
	{
		
	}
	
	/**
	 \brief the listen-work loop for workers in an MPI ring.
	 
	 \return SUCCESSFUL flag from process.
	 */
	int main_loop();

};




/**
 \brief reads in projection from file if user specified, creates one otherwise.

 currently defaults to create a random real projection with homogeneous value 0;

 \param pi the projection vectors to fill.  must be initted already, but not necessarily the correct size.
 \param program_options The current state of Bertini_real.
 \param num_vars how many variables to set up, including the homogenizing variable.
 \param num_projections how many proj vectors to set up.  again, these must already be allocated outside this call.
 */
void get_projection(vec_mp *pi,
					BertiniRealConfig program_options,
					int num_vars,
					int num_projections);




#endif





