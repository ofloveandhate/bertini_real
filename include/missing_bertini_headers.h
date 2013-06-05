

#include "cascade.h"
#include "polysolve.h"  // the bertini  eval_funcs


#ifndef MISSING_BERTINI_HEADERS_H
#define MISSING_BERTINI_HEADERS_H


void start_system_eval_data_clear_d(start_system_eval_data_d *SSED);//actually lives in bertini library...  testing if this works.

void patch_eval_data_clear_d(patch_eval_data_d *PED);//another which lives in bertini
void patch_eval_data_clear_mp(patch_eval_data_mp *PED);//another which lives in bertini
void changePatchPrec_mp(int new_prec, patch_eval_data_mp *PED); // in bertini


/**
 from the bertini library.  the prototype is not in any header file.
 */
int checkForReal_d(point_d Pt, double realTol);
/**
 from the bertini library.  the prototype is not in any header file.
 */
int checkForReal_mp(point_mp Pt, double realTol);


/**
 from the bertini library.  the prototype is not in any header file.
 */
void findMultSol(post_process_t *endPoints, int num_sols, int num_vars, preproc_data *PPD, double finalTol);



#endif

