#include "curve.hpp"



void curve_decomposition::main(vertex_set & V,
                               witness_set & W, // not const, will be modified
                               vec_mp *projections,
                               BR_configuration & program_options,
                               solver_configuration & solve_options)
{
	
#ifdef functionentry_output
	std::cout << "curve::main" << std::endl;
#endif
	
	
		
	
	
	component_num = W.component_number();
	dimension = W.dimension();
	num_variables = W.num_variables();
	
	
	if (1) {
		// perform an isosingular deflation
		// note: somehow, you do not need witness_data to perform isosingular deflation
		if (program_options.verbose_level>=2)
			printf("performing isosingular deflation\n");
		
		
		program_options.input_deflated_filename = W.input_filename();
		
		std::stringstream converter;
		converter << "_dim_" << W.dimension() << "_comp_" << W.component_number() << "_deflated";
		program_options.input_deflated_filename += converter.str();
		converter.clear(); converter.str("");
		
		std::set<unsigned int> zeroonly;
		zeroonly.insert(0);
		W.write_dehomogenized_coordinates("witness_points_dehomogenized",zeroonly); // write the points to file
		int num_deflations, *deflation_sequence = NULL;
		isosingular_deflation(&num_deflations, &deflation_sequence,
							  program_options, W.input_filename(),
							  "witness_points_dehomogenized",
							  program_options.input_deflated_filename,
							  program_options.max_deflations);
		free(deflation_sequence);
		
		
		
		
		
		W.set_input_filename(program_options.input_deflated_filename);
	}
	else {
		program_options.input_deflated_filename = program_options.input_filename;
		//nothing
	}
	
	
	input_filename = W.input_filename();
	
	
	// this wraps around a bertini routine
	parse_input_file(W.input_filename());
	
	preproc_data_clear(&solve_options.PPD);
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	
	int self_conjugate = 1;
	if (W.num_synth_vars()==0) {
		
		if (program_options.verbose_level>=2) {
			printf("checking if component is self-conjugate\n");
		}
		self_conjugate = checkSelfConjugate( *W.point(0), program_options, W.input_filename());  //later:  could be passed in from user, if we want
		
		
		
		if ( verify_projection_ok(W,projections,solve_options) == 1 ){
			if (program_options.verbose_level>=1) {
				printf("verified projection is ok\n");
			}
		}
		else{
			printf("the projection is invalid, in that the jacobian of the randomized system\nbecomes singular at a random point, when the projection is concatenated\n");
			
			br_exit(196);
		}
		
		
	}
	
	
	
	
	if (program_options.user_sphere) {
		read_sphere(program_options.bounding_sphere_filename);
	}
	
	
	
	
	
	
	add_projection(projections[0]);
	
	if (self_conjugate==0)  //C is not self-conjugate
	{
		//Call non-self-conjugate case code
		
		computeCurveNotSelfConj(W, projections[0], V, num_variables,
                                program_options, solve_options);
		
	}
	else
	{
		//Call self-conjugate case code
		
		computeCurveSelfConj(W,
                             projections,
                             V,
                             program_options, solve_options);
	}
	
}











void curve_decomposition::computeCurveSelfConj(const witness_set & W_curve,
                                               vec_mp *projections,
                                               vertex_set &V,
                                               BR_configuration & program_options,
                                               solver_configuration & solve_options)
{
	
#ifdef functionentry_output
	std::cout << "curve::computeCurveSelfConj" << std::endl;
#endif
	
	
    

	
	// 2) randomize down to N-1 equations
	// to get a square system for the homotopies in the following steps.
	
	this->randomizer->setup(W_curve.num_variables()-W_curve.num_patches()-1,solve_options.PPD.num_funcs);
	
	
	
    // 4) solve for critical conditions for random complex projection
	witness_set W_crit_real;
	
	
	compute_critical_points(W_curve,
                            projections,
                            program_options,
                            solve_options,
                            W_crit_real);
	
	
	
	interslice(W_curve,
               W_crit_real,
               projections,
               program_options,
               solve_options,
               V);
	
	return;
} // re: computeCurveSelfConj







//subfunctions
int curve_decomposition::compute_critical_points(const witness_set & W_curve,
                                                 vec_mp *projections,
                                                 BR_configuration & program_options,
                                                 solver_configuration & solve_options,
                                                 witness_set & W_crit_real)
{
#ifdef functionentry_output
	std::cout << "curve::compute_critical_points" << std::endl;
#endif
	
	if (!this->randomizer->is_ready()) {
		std::cout << "randomizer is not setup at compute_critical_points" << std::endl;
		br_exit(29889);
	}
	
	W_crit_real.set_input_filename(W_curve.input_filename());
	
	
	solver_output solve_out;
	
	nullspace_config ns_config;
	compute_crit_nullspace(solve_out, // the returned value
                           W_curve,            // input the original witness set
                           this->randomizer,
                           projections,
                           1,  // dimension of ambient complex object
                           1,   //  target dimension to find
                           1,   // COdimension of the critical set to find.
                           program_options,
                           solve_options,
                           &ns_config);
	ns_config.clear();
	
    
	solve_out.get_noninfinite_w_mult_full(W_crit_real);
	
	
	W_crit_real.only_first_vars(W_curve.num_variables()); // trim the fat, since we are at the lowest level.
	W_crit_real.sort_for_real(&solve_options.T);
	W_crit_real.sort_for_unique(&solve_options.T);
	
	
	if (have_sphere_radius) {
		W_crit_real.sort_for_inside_sphere(sphere_radius, sphere_center);
	}
	else
	{
		std::cout << color::red() << "computing sphere bounds..." << color::console_default() << std::endl;
		compute_sphere_bounds(W_crit_real);
	}
	
	
    witness_set W_additional;
	// now get the sphere intersection critical points and ends of the interval
	get_additional_critpts(&W_additional,  // the returned value
                           W_curve,       // all else here is input
                           program_options,
                           solve_options);
	
    W_additional.sort_for_real(&solve_options.T);
	W_additional.sort_for_unique(&solve_options.T);
	
	
    
    W_crit_real.merge(W_additional);
	
	
	return SUCCESSFUL;
}




int curve_decomposition::get_additional_critpts(witness_set *W_additional,
                                                const witness_set & W_curve,
                                                BR_configuration & program_options,
                                                solver_configuration & solve_options)
{
#ifdef functionentry_output
	std::cout << "curve::get_additional_critpts" << std::endl;
#endif
	
	if (!this->randomizer->is_ready()) {
		std::cout << "randomizer is not ready to go" << std::endl;
		br_exit(13091270);
	}
    
    if (W_curve.num_linears()!=1) {
        std::cout << color::red() << "the input witness set to get_additional_critpts had an incorrect number of linears: " << W_curve.num_linears() << color::console_default() << std::endl;
        br_exit(-518);
    }
    
    
	//build up the start system
	if (program_options.quick_run<=1)
		solve_options.robust = true;
	else
		solve_options.robust = false;
	
	
	
	
	int blabla;
	
	
	
	parse_input_file(W_curve.input_filename(), &blabla);
	preproc_data_clear(&solve_options.PPD); // ugh this sucks
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	multilin_config ml_config(solve_options,this->randomizer); // copies in the randomizer matrix and sets up the SLP & globals.
	
	
	vec_mp *multilin_linears = (vec_mp *) br_malloc(1*sizeof(vec_mp));
	init_vec_mp2(multilin_linears[0],W_curve.num_variables(),solve_options.T.AMP_max_prec);
	multilin_linears[0]->size = W_curve.num_variables();
    
	
	witness_set W_sphere = W_curve;
    //grab just the shell of the input witness set
    W_sphere.reset_points();
    W_sphere.reset_linears();
    W_sphere.reset_patches();
    
    
	sphere_config sp_config(this->randomizer);
	for (int jj=0; jj<W_curve.num_variables(); jj++) {
        set_zero_mp(&multilin_linears[0]->coord[jj]);
    }
	
	for (int ii=0; ii<2; ii++) {
		
		for (int jj=0; jj<W_curve.num_natural_variables(); jj++) {
			get_comp_rand_mp(&multilin_linears[0]->coord[jj]);
		}
		
		vec_cp_mp(sp_config.starting_linear[ii], multilin_linears[0]);
		
		witness_set W_temp;
		
		
		
		
		
		solver_output fillme;
		multilin_solver_master_entry_point(W_curve,         // witness_set
                                           fillme, // the new data is put here!
                                           multilin_linears,
                                           ml_config,
                                           solve_options);
		
		fillme.get_noninfinite_w_mult(W_temp); // should be ordered
		
		W_sphere.merge(W_temp); // copy in the points
		
	}
	
	clear_vec_mp(multilin_linears[0]);
	free(multilin_linears);
	
	// no need to copy in any linears, because the following solve is 0-dimensional.
    // we DO need to copy all the patches from the originating witness set, though.
	W_sphere.copy_patches(W_curve);
	
	
	
	
	
	
	// need to actually move to the sphere system now.
    if (program_options.verbose_level>=1) {
        std::cout << "sphere intersection computation" << std::endl;
    }
    
	
	
	sp_config.set_memory(solve_options); // gets the SLP in memory, and sets up the global memory structures used for evaluation
	sp_config.set_center(this->sphere_center);
	sp_config.set_radius(this->sphere_radius);
	
	
	
	
	
	
	solver_output fillme;
	sphere_solver_master_entry_point(W_sphere,
                                     fillme, // returned value
                                     sp_config,
                                     solve_options);
	
	//get stuff into W_additional from fillme.
	fillme.get_noninfinite_w_mult_full(*W_additional);
	
	
	return 0;
}







int curve_decomposition::interslice(const witness_set & W_curve,
                                    const witness_set & W_crit_real,
                                    vec_mp *projections,
                                    BR_configuration & program_options,
                                    solver_configuration & solve_options,
                                    vertex_set & V)
{
	
#ifdef functionentry_output
	std::cout << "curve::interslice" << std::endl;
#endif
	
	if (!this->randomizer->is_ready()) {
		std::cout << "in interslice, randomizer is not set up properly." << std::endl;
		br_exit(5023);
	}
	
	V.set_curr_projection(projections[0]);
	V.set_curr_input(W_crit_real.input_filename());
	
	this->W = W_curve; // copy in the witness set
	
	this->copy_patches(W_curve);
	
	
	this->num_variables = W_crit_real.num_variables();
	input_filename = W_curve.input_filename();
	
	int blabla;
	parse_input_file(W_curve.input_filename(), &blabla);
	solve_options.get_PPD();
	
	
	
	

	
	///////
	//
	//   actually form crit.
	//
	/////////
	
	
	vertex temp_vertex;
	
	
	std::map<int, int> crit_point_counter;
	
	
	
	
	
	for (unsigned int ii=0; ii<W_crit_real.num_points(); ii++){
		
		if (program_options.verbose_level>=8)
			printf("adding point %u of %zu from W_crit_real to vertices\n",ii,W_crit_real.num_points());
		temp_vertex.set_point( *W_crit_real.point(ii));
		temp_vertex.set_type(CRITICAL); // set type
		
		int I = index_in_vertices_with_add(V, temp_vertex);
		crit_point_counter[I] = 0;
	}
    
    
    
    
    
    
    
	vec_mp crit_downstairs; init_vec_mp(crit_downstairs,0);
	vec_mp midpoints_downstairs; init_vec_mp(midpoints_downstairs,0);
	std::vector< int > index_tracker;
	
    V.compute_downstairs_crit_midpts(W_crit_real,
                                     crit_downstairs,
                                     midpoints_downstairs,
                                     index_tracker, projections[0]);
	

	
	
	int num_midpoints = midpoints_downstairs->size;
	
	if (program_options.verbose_level>=-1) {
		print_point_to_screen_matlab(crit_downstairs,"crit_downstairs");
//		print_point_to_screen_matlab(midpoints_downstairs,"midpoints_downstairs");
	}
	
	
	
	
	
    
	V.set_curr_input(W_curve.input_filename());
	
	int edge_counter = 0; // set the counter
	
	
	std::vector< witness_set> midpoint_witness_sets;
	midpoint_witness_sets.resize(num_midpoints);
	
	
    vec_mp particular_projection;  init_vec_mp(particular_projection,W_curve.num_variables());
		particular_projection->size = W_curve.num_variables();
	vec_cp_mp(particular_projection,projections[0]);
	
	
	

	
    
	
	multilin_config ml_config(solve_options,this->randomizer);
	
	
	
	if (program_options.quick_run<=1)
		solve_options.robust = true;
	else
		solve_options.robust = false;
	
	
	
	
	solve_options.backup_tracker_config();
	
	for (int ii=0; ii<num_midpoints; ii++) {
		
		if (program_options.verbose_level>=2) {
			
			printf("solving midpoints upstairs %d, projection value %lf\n",ii,mpf_get_d(midpoints_downstairs->coord[ii].r));
		}
		
		neg_mp(&particular_projection->coord[0], &midpoints_downstairs->coord[ii]);
		
		
		
		solver_output fillme;
		multilin_solver_master_entry_point(W_curve,         // witness_set
                                           fillme, // the new data is put here!
                                           &particular_projection,
                                           ml_config,
                                           solve_options);
		
		fillme.get_noninfinite_w_mult_full(midpoint_witness_sets[ii]); // is ordered


		midpoint_witness_sets[ii].sort_for_real(&solve_options.T);
		midpoint_witness_sets[ii].sort_for_inside_sphere(sphere_radius, sphere_center);
		
		

		if (program_options.verbose_level>=2) {
			midpoint_witness_sets[ii].print_to_screen();
            std::cout << "midpoint_downstairs " << ii << " had " << midpoint_witness_sets[ii].num_points() << " real points" << std::endl;
		}
		edge_counter += midpoint_witness_sets[ii].num_points();
	}
	
	solve_options.reset_tracker_config();
	
	
    // 7) find edge endpoints
    
	
	
	witness_set Wleft, Wright;
	
	edge temp_edge;
	comp_mp left_proj_val; init_mp(left_proj_val);
	comp_mp right_proj_val; init_mp(right_proj_val);
	
	
	
	std::set<int> found_indices_left;
	std::set<int> found_indices_right;
    
	std::map<int, std::vector< int > > edge_occurence_tracker_left;
	std::map<int, std::vector< int > > edge_occurence_tracker_right;
	
    std::vector< std::set< int > > found_indices_crit;
    std::vector< std::set< int > > found_indices_mid;
    
    found_indices_mid.resize(num_midpoints);
    found_indices_crit.resize(num_midpoints+1);
    
    solve_options.use_gamma_trick = 0;
	
	for (int ii=0; ii<num_midpoints; ++ii) {
		std::cout << color::brown() << "connecting midpoint downstairs " << ii << " of " << num_midpoints << color::console_default() << std::endl;
        
        
        solve_options.backup_tracker_config();
        

        if (program_options.quick_run<=1)
			solve_options.robust = true;
		else
			solve_options.robust = false;
		
		
		
        int keep_going = 1;
        int iterations = 0;
		int maxits = 2;
        while (keep_going==1 && (iterations<maxits))
        {
            
            iterations++;
            keep_going = 0; // we would like to stop computing
            
            
            if (program_options.verbose_level>=2)
			{
                print_comp_matlab(&crit_downstairs->coord[ii],"left ");
				print_comp_matlab(&crit_downstairs->coord[ii+1],"right ");
			}
			
	
			solver_output fillme;
			
			neg_mp(&particular_projection->coord[0], &crit_downstairs->coord[ii]);
            multilin_solver_master_entry_point(midpoint_witness_sets[ii],         // input witness_set
                                               fillme, // the new data is put here!
                                               &particular_projection,
                                               ml_config,
                                               solve_options);
			
			fillme.get_noninfinite_w_mult_full(Wleft); // should be ordered
			
            
			
			fillme.reset();
			
            neg_mp(&particular_projection->coord[0], &crit_downstairs->coord[ii+1]);
            multilin_solver_master_entry_point(midpoint_witness_sets[ii],         // witness_set
                                               fillme, // the new data is put here!
                                               &particular_projection,
                                               ml_config,
                                               solve_options);
			
			fillme.get_noninfinite_w_mult_full(Wright); // should be ordered
			
			
			witness_set Wright_real = Wright; // this is unnecessary
			witness_set Wleft_real = Wleft;   // this is unnecessary
            
			Wright_real.sort_for_real(&solve_options.T);
			Wleft_real.sort_for_real(&solve_options.T);
            
            if (Wleft_real.num_points()!=midpoint_witness_sets[ii].num_points()) {
                std::cout << color::red() << "had a critical failure\n moving left was deficient " << midpoint_witness_sets[ii].num_points()-Wleft_real.num_points() << " points" << color::console_default() << std::endl;
                keep_going = 1;
            }
            
            if (Wright_real.num_points()!=midpoint_witness_sets[ii].num_points()) {
				std::cout << color::red() << "had a critical failure\n moving right was deficient " << midpoint_witness_sets[ii].num_points()-Wright_real.num_points() << " points" << color::console_default() << std::endl;
				keep_going = 1;
            }
            
            if (!keep_going) {
                // this is good, it means we have same number out as in, so we can do a full mapping.
                break;
            }
            else if (iterations<maxits){
              //tighten some tolerances, change it up.
                Wleft.reset();
                Wright.reset();
                std::cout << "trying to recover the failure..." << std::endl;
				
                solve_options.T.endgameNumber = 2;
                // what else can i do here to improve the probability of success?
                solve_options.T.basicNewtonTol   *= 1e-1; // tracktolbeforeeg
                solve_options.T.endgameNewtonTol *= 1e-1; // tracktolduringeg
//				solve_options.T.maxStepSize *= 0.5;
//				solve_options.T.odePredictor = 7;
				std::cout << "tracktolBEFOREeg: "	<< solve_options.T.basicNewtonTol << " tracktolDURINGeg: "	<< solve_options.T.endgameNewtonTol << std::endl;
            }
			else
			{
				Wleft.reset_points();
                Wright.reset_points();
				
				witness_set W_single = midpoint_witness_sets[ii];
				witness_set W_single_sharpened;
				
				
				
				witness_set W_single_right,W_single_left,W_midpoint_replacement = midpoint_witness_sets[ii];
				
				W_midpoint_replacement.reset_points();
				
				
				for (unsigned int kk=0; kk<midpoint_witness_sets[ii].num_points(); kk++) {
					
					W_single.reset_points();
					W_single_sharpened.reset();
					
					
					
					//sharpen up the initial point.
					
					W_single.add_point( *(midpoint_witness_sets[ii].point(kk)));
					
					
					int prev_sharpen_digits = solve_options.T.sharpenDigits;
					solve_options.T.sharpenDigits = MIN(4*solve_options.T.sharpenDigits,300);
					
					neg_mp(&particular_projection->coord[0], &midpoints_downstairs->coord[ii]);
					
					solver_output fillme;
					multilin_solver_master_entry_point(W_single,         // input witness_set
													   fillme, // the new data is put here!
													   &particular_projection,
													   ml_config,
													   solve_options);
					
					fillme.get_noninfinite_w_mult_full(W_single_sharpened);
					fillme.reset();

					
					
					solve_options.T.sharpenDigits = prev_sharpen_digits;
					
					
					
					
					
					//go left and right
					
					comp_mp one_e_minus_four;  init_mp2(one_e_minus_four,1024);
					set_zero_mp(one_e_minus_four);  mpf_set_str(one_e_minus_four->r,"1e-4",10);
					
					comp_mp temp;  init_mp2(temp,1024);
					sub_mp(temp, &crit_downstairs->coord[num_midpoints], &crit_downstairs->coord[0]);
					
					neg_mp(&particular_projection->coord[0], &crit_downstairs->coord[ii]);
					int num_its = 0;
					while (num_its < 1 && W_single_left.num_points()==0) {
						W_single_left.reset();
						
						std::cout << num_its << "th iteration, going left, midpoint " << ii << std::endl;
						
						if (num_its == 0) { // on first try, go default
							solve_options.T.maxNewtonIts = 2;
						}
						else if (num_its == 1) { // on second try, go to straight-up predictor
							solve_options.T.maxNewtonIts = 1;
						}
						else if (num_its == 2){
							solve_options.T.maxNewtonIts = 2;
							sub_mp(&particular_projection->coord[0], &particular_projection->coord[0], temp);
						}
						else{ //try many steps for correction.
							solve_options.T.maxNewtonIts = 4;
							solve_options.T.outputLevel = 3;
						}
						
						solver_output fillme;
						multilin_solver_master_entry_point(W_single_sharpened,         // witness_set
														   fillme, // the new data is put here!
														   &particular_projection,
														   ml_config,
														   solve_options);
						// get stuff from fillme
						fillme.get_noninfinite_w_mult_full(W_single_left);
						fillme.reset();
						
						W_single_left.sort_for_real(&solve_options.T);
						
						num_its++;
					}
					
					
					
					neg_mp(&particular_projection->coord[0], &crit_downstairs->coord[ii+1]);
					num_its = 0;
					
					while (num_its < 1 && W_single_right.num_points()==0) {
						W_single_right.reset();
						
						std::cout << num_its << "th iteration, going right, midpoint " << ii << std::endl;
						
						if (num_its == 0) { // on first try, go default
							solve_options.T.maxNewtonIts = 2;
						}
						else if (num_its == 1) { // on second try, go to straight-up predictor
							solve_options.T.maxNewtonIts = 1;
						}
						else if (num_its == 2){
							solve_options.T.maxNewtonIts = 2;
							add_mp(&particular_projection->coord[0], &particular_projection->coord[0], temp);
						}
						else{ //try many steps for correction.
							solve_options.T.maxNewtonIts = 4;
							solve_options.T.outputLevel = 3;
						}
						solver_output fillme;
						multilin_solver_master_entry_point(W_single_sharpened,         // witness_set
														   fillme, // the new data is put here!
														   &particular_projection,
														   ml_config,
														   solve_options);
						//get stuff from fillme
						fillme.get_noninfinite_w_mult_full(W_single_right);
						fillme.reset();
						
						
						if (num_its==3) {
//							std::cout << "paused for inspecting output file" << std::endl;
//							sleep(600);
						}
						W_single_right.sort_for_real(&solve_options.T);
						
						num_its++;
					}
					
					


					
					
					if (W_single_right.num_points()==1 && W_single_left.num_points()==1) {
						W_midpoint_replacement.add_point( *(midpoint_witness_sets[ii].point(kk)));
						Wleft.add_point(*(W_single_left.point(0)));
						Wright.add_point(*(W_single_right.point(0)));
					}
					else{
						temp_vertex.set_point( *(midpoint_witness_sets[ii].point(kk)) ) ;
						temp_vertex.set_type(PROBLEMATIC); // set type
						index_in_vertices_with_add(V, temp_vertex);
					}
					
					
				}
				
				
				
				midpoint_witness_sets[ii].reset_points();
				midpoint_witness_sets[ii].copy_points(W_midpoint_replacement);
				
			}
			
		}
		
		
        solve_options.reset_tracker_config();
        
		
		for (unsigned int kk=0; kk<midpoint_witness_sets[ii].num_points(); kk++) {
			temp_vertex.set_point( *(midpoint_witness_sets[ii].point(kk)) );
			temp_vertex.set_type(MIDPOINT); // set type
			
			temp_edge.midpt(index_in_vertices_with_add(V, temp_vertex)); // gets the index of the new midpoint as it is added
			
			temp_vertex.set_point( *(Wleft.point(kk)) );
			temp_vertex.set_type(NEW); // set type
			
			temp_edge.left(index_in_vertices_with_add(V, temp_vertex));
			
			
			temp_vertex.set_point( *(Wright.point(kk)) );
			temp_vertex.set_type(NEW); // set type
			
			temp_edge.right(index_in_vertices_with_add(V, temp_vertex));
			
			// keep track of those indices we found.
			
			found_indices_left.insert(temp_edge.left());
			found_indices_right.insert(temp_edge.right());
			
            
            found_indices_crit[ii].insert(temp_edge.left());
            found_indices_crit[ii+1].insert(temp_edge.right());
            found_indices_mid[ii].insert(temp_edge.midpt());

            
			int edge_num = add_edge(temp_edge);
			edge_occurence_tracker_left[temp_edge.left()].push_back(edge_num);
			edge_occurence_tracker_right[temp_edge.right()].push_back(edge_num);
			
			
			// count the number of times a critical point is tracked *to*.  this is for degenerate edge testing.  those which never get tracked to, make degenerate edges (isolated points are included in this).
			if (crit_point_counter.find(temp_edge.left())==crit_point_counter.end()) {
				crit_point_counter[temp_edge.left()] = 1;
			}
			else{
				crit_point_counter[temp_edge.left()] ++;
			}
			
			
			
			if (crit_point_counter.find(temp_edge.right())==crit_point_counter.end()) {
				crit_point_counter[temp_edge.right()] = 1;
			}
			else{
				crit_point_counter[temp_edge.right()] ++;
			}
			
			
			
			
			
			
			
			if (program_options.verbose_level>=2) {
				printf("done connecting upstairs midpoint %d (downstairs midpoint %d)\n",kk,ii);
                
				printf("indices of left, mid, right: %d %d %d\n",temp_edge.left(),temp_edge.midpt(),temp_edge.right());
				printf("\n\n");
			}
		}
		Wleft.reset();
		Wright.reset();
		
        
        
	}//re: for ii
	clear_mp(left_proj_val); clear_mp(right_proj_val);
	
	
    for (int ii=0; ii<num_midpoints; ii++) {
		std::vector<int> bad_crit = V.assert_projection_value(found_indices_crit[ii], &crit_downstairs->coord[ii]);
        std::vector<int> bad_mid = V.assert_projection_value(found_indices_mid[ii], &midpoints_downstairs->coord[ii]);
    }
	std::vector<int> bad_crit = V.assert_projection_value(found_indices_crit[num_midpoints], &crit_downstairs->coord[num_midpoints]);
	
	
	
	
	
	// add degenerate edges for any critical point which is not mapped to by an edge.
	std::map<int,int>::iterator crit_pt_iterator;
	for (crit_pt_iterator = crit_point_counter.begin(); crit_pt_iterator != crit_point_counter.end(); crit_pt_iterator++) {
		int curr_index = crit_pt_iterator->first;
//		int num_occurrences_local = crit_pt_iterator->second;
//		
//		if (num_occurrences_local==0) {
			
            //			vec_cp_mp(temp_vertex.pt_mp, V[curr_index].pt_mp);// set point
            //			projection_value_homogeneous_input(temp_vertex.projVal_mp,  V[curr_index].pt_mp,projections[0]);
            //			temp_vertex.type = ISOLATED; // set type
			
			edge E(curr_index,curr_index,curr_index);
			
			add_edge(E);
//		}
	}
	
	
	
	
	
	
	
	
	
	
	
	/*
	//              \\              //
	//\\\\\\\\\\\\\\\\\            /////////////
	///////////////////  merge    /////////////
	//////////////////            \\\\\\\\\\\\\\
	//             //              \\
	*/
	
	if (program_options.merge_edges==true) {
		this->merge(midpoint_witness_sets[0],V,projections,program_options,solve_options);
	}// re: if merge_edges==true
	else
	{
		
		// since we are not merging, we need to NOT leave the type indicator as NEW, because it may throw off later merges.
		for (std::set<int>::iterator setiter = found_indices_right.begin(); setiter != found_indices_right.end(); setiter++) {
			int curr_index = *setiter;
			if (V[curr_index].type()==NEW) { // only need to look at one of right and left here.
				V[curr_index].set_type(SEMICRITICAL);
			}
		}
		
		for (std::set<int>::iterator setiter = found_indices_left.begin(); setiter != found_indices_left.end(); setiter++) {
			int curr_index = *setiter;
			if (V[curr_index].type()==NEW) { // only need to look at one of right and left here.
				V[curr_index].set_type(SEMICRITICAL);
			}
		}
		
	}//re: else merge_edges==false
	
	
	
	
	
	
	//done
	
	if (program_options.verbose_level>=0) {
		printf("num_edges = %d\n",num_edges);
	}
	
	
	
	clear_vec_mp(particular_projection);
	
	clear_vec_mp(crit_downstairs);
	clear_vec_mp(midpoints_downstairs);
	
	return SUCCESSFUL;
} // re: interslice












//returns <-1> if no candidate found
std::vector<int> curve_decomposition::get_merge_candidate(const vertex_set & V){
	
#ifdef functionentry_output
	std::cout << "curve::get_merge_candidate" << std::endl;
#endif
	
	
	
	std::vector< int > default_found_edges;
	default_found_edges.push_back(-1);
    
	
	// looking for edges with the type NEW, by looking at the left endpoint
	for (int tentative_right_edge=0; tentative_right_edge < this->num_edges; tentative_right_edge++) {
//		std::cout << "looking at edge " << tentative_right_edge << " for merge candidate" << std::endl;
		
		if (V[edges[tentative_right_edge].left()].type() == NEW && V[edges[tentative_right_edge].right()].type() != NEW) {
			// found a starting point for the merges
			
			if (edges[tentative_right_edge].is_degenerate())
				continue; // degenerate edge, should not blabla, but i think hypothetically this will never happen?
			
			std::vector<int> tentative_edge_list;
			tentative_edge_list.push_back(tentative_right_edge);
			
			while (1) {

				
				
				// this goes into an infinite loop if it finds a degenerate edge with the point...
				
//				std::cout << tentative_edge_list.back() << " " << edges[tentative_edge_list.back()].left << " " << edges[tentative_edge_list.back()].midpt << " " << edges[tentative_edge_list.back()].right << std::endl;
				
				int tentative_left_edge = this->nondegenerate_edge_w_right(edges[tentative_edge_list.back()].left());
				
				
				if (tentative_left_edge < 0) {
					std::cout << "found that edge " << tentative_edge_list.back() << " has NEW leftpoint, but \\nexists edge w point " << edges[tentative_edge_list.back()].left() << " as right point." << std::endl;
					break;
					//gotta do something careful here?   i suspect that this happens when two points are very near to each other...
				}
				
//				std::cout << tentative_left_edge << " tent left, with left point type " << V[edges[tentative_left_edge].left].type << std::endl;
				
				tentative_edge_list.push_back(tentative_left_edge);
				
				
				if (V[edges[tentative_left_edge].left()].type() != NEW) {
					break;
				}
				
			}
			
			if (tentative_edge_list.size()>1) {
				return tentative_edge_list;
			}
			else{
				continue;
			}
			
		}
	}
	
	return default_found_edges;
}






void curve_decomposition::merge(witness_set & W_midpt,
                                vertex_set & V,
                                vec_mp * projections,
								BR_configuration & program_options,
                                solver_configuration & solve_options)
{
#ifdef functionentry_output
	std::cout << "curve::merge" << std::endl;
#endif
	
	
	
	vec_mp particular_projection; init_vec_mp(particular_projection,0);
	vec_cp_mp(particular_projection, projections[0]);
	
	comp_mp half;  init_mp2(half,1024);  comp_mp temp;  init_mp2(temp,1024);  comp_mp temp2;  init_mp2(temp2,1024);
	mpf_set_str(half->r, "0.5", 10); mpf_set_str(half->i, "0.0", 10);
	
	
	multilin_config ml_config(solve_options);
	
	
	std::vector< int > edges_to_merge = this->get_merge_candidate(V);
	
	
	comp_mp new_proj_val;  init_mp2(new_proj_val,1024);
	
	
//TODO: parallelize this loop
	
	
	while (edges_to_merge.back()!=-1) { // this value is updated at the end of the loop
		// then there are edges superfluous and need to be merged
		
		int rightmost_edge = edges_to_merge.front();
		int leftmost_edge = edges_to_merge.back();
		
		int moving_edge = edges_to_merge[int(edges_to_merge.size())/2];
		
		if (solve_options.verbose_level>=1) {
			std::cout << color::cyan() << "merging edges: ";
			for (int zz=edges_to_merge.size()-1; zz>=0; zz--) {
				std::cout << edges_to_merge[zz] << " ";
			}
			std::cout << color::console_default() <<  std::endl;
		}
		
		if (edges_to_merge.back() < 0) {
			std::cout << "error: attemping to merge an edge with negative index!" << std::endl;
//			
//			std::cout << "<" << edges[left_edge_w_pt].left << " " << edges[left_edge_w_pt].midpt << " " << edges[left_edge_w_pt].right << "> <";
//			std::cout << edges[right_edge_w_pt].left << " " << edges[right_edge_w_pt].midpt << " " << edges[right_edge_w_pt].right << ">" << std::endl;
//			
			// do something better than break!
			break;                  
		}
		
		
		
		
		
		

		
		witness_set W_temp;
		
		
		// get the projection value of the midpoint we will be moving from.
		//arbitrarily chose to move from the midpoint of the left edge.
		projection_value_homogeneous_input(&particular_projection->coord[0], *V[edges[moving_edge].midpt()].point(), projections[0]);
		neg_mp(&particular_projection->coord[0],&particular_projection->coord[0]);
		
		
		W_midpt.reset_linears(); // this witness set should only have a single linear
		W_midpt.add_linear(particular_projection);
		
		W_midpt.reset_points();
		W_midpt.add_point(*V[edges[moving_edge].midpt()].point());
		// I arbitrarily chose the left edge's midpoint as source to track to new midpoint.
		
		projection_value_homogeneous_input(temp,*V[edges[leftmost_edge].left()].point(),projections[0]);
		projection_value_homogeneous_input(temp2,*V[edges[rightmost_edge].right()].point(),projections[0]);
		
		
		add_mp(new_proj_val, temp, temp2);
		mul_mp(new_proj_val, new_proj_val, half); // now it is the average value
		
		neg_mp(&particular_projection->coord[0], new_proj_val); // set it in the linear for tracking
		
		
		if (program_options.quick_run<=1)
			solve_options.robust = true;
		else
			solve_options.robust = false;
		
		
		
		
		ml_config.set_randomizer(this->randomizer);
		solver_output fillme;
		multilin_solver_master_entry_point(W_midpt,         // witness_set
                                           fillme, // the new data is put here!
                                           &particular_projection,
                                           ml_config,
                                           solve_options);
		
		fillme.get_noninfinite_w_mult_full(W_temp); // should be ordered

		if (W_temp.num_points()==0) {
			std::cout << "merging multilin solver returned NO POINTS!!!" << std::endl;
			continue;
//TODO:  IMMEDIATELY, insert some catch code for when this returns 0 points.
		}

		
        // each member of W_temp should real.  if a member of V already, mark index.  else, add to V, and mark.
		vertex temp_vertex;
		temp_vertex.set_point( *(W_temp.point(0)));
		temp_vertex.set_type(MIDPOINT);
		
		edge temp_edge; // create new empty edge
		
		//set the left, mid and right points
		temp_edge.left(edges[leftmost_edge].left());
		temp_edge.midpt(index_in_vertices_with_add(V, temp_vertex));
		temp_edge.right(edges[rightmost_edge].right());
		
		
		// copy over the removed points for all the edges we are going to merge.

		for (unsigned int zz=0; zz!=edges_to_merge.size(); zz++) {
			int merge_me_away = edges_to_merge[zz];  //set an index into the merge edges
			for (auto vec_iter = edges[merge_me_away].removed_begin(); vec_iter!=edges[merge_me_away].removed_end(); vec_iter++)
			{
				temp_edge.add_removed_point( *vec_iter );
			}
			
			if (zz==0){ // rightmost edge
				temp_edge.add_removed_point(edges[merge_me_away].left());
				temp_edge.add_removed_point(edges[merge_me_away].midpt());
				V[edges[merge_me_away].midpt()].set_removed(true);
				V[edges[merge_me_away].left()].set_removed(true);
			}
			else if (zz==edges_to_merge.size()-1){ // leftmost edge
//				temp_edge.removed_points.push_back(edges[merge_me_away].right());
				temp_edge.add_removed_point(edges[merge_me_away].midpt());
				V[edges[merge_me_away].midpt()].set_removed(true);
//				V[edges[merge_me_away].right()].removed = 1;
			}
			else {
				temp_edge.add_removed_point(edges[merge_me_away].left());
				temp_edge.add_removed_point(edges[merge_me_away].midpt());
				V[edges[merge_me_away].midpt()].set_removed(true);
				V[edges[merge_me_away].left()].set_removed(true);
			}
		}

		add_edge(temp_edge);
		// tacks this onto the end of the edge vector
		
        //		std::cout << "adding edge " << temp_edge.left << " " << temp_edge.midpt << " " << temp_edge.right << " " << std::endl;
		//add the new_edge
		
		
		// delete the old edges // note that we can't do this *IN* the loop because we were using indexes into the edge set.  have to do it after added new edge
		std::vector< edge > post_merge_edges;
		// this should be changed.
		int num_removed_edges = 0;
		for (int ii = 0; ii<this->num_edges; ii++) {
			bool remove_flag = false;
			for (unsigned int zz=0; zz!=edges_to_merge.size(); zz++) {
				if (edges_to_merge[zz] == ii) {
					remove_flag = true;
				}
			}
			
			if (remove_flag==false) { // if don't want to remove the edge // (ii!=left_edge_w_pt) && (ii!=right_edge_w_pt)
				post_merge_edges.push_back( this->edges[ii]);
			}
			else{
				num_removed_edges++;
			}
			//otherwise skip it.
		}
		
		
		if (num_removed_edges!=int(edges_to_merge.size())) {
			std::cout << "claiming to have merged away " << num_removed_edges << " edges, but had " << edges_to_merge.size() << " in the list to merge" << std::endl;
			br_exit(-524);
		}
		
		//swap this's edge info to post-merge info.
		this->edges.swap(post_merge_edges);
		this->num_edges = int(this->edges.size());
		
		edges_to_merge = this->get_merge_candidate(V);
	}// re: while
	clear_mp(half); clear_mp(temp); clear_mp(temp2);
	
	
	clear_mp(new_proj_val);
	clear_vec_mp(particular_projection);
	
} // re: merge






































// will compute a randomizer matrix since you don't provide one. must have current PPD in solve_options for this to work correctly
int verify_projection_ok(const witness_set & W,
                         vec_mp * projection,
                         solver_configuration & solve_options)
{
	
	
	system_randomizer randomizer;
	randomizer.setup(W.num_variables()-W.num_patches()-W.dimension(), solve_options.PPD.num_funcs);
	
	
	int invalid_flag = verify_projection_ok(W, &randomizer, projection, solve_options);

	
	return invalid_flag;
}








int verify_projection_ok(const witness_set & W,
                         system_randomizer * randomizer,
                         vec_mp * projection,
                         solver_configuration & solve_options)
{
	
	int invalid_flag;
	
	parse_input_file(W.input_filename());
	
	
	vec_mp temp_rand_point;  init_vec_mp(temp_rand_point,W.num_variables()); temp_rand_point->size = W.num_variables();
	set_one_mp(&temp_rand_point->coord[0]); // first coordinate must be 1
	for (int ii=1; ii<W.num_variables(); ++ii) {
		get_comp_rand_mp(&temp_rand_point->coord[ii]);
	}
	
	prog_t SLP;
	setupProg(&SLP, solve_options.T.Precision, 2);
	
	
	comp_mp zerotime; init_mp(zerotime);
	set_zero_mp(zerotime);
	
	
	
	eval_struct_mp ED; init_eval_struct_mp(ED, 0, 0, 0);
	evalProg_mp(ED.funcVals, ED.parVals, ED.parDer, ED.Jv, ED.Jp, temp_rand_point, zerotime, &SLP);
	
	mat_mp AtimesJ; init_mat_mp(AtimesJ, 1, 1); AtimesJ->rows = AtimesJ->cols = 1;
		
	randomizer->randomize(temp_rand_point,AtimesJ,ED.funcVals,ED.Jv, &temp_rand_point->coord[0]); // temp_rand_point is the first argument simply to save the creation/deletion of a vec_mp to hold the randomized function values.
	

	mat_mp detme;  init_mat_mp(detme, W.num_variables()-1, W.num_variables()-1);
	detme->cols = detme->rows= W.num_variables()-1;
	
	
	//set the matrix
	for (int ii=0; ii< AtimesJ->rows; ++ii) {
		for (int jj=0; jj<AtimesJ->cols - 1; ++jj) {
			set_mp(&detme->entry[ii][jj],&AtimesJ->entry[ii][jj+1]); // omit the homogeneous coordinate's columns
		}
	}
	
	int offset = W.num_variables()-1 - W.dimension();
	for (int jj=0; jj < W.dimension(); jj++){
		for (int ii=0; ii<W.num_variables()-1; ii++) {
			set_mp(&detme->entry[offset+jj][ii], &projection[jj]->coord[ii+1]);
		}
	}
	
	
	comp_mp determinant; init_mp(determinant);
	take_determinant_mp(determinant,detme);
	
	if ( d_abs_mp(determinant) < 1e-2){
		invalid_flag = 0;
		std::cout << d_abs_mp(determinant) << "\n";
		print_matrix_to_screen_matlab(ED.Jv,"Jv");
		print_matrix_to_screen_matlab(detme,"detme");
	}
	else
		invalid_flag = 1;
	
	clear_mat_mp(detme);
	clear_mat_mp(AtimesJ);
	clear_vec_mp(temp_rand_point);
	clear_mp(determinant);
	clear_mp(zerotime);
	
	clear_eval_struct_mp(ED);
	
	clearProg(&SLP, solve_options.T.MPType, 1);
	
	return invalid_flag;
	
}









void curve_decomposition::send(int target, parallelism_config & mpi_config)
{
#ifdef functionentry_output
	std::cout << "curve::send" << std::endl;
#endif
	
	
	
	decomposition::send(target, mpi_config);
	
	MPI_Send(&num_edges, 1, MPI_INT, target, CURVE, MPI_COMM_WORLD);
	for (int ii=0; ii<num_edges; ii++) {
		edges[ii].send(target, mpi_config);
	}
	
}



void curve_decomposition::receive(int source, parallelism_config & mpi_config)
{
#ifdef functionentry_output
	std::cout << "curve::receive" << std::endl;
#endif
	
	
	
	decomposition::receive(source, mpi_config);
	
	MPI_Status statty_mc_gatty;
	
	int temp_num_edges;
	MPI_Recv(&temp_num_edges, 1, MPI_INT, source, CURVE, MPI_COMM_WORLD, &statty_mc_gatty);
	
	
	for (int ii=0; ii<temp_num_edges; ii++) {
		edge E;
		E.receive(source, mpi_config);
		add_edge(E);
	}
	
}






int curve_decomposition::setup(boost::filesystem::path containing_folder){
	decomposition::setup(containing_folder / "decomp");
	
	setup_edges(containing_folder / "E.edge");
	
	
	return 1;
}





int curve_decomposition::setup_edges(boost::filesystem::path INfile)
//setup the vertex structure
{
#ifdef functionentry_output
	std::cout << "curve::setup_edges" << std::endl;
#endif
	
	
	FILE *IN = safe_fopen_read(INfile);
	
	fscanf(IN, "%d\n", &this->num_edges);
	int left, midpt, right;
	for(int ii=0;ii<this->num_edges;ii++) {
		fscanf(IN,"%d %d %d",&left, &midpt, &right); scanRestOfLine(IN);
		this->edges.push_back(edge(left, midpt, right));
	}
	
	fclose(IN);
	return this->num_edges;
}








void curve_decomposition::print(boost::filesystem::path base)
{
#ifdef functionentry_output
	std::cout << "curve::print" << std::endl;
#endif
	
    //	std::cout << "printing curve decomposition to folder " << base << std::endl;
	
	decomposition::print(base);
	
	boost::filesystem::path edgefile = base / "E.edge";
	
	curve_decomposition::print_edges(edgefile);
	
}







void curve_decomposition::print_edges(boost::filesystem::path outputfile)
{
#ifdef functionentry_output
	std::cout << "curve::print_edges" << std::endl;
#endif
	
	int ii;
	FILE *OUT = safe_fopen_write(outputfile);
	
	// output the number of vertices
	fprintf(OUT,"%d\n\n",num_edges);
	
	for(ii=0;ii<num_edges;ii++)
		fprintf(OUT,"%d %d %d \n",
                edges[ii].left(),
                edges[ii].midpt(),
                edges[ii].right() );
	fclose(OUT);
}










void curve_decomposition::computeCurveNotSelfConj(const witness_set & W_in,
                                                  vec_mp         pi_mp,
                                                  vertex_set			&V,
                                                  int           num_vars,
                                                  BR_configuration & program_options,
                                                  solver_configuration & solve_options)

{
#ifdef functionentry_output
	std::cout << "curve::computeCurveNotSelfConj" << std::endl;
#endif
	
	
	
	// num_vars includes the homogeneous variable
	
	
    
	std::string bertini_system_command = program_options.bertini_command;
	bertini_system_command.append(" input_NSC start_NSC");
	
    FILE *IN = NULL;
    
	
	int *declarations = NULL;
    partition_parse(&declarations, W_in.input_filename(), "func_input_nsc", "config_nsc",1);
	// optional:  here perform sanity checks
	free(declarations);
	
	
	
	
    //generate input file
    diag_homotopy_input_file("input_NSC", "func_input_nsc","func_inputbar","config_nsc",*W_in.linear(0),num_vars-1);
    //generate start file
	diag_homotopy_start_file("start_NSC",  W_in);
	
	
    //run bertini
	
	copyfile("witness_data","witness_data_0");
	
    int retVal = system(bertini_system_command.c_str());
	if (retVal!=0) {
		std::cout << "system calling bertini for diagonal intersection returned nonzero value!" << std::endl;
		br_exit(retVal);
	}
	
	rename("witness_data_0","witness_data");
	
	
    //read the real solutions
    IN = safe_fopen_read("real_solutions");
    
	int num_sols;
    fscanf(IN, "%d\n\n", &num_sols);
	
	
	vertex temp_vertex;
	change_size_vec_mp(*temp_vertex.point(),num_vars); (*temp_vertex.point())->size = num_vars;
	temp_vertex.set_type(ISOLATED);
	
	
	vec_mp cur_sol,cur_sol_bar;
    init_vec_mp(cur_sol,num_vars); cur_sol->size = num_vars;
	set_one_mp(&cur_sol->coord[0]);
	
    init_vec_mp(cur_sol_bar,num_vars); cur_sol_bar->size = num_vars;
	set_one_mp(&cur_sol_bar->coord[0]);
	
	
	for (int ii=0; ii<num_sols; ii++) {
		for (int jj=0; jj<num_vars-1; jj++){
			mpf_inp_str(cur_sol->coord[jj+1].r, IN, 10);
			mpf_inp_str(cur_sol->coord[jj+1].i, IN, 10);
			
			mpf_inp_str(cur_sol_bar->coord[jj+1].r, IN, 10);
			mpf_inp_str(cur_sol_bar->coord[jj+1].i, IN, 10);
		}
        
        //check if x=x_bar
		
        if (isSamePoint_homogeneous_input(cur_sol,cur_sol_bar)) { // x=x_bar
			temp_vertex.set_point(cur_sol);
			
			index_in_vertices_with_add(V, temp_vertex);
			
		}
	}
	
    //set the number of vertices
    fclose(IN);
    //clear
    clear_vec_mp(cur_sol);
    clear_vec_mp(cur_sol_bar);
    
    
	// delete temporary files
    remove("func_input_nsc");
    remove("config_nsc");
    remove("func_inputbar");
    remove("var_names");
	
	
}




void diag_homotopy_input_file(boost::filesystem::path outputFile,
                              boost::filesystem::path funcInputx,
                              boost::filesystem::path funcInputy,
                              boost::filesystem::path configInput,
                              vec_mp L,
                              int   num_vars)
/***************************************************************\
 * USAGE: setup input file to do diagonal homotopy             *
 * ARGUMENTS: name of output file, function & configuration input*
 * RETURN VALUES: none                                           *
 * NOTES:                                                        *
 \***************************************************************/
{
    char ch,**str,*fmt = NULL;
    int ii,jj,size;
    mat_d A;
	
	FILE *IN = NULL;
	
	
    str=(char **)br_malloc(num_vars*sizeof(char *));
    for(ii=0;ii<num_vars;ii++)
        str[ii]=(char*)br_malloc(sizeof(char)*256);
    
	FILE *OUT = safe_fopen_write(outputFile);
	
    // setup configurations in OUT
    fprintf(OUT, "CONFIG\n");
    IN = safe_fopen_read(configInput);
    fclose(IN);
	
    fprintf(OUT, "USERHOMOTOPY: 1;\nDeleteTempFiles: 0;\nEND;\nINPUT\n");
	
    // setup variables in OUT
    IN = safe_fopen_read(funcInputx);
    while ((ch = fgetc(IN)) != EOF )
        fprintf(OUT, "%c", ch);
    fclose(IN);
	
	//setup the function name in OUT
	IN = safe_fopen_read(funcInputy);
    while ((ch = fgetc(IN)) != EOF )
        fprintf(OUT, "%c", ch);
    fclose(IN);
	
	
    IN = safe_fopen_read("var_names");
    ii=0;jj=0;
    while ((ch = fgetc(IN)) != EOF)
    {
        if(ch!='\n')
            str[ii][jj++]=ch;
        else
        {
            str[ii++][jj]='\0';
            jj=0;
        }
    }
	
    //setup the linear equations
    // find the size needed
    size = 1 + snprintf(NULL, 0, "%%.%dlf+%%.%dlf*I", 15, 15);
    // allocate size
    fmt = (char *)br_malloc(size * sizeof(char));
    // setup fmt
    sprintf(fmt, "%%.%dlf+%%.%dlf*I", 15, 15);
    // output the linear function L and L_bar
	
	comp_mp temp;  init_mp2(temp,L->curr_prec);
    for (ii = 0; ii < L->size; ii++)
    {
        fprintf(OUT, "bertini_real_L%d = ",ii);
        // print output
		mpf_out_str(OUT,10,0,L->coord[ii].r);
		fprintf(OUT,"+I*");
		mpf_out_str(OUT,10,0,L->coord[ii].i);
		fprintf(OUT, ";\n");
		
        fprintf(OUT, "bertini_real_Lbar%d = ",ii);
		conjugate_mp(temp, &L->coord[ii]);
        // print output
		mpf_out_str(OUT,10,0,temp->r);
		fprintf(OUT,"+I*");
		mpf_out_str(OUT,10,0,temp->i);
        fprintf(OUT, ";\n");
    }
    fprintf(OUT, "\n\n");
	clear_mp(temp);
	
    //Generate a random matrix A and output to input file.
    init_mat_d(A, 2, num_vars);
	make_matrix_random_d(A, 2, num_vars);
    for (ii = 0; ii < 2; ii++)
        for(jj=0;jj<num_vars;jj++)
        {
            fprintf(OUT, "A%d%d = ",ii,jj);
            // print output
            fprintf(OUT, fmt, A->entry[ii][jj].r, A->entry[ii][jj].i);
            fprintf(OUT, ";\n");
        }
    //setup the diagonal homotopy functions
    fprintf(OUT, "\nbertini_real_L=t*(");
    //(Lx-1)*t+(1-t)*A[0]*(x-x_bar)
    for(ii=0;ii<num_vars;ii++)
    {
        fprintf(OUT, "bertini_real_L%d*", ii);
        fprintf(OUT, "%s", str[ii]);
        fprintf(OUT, "+");
    }
    fprintf(OUT, "-1)+(1-t)*(");
    for(ii=0;ii<num_vars;ii++)
    {
        fprintf(OUT, "A0%d*", ii);
        fprintf(OUT, "(%s-%s", str[ii],str[ii]);
        fprintf(OUT, "bar)+");
    }
    //(L_bar x_bar-1)*t+(1-t)*A[1]*(x-x_bar)
    fprintf(OUT, "0);\nbertini_real_Lbar=t*(");
    for(ii=0;ii<num_vars;ii++)
    {
        fprintf(OUT, "bertini_real_Lbar%d*", ii);
        fprintf(OUT, "%sbar", str[ii]);
        fprintf(OUT, "+");
    }
    fprintf(OUT, "-1)+(1-t)*(");
    for(ii=0;ii<num_vars;ii++)
    {
        fprintf(OUT, "A1%d*", ii);
        fprintf(OUT, "(%s-%s", str[ii],str[ii]);
        fprintf(OUT, "bar)+");
    }
    fprintf(OUT, "0);\nEND;");
    fclose(OUT);
    //free
    for(ii=0;ii<num_vars;ii++)
        free(str[ii]);
    free(str);
    free(fmt);
    clear_mat_d(A);
    return;
}



void diag_homotopy_start_file(boost::filesystem::path startFile,
                              const witness_set & W)
/***************************************************************\
 * USAGE: setup start file to do diagonal homotopy             *
 * ARGUMENTS: name of output file, start points & number of variables*
 * RETURN VALUES: none                                           *
 * NOTES:                                                        *
 \***************************************************************/
{
    

	FILE *OUT = safe_fopen_write(startFile);
    

	// output the number of start points
    fprintf(OUT,"%zu\n\n",W.num_points()*W.num_points());
	
	
	comp_mp temp; init_mp(temp);
	vec_mp result; init_vec_mp(result,0);
	vec_mp result2; init_vec_mp(result2,0);
	
	for (unsigned int ii=0; ii<W.num_points(); ii++){
		
		vec_mp * outer_point = W.point(ii);
		
		change_prec_vec_mp(result,(*outer_point)->curr_prec);
		dehomogenize(&result,*outer_point);
		
		for (unsigned int jj=0; jj<W.num_points(); jj++) { // output {w \bar{w}}'
			
			vec_mp * inner_point = W.point(jj);
			
			change_prec_vec_mp(result,(*inner_point)->curr_prec);
			dehomogenize(&result2,*inner_point);
			
			change_prec_mp(temp,(*inner_point)->curr_prec);
			
			for(int kk=0; kk<W.num_variables()-1;kk++) {
				print_mp(OUT, 0, &result->coord[kk]); fprintf(OUT, "\n");
				
				conjugate_mp(temp, &result2->coord[kk] )
				print_mp(OUT, 0, temp); fprintf(OUT, "\n");
			} // re: kk
			fprintf(OUT,"\n");

		}// re: jj
	}// re: ii
	
	clear_vec_mp(result2);
	clear_vec_mp(result);
	clear_mp(temp);

    fclose(OUT);
}































