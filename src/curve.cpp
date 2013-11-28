#include "curve.hpp"



void curve_decomposition::main(vertex_set & V,
                               witness_set & W, // not const, will be modified
                               vec_mp *projections,
                               BR_configuration & program_options,
                               solver_configuration & solve_options)
{
	
	int num_vars = W.num_variables;
	
	solve_options.robust = false;
	
	
	
	// perform an isosingular deflation
	// note: you do not need witness_data to perform isosingular deflation
	if (program_options.verbose_level>=2)
		printf("performing isosingular deflation\n");
	
	W.write_dehomogenized_coordinates("witness_points_dehomogenized"); // write the points to file
	int num_deflations, *deflation_sequence = NULL;
	isosingular_deflation(&num_deflations, &deflation_sequence,
                          program_options, W.input_filename,
                          "witness_points_dehomogenized",
                          program_options.max_deflations,
                          W.dim, W.comp_num);
	free(deflation_sequence);
	
	
	
	program_options.input_deflated_filename = W.input_filename;
	
	std::stringstream converter;
	converter << "_dim_" << W.dim << "_comp_" << W.comp_num << "_deflated";
	program_options.input_deflated_filename += converter.str();
	converter.clear(); converter.str("");
	
	W.input_filename = program_options.input_deflated_filename;
	
	// this wraps around a bertini routine
	parse_input_file(W.input_filename);
	
	preproc_data_clear(&solve_options.PPD);
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	
	int self_conjugate = 1;
	if (W.num_synth_vars==0) {
		
		if (program_options.verbose_level>=2) {
			printf("checking if component is self-conjugate\n");
		}
		self_conjugate = checkSelfConjugate(W, num_vars, program_options, W.input_filename);  //later:  could be passed in from user, if we want
		
		//regenerate the various files, since we ran bertini since then.
		parse_input_file(W.input_filename);
		
		if (verify_projection_ok(W,
                                 projections,
                                 solve_options)==1){
			if (program_options.verbose_level>=1) {
				printf("verified projection is ok\n");
			}
		}
		else{
			printf("the projection is invalid, in that the jacobian of the randomized system\nbecomes singular at a random point, when the projection is concatenated\n");
			
			print_point_to_screen_matlab(projections[0], "projections[0]");
			
			br_exit(196);
		}
		
	}
	
	
	
	
	input_filename = W.input_filename;
	component_num = W.comp_num;
	dimension = W.dim;
	num_variables = num_vars;
	add_projection(projections[0]);
	
	if (self_conjugate==0)  //C is not self-conjugate
	{
		//Call non-self-conjugate case code
		
		computeCurveNotSelfConj(W, projections[0], V, num_vars,
                                program_options, solve_options);
		
	}
	else
	{
		//Call self-conjugate case code
		
		computeCurveSelfConj(W,
                             projections,
                             V,
                             num_vars,
                             program_options, solve_options);
	}
	
}











void curve_decomposition::computeCurveSelfConj(const witness_set & W_curve,
                                               vec_mp *projections,
                                               vertex_set &V,
                                               int num_vars,
                                               BR_configuration & program_options,
                                               solver_configuration & solve_options)
{
	//IN DEVELOPMENT
    
	
	
	int ambient_dim = 1;
	
	witness_set Wtemp;
	
	
	
	
	
	// 2) randomize down to N-1 equations
	// to get a square system for the homotopies in the following steps.
	
	
	//create the matrix
	init_mat_mp2(this->randomizer_matrix,
                 W_curve.num_variables-W_curve.num_patches-ambient_dim,solve_options.PPD.num_funcs,
                 solve_options.T.AMP_max_prec);
	
	
	//get the matrix and the degrees of the resulting randomized functions.
	make_randomization_matrix_based_on_degrees(this->randomizer_matrix, this->randomized_degrees, W_curve.num_variables-W_curve.num_patches-ambient_dim, solve_options.PPD.num_funcs);
	
	if (program_options.verbose_level>=4)
		print_matrix_to_screen_matlab(randomizer_matrix,"randomization");
	
	FILE *OUT = safe_fopen_write("Rand_matrix");
	fprintf(OUT,"%d %d\n",randomizer_matrix->rows,randomizer_matrix->cols);
	fprintf(OUT,"%d\n",solve_options.T.AMP_max_prec);
	fprintf(OUT,"\n");
	print_matrix_to_file_mp(OUT, 0, randomizer_matrix);
	fclose(OUT);
	
	
	
    // 4) solve for critical conditions for random complex projection
	witness_set W_crit_real;
	
	
	compute_critical_points(W_curve,
                            this->randomizer_matrix,
                            this->randomized_degrees,
                            projections,
                            program_options,
                            solve_options,
                            W_crit_real);
	
	
	
	interslice(W_curve,
               W_crit_real,
               randomizer_matrix,
               projections,
               program_options,
               solve_options,
               V);
	
	return;
} // re: computeCurveSelfConj








int curve_decomposition::interslice(const witness_set & W_curve,
                                    const witness_set & W_crit_real,
                                    mat_mp randomizer_matrix,
                                    vec_mp *projections,
                                    BR_configuration & program_options,
                                    solver_configuration & solve_options,
                                    vertex_set & V)
{
	
	V.set_curr_projection(projections[0]);
	
	
	
	
	for (int ii=0; ii<W_curve.num_patches; ii++)
		this->add_patch(W_curve.patch_mp[ii]);
	
	this->num_variables = W_crit_real.num_variables;
	input_filename = W_curve.input_filename;
	
	int blabla; int *declarations = NULL;
	
	parse_input_file(W_curve.input_filename, &blabla);
	solve_options.get_PPD();
	
	
	std::vector< int > randomized_degrees;
	//get the matrix and the degrees of the resulting randomized functions.
	make_randomization_matrix_based_on_degrees(randomizer_matrix, randomized_degrees, W_curve.num_variables-W_curve.num_patches-W_curve.num_linears, solve_options.PPD.num_funcs);
	
	
	///////
	//
	//   actually form crit.
	//
	/////////
	
	
	vertex temp_vertex;
	
	
	std::map<int, int> crit_point_counter;
	
	
	V.set_curr_input(W_crit_real.input_filename);
	
	
	for (int ii=0; ii<W_crit_real.num_pts; ii++){
		
		if (program_options.verbose_level>=8)
			printf("adding point %d of %d from W_crit_real to vertices\n",ii,W_crit_real.num_pts);
		
		vec_cp_mp(temp_vertex.pt_mp, W_crit_real.pts_mp[ii]);// set point
		temp_vertex.type = CRITICAL; // set type
		
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
		print_point_to_screen_matlab(midpoints_downstairs,"midpoints_downstairs");
	}
	
	
	
	
	
    
	V.set_curr_input(W_curve.input_filename);
	
	int edge_counter = 0; // set the counter
	
	
	std::vector< witness_set> midpoint_witness_sets;
	midpoint_witness_sets.resize(num_midpoints);
	
	
    vec_mp particular_projection;  init_vec_mp(particular_projection,W_curve.num_variables); particular_projection->size = W_curve.num_variables;
	vec_cp_mp(particular_projection,projections[0]);
	
	
	
	
	solve_options.allow_multiplicity = 0;
	solve_options.allow_singular = 0;
	solve_options.complete_witness_set = 1;
	
    
	
    //	print_matrix_to_screen_matlab(randomizer_matrix,"rand_interslice");
	multilin_config ml_config(solve_options,randomizer_matrix);
	
	solve_options.robust = true;
	solve_options.backup_tracker_config();
	
	for (int ii=0; ii<num_midpoints; ii++) {
		
		if (program_options.verbose_level>=2) {
			
			printf("solving midpoints upstairs %d, projection value %lf\n",ii,mpf_get_d(midpoints_downstairs->coord[ii].r));
		}
		
		neg_mp(&particular_projection->coord[0], &midpoints_downstairs->coord[ii]);
		
		
		solve_options.complete_witness_set = 1;
		
		multilin_solver_master_entry_point(W_curve,         // witness_set
                                           &midpoint_witness_sets[ii], // the new data is put here!
                                           &particular_projection,
                                           ml_config,
                                           solve_options);

		
		midpoint_witness_sets[ii].sort_for_real(solve_options.T);
		midpoint_witness_sets[ii].sort_for_unique(solve_options.T);
		midpoint_witness_sets[ii].sort_for_inside_sphere(sphere_radius, sphere_center);
		

		if (program_options.verbose_level>=2) {
			midpoint_witness_sets[ii].print_to_screen();
            std::cout << "midpoint_downstairs " << ii << " had " << midpoint_witness_sets[ii].num_pts << " real points" << std::endl;
		}
		edge_counter += midpoint_witness_sets[ii].num_pts;
	}
	
	solve_options.reset_tracker_config();
	
	
    // 7) find edge endpoints
    
	solve_options.allow_multiplicity = 1;
	solve_options.allow_singular = 1;
	
	
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
    
    
	for (int ii=0; ii<num_midpoints; ++ii) {
		std::cout << color::brown() << "connecting midpoint downstairs " << ii << color::console_default() << std::endl;
        
        
        solve_options.backup_tracker_config();
        

        solve_options.robust = true;
        int keep_going = 1;
        int iterations = 0;
        while (keep_going==1 && (iterations<3))
        {
            
            iterations++;
            keep_going = 0; // we would like to stop computing
            
            neg_mp(&particular_projection->coord[0], &crit_downstairs->coord[ii]);
            if (program_options.verbose_level>=2)
                print_comp_matlab(&crit_downstairs->coord[ii],"left ");
            multilin_solver_master_entry_point(midpoint_witness_sets[ii],         // witness_set
                                               &Wleft, // the new data is put here!
                                               &particular_projection,
                                               ml_config,
                                               solve_options);
            
            
            
            neg_mp(&particular_projection->coord[0], &crit_downstairs->coord[ii+1]);
            
            if (program_options.verbose_level>=2)
                print_comp_matlab(&crit_downstairs->coord[ii+1],"right ");
            
            multilin_solver_master_entry_point(midpoint_witness_sets[ii],         // witness_set
                                               &Wright, // the new data is put here!
                                               &particular_projection,
                                               ml_config,
                                               solve_options);
            
//            if (iterations<=1) {
//                Wleft.sort_for_real(solve_options.T);
//                Wright.sort_for_real(solve_options.T);
//            }
            
            
            
            if (Wleft.num_pts!=midpoint_witness_sets[ii].num_pts) {
                std::cout << color::red() << "had a critical failure\n moving left was deficient a point" << color::console_default() << std::endl;
                keep_going = 1;
            }
            
            if (Wright.num_pts!=midpoint_witness_sets[ii].num_pts) {
                std::cout << color::red() << "had a critical failure\n moving right was deficient a point" << color::console_default() << std::endl;
                keep_going = 1;
            }
            
            if (!keep_going) {
                
                break;
            }
            else if (iterations<3) {
              //tighten some tolerances, change it up.
                Wleft.reset();
                Wright.reset();
                std::cout << "trying to recover the failure...";
                solve_options.T.endgameNumber = 2;
                // what else can i do here to improve the probability of success?
                solve_options.T.basicNewtonTol   *= 1e-6;
                solve_options.T.endgameNewtonTol *= 1e-6;
            }
			else // iterations == 3
			{

			}


        }
        solve_options.reset_tracker_config();
        
		
		for (int kk=0; kk<midpoint_witness_sets[ii].num_pts; kk++) {
			
			if (Wleft.num_pts!=midpoint_witness_sets[ii].num_pts || Wright.num_pts!=midpoint_witness_sets[ii].num_pts)
			{
				break;
			}
			
			
			bool ok = true;
			
			vec_mp result; init_vec_mp(result,W_curve.num_variables-W_curve.num_synth_vars-1);
			result->size = W_curve.num_variables-W_curve.num_synth_vars-1;
			
			
			for (int jj=1; jj<W_curve.num_variables-W_curve.num_synth_vars; jj++) {
				div_mp(&result->coord[jj-1], &Wleft.pts_mp[kk]->coord[jj], &Wleft.pts_mp[kk]->coord[0]);
			}
            
			if (!checkForReal_mp(result, solve_options.T.real_threshold)) {
				std::cout << color::red();
				print_point_to_screen_matlab(result,"not_real_but_should_be_left");
                print_point_to_screen_matlab(Wleft.pts_mp[kk],"nat");
				std::cout << color::console_default();
				ok = false;
			}
			
			
			for (int jj=1; jj<W_curve.num_variables-W_curve.num_synth_vars; jj++) {
				div_mp(&result->coord[jj-1], &Wright.pts_mp[kk]->coord[jj], &Wright.pts_mp[kk]->coord[0]);
			}
			
			if (!checkForReal_mp(result, solve_options.T.real_threshold)) {
				std::cout << color::red();
				print_point_to_screen_matlab(result,"not_real_but_should_be_right");
                print_point_to_screen_matlab(Wright.pts_mp[kk],"nat");
				std::cout << color::console_default();
				ok = false;
			}
			
            clear_vec_mp(result);
            
            
			if (!ok) {
				std::cout << color::magenta();
				print_point_to_screen_matlab(midpoint_witness_sets[ii].pts_mp[kk],"generating_point");
				
				print_comp_matlab(&crit_downstairs->coord[ii],"left_proj_val");
				print_comp_matlab(&midpoints_downstairs->coord[ii],"midpt_proj_val");
				print_comp_matlab(&crit_downstairs->coord[ii+1],"right_proj_val");
				print_point_to_screen_matlab(particular_projection,"particular_projection");
				std::cout << color::console_default();
                
//                midpoint_witness_sets[ii].print_to_screen();
//                mypause();
                
                
                continue;
                
                
//                ii--;
                
                //
			}
			
			
			
			
			
			
			vec_cp_mp(temp_vertex.pt_mp,midpoint_witness_sets[ii].pts_mp[kk]);// set point
			temp_vertex.type = MIDPOINT; // set type
			
			temp_edge.midpt = index_in_vertices_with_add(V, temp_vertex); // gets the index of the new midpoint as it is added
			
			
			vec_cp_mp(temp_vertex.pt_mp,Wleft.pts_mp[kk]);// set point
			temp_vertex.type = NEW; // set type
			
			temp_edge.left  = index_in_vertices_with_add(V, temp_vertex);
			
			
			vec_cp_mp(temp_vertex.pt_mp,Wright.pts_mp[kk]);// set point
			temp_vertex.type = NEW; // set type
			
			temp_edge.right = index_in_vertices_with_add(V, temp_vertex);
			
			// keep track of those indices we found.
			
			found_indices_left.insert(temp_edge.left);
			found_indices_right.insert(temp_edge.right);
			
            
            found_indices_crit[ii].insert(temp_edge.left);
            found_indices_crit[ii+1].insert(temp_edge.right);
            found_indices_mid[ii].insert(temp_edge.midpt);

            
			int edge_num = add_edge(temp_edge);
			edge_occurence_tracker_left[temp_edge.left].push_back(edge_num);
			edge_occurence_tracker_right[temp_edge.right].push_back(edge_num);
			
			
			// count the number of times a critical point is tracked *to*.  this is for degenerate edge testing.  those which never get tracked to, make degenerate edges (isolated points are included in this).
			if (crit_point_counter.find(temp_edge.left)==crit_point_counter.end()) {
				crit_point_counter[temp_edge.left] = 1;
			}
			else{
				crit_point_counter[temp_edge.left] ++;
			}
			
			
			
			if (crit_point_counter.find(temp_edge.right)==crit_point_counter.end()) {
				crit_point_counter[temp_edge.right] = 1;
			}
			else{
				crit_point_counter[temp_edge.right] ++;
			}
			
			
			
			
			
			
			
			if (program_options.verbose_level>=2) {
				printf("done connecting upstairs midpoint %d (downstairs midpoint %d)\n",kk,ii);
                
				printf("indices of left, mid, right: %d %d %d\n",temp_edge.left,temp_edge.midpt,temp_edge.right);
				printf("\n\n");
			}
		}
		Wleft.reset();
		Wright.reset();
		
        
        
	}//re: for ii
	clear_mp(left_proj_val); clear_mp(right_proj_val);
	
	
    for (int ii=0; ii<num_midpoints; ii++) {
        V.assert_projection_value(found_indices_crit[ii], &crit_downstairs->coord[ii]);
        V.assert_projection_value(found_indices_mid[ii], &midpoints_downstairs->coord[ii]);
    }
	V.assert_projection_value(found_indices_crit[num_midpoints], &crit_downstairs->coord[num_midpoints]);
	
	
	
	
	
	// add degenerate edges for any critical point which is not mapped to by an edge.
	std::map<int,int>::iterator crit_pt_iterator;
	for (crit_pt_iterator = crit_point_counter.begin(); crit_pt_iterator != crit_point_counter.end(); crit_pt_iterator++) {
		int curr_index = crit_pt_iterator->first;
		int num_occurrences_local = crit_pt_iterator->second;
		
		if (num_occurrences_local==0) {
			
            //			vec_cp_mp(temp_vertex.pt_mp, V.vertices[curr_index].pt_mp);// set point
            //			projection_value_homogeneous_input(temp_vertex.projVal_mp,  V.vertices[curr_index].pt_mp,projections[0]);
            //			temp_vertex.type = ISOLATED; // set type
			
			edge E(curr_index,curr_index,curr_index);
			
			add_edge(E);
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	
	//              \\              //
	//\\\\\\\\\\\\\\\\\            /////////////
	///////////////////  merge    /////////////
	//////////////////            \\\\\\\\\\\\\\
	//             //              \\
	
	if (program_options.merge_edges==true) {
		mat_cp_mp(this->randomizer_matrix, randomizer_matrix);
		this->merge(midpoint_witness_sets[0],V,projections,solve_options);
	}// re: if merge_edges==true
	else
	{
		
		// since we are not merging, we need to NOT leave the type indicator as NEW, because it may throw off later merges.
		for (std::set<int>::iterator setiter = found_indices_right.begin(); setiter != found_indices_right.end(); setiter++) {
			int curr_index = *setiter;
			if (V.vertices[curr_index].type==NEW) { // only need to look at one of right and left here.
				V.vertices[curr_index].type = SEMICRITICAL;
			}
		}
		
		for (std::set<int>::iterator setiter = found_indices_left.begin(); setiter != found_indices_left.end(); setiter++) {
			int curr_index = *setiter;
			if (V.vertices[curr_index].type==NEW) { // only need to look at one of right and left here.
				V.vertices[curr_index].type = SEMICRITICAL;
			}
		}
		
	}//re: else merge_edges==false
	
	
	
	
	
	
	//done
	
	if (program_options.verbose_level>=0) {
		printf("num_edges = %d\n",num_edges);
	}
	
	
	
	
	
	
	
	return SUCCESSFUL;
} // re: interslice












//returns (-1,-1) if no candidate found
std::pair<int,int> curve_decomposition::get_merge_candidate(const vertex_set & V){
	
	std::pair< int, int> found_edges = std::pair<int,int>(-1,-1);
    
	
	// looking for edges with the type NEW, by looking at the left endpoint
	for (int tentative_right_edge=0; tentative_right_edge<this->num_edges; tentative_right_edge++) {
		if (V.vertices[edges[tentative_right_edge].left].type == NEW) {
			
			if (edges[tentative_right_edge].left == edges[tentative_right_edge].right)
				continue;
			
            
			int tentative_left_edge = this->edge_w_right(edges[tentative_right_edge].left);
			
			
			if (tentative_left_edge < 0) {
				std::cout << "found that edge " << tentative_right_edge << " has NEW leftpoint, but \\nexists edge w point " << edges[tentative_right_edge].left << " as right point." << std::endl;
				continue;
			}
			
			found_edges.second = tentative_right_edge;
			found_edges.first = tentative_left_edge;
			
			return found_edges;
		}
	}
	
	return found_edges;
}






void curve_decomposition::merge(witness_set & W_midpt,
                                vertex_set & V,
                                vec_mp * projections,
                                solver_configuration & solve_options)
{
	
	vec_mp particular_projection; init_vec_mp(particular_projection,0);
	vec_cp_mp(particular_projection, projections[0]);
	
	comp_mp half;  init_mp(half);  comp_mp temp;  init_mp(temp);  comp_mp temp2;  init_mp(temp2);
	mpf_set_str(half->r, "0.5", 10); mpf_set_str(half->i, "0.0", 10);
	
	
	multilin_config ml_config(solve_options);
	
	
	std::pair< int, int > edges_to_merge = this->get_merge_candidate(V);
	
	
	comp_mp new_proj_val;  init_mp2(new_proj_val,1024);
	
	
	while (edges_to_merge.first!=-1) { // this value is updated at the end of the loop
		// then there are edges superfluous and need to be merged
		int left_edge_w_pt = edges_to_merge.first;
		int right_edge_w_pt  = edges_to_merge.second;
		std::cout << color::cyan() << "merging edges " << left_edge_w_pt << " " << right_edge_w_pt << color::console_default() <<  std::endl;
		
		
		if (edges_to_merge.first < 0 || edges_to_merge.second < 0) {
			std::cout << "error: attemping to merge edges with negative index!" << std::endl;
			
			std::cout << "<" << edges[left_edge_w_pt].left << " " << edges[left_edge_w_pt].midpt << " " << edges[left_edge_w_pt].right << "> <";
			std::cout << edges[right_edge_w_pt].left << " " << edges[right_edge_w_pt].midpt << " " << edges[right_edge_w_pt].right << ">" << std::endl;
			
			
			br_exit(-8009);
		}
		
		
		
		
		
		
		
		std::cout << "<" << edges[left_edge_w_pt].left << " " << edges[left_edge_w_pt].midpt << " " << edges[left_edge_w_pt].right << "> <";
		std::cout << edges[right_edge_w_pt].left << " " << edges[right_edge_w_pt].midpt << " " << edges[right_edge_w_pt].right << ">" << std::endl;
		
		
		witness_set W_temp;
		
		
		// get the projection value of the midpoint we will be moving from.
		//arbitrarily chose to move from the midpoint of the left edge.
		projection_value_homogeneous_input(&particular_projection->coord[0], V.vertices[edges[left_edge_w_pt].midpt].pt_mp, projections[0]);
		neg_mp(&particular_projection->coord[0],&particular_projection->coord[0]);
		
		
		W_midpt.reset_linears(); // this witness set should only have a single linear
		W_midpt.add_linear(particular_projection);
		
		W_midpt.reset_points();
		W_midpt.add_point(V.vertices[edges[left_edge_w_pt].midpt].pt_mp);
		// I arbitrarily chose the left edge's midpoint as source to track to new midpoint.
		
		projection_value_homogeneous_input(temp,V.vertices[edges[left_edge_w_pt].left].pt_mp,projections[0]);
		projection_value_homogeneous_input(temp2,V.vertices[edges[right_edge_w_pt].right].pt_mp,projections[0]);
		
		vertex temp_vertex;
		add_mp(new_proj_val, temp, temp2);
		
		mul_mp(new_proj_val, new_proj_val, half); // now it is the average value
		
		neg_mp(&particular_projection->coord[0], new_proj_val); // set it in the linear for tracking
		
		
		solve_options.allow_multiplicity = 0;
		solve_options.allow_singular = 0;
		solve_options.complete_witness_set = 0;
		solve_options.robust = true;
		
		ml_config.set_randomizer(this->randomizer_matrix);
		
		multilin_solver_master_entry_point(W_midpt,         // witness_set
                                           &W_temp, // the new data is put here!
                                           &particular_projection,
                                           ml_config,
                                           solve_options);
		
        // each member of W_temp should real.  if a member of V already, mark index.  else, add to V, and mark.
		
		vec_cp_mp(temp_vertex.pt_mp, W_temp.pts_mp[0]);
		temp_vertex.type = MIDPOINT;
		
		edge temp_edge;
		temp_edge.left = edges[left_edge_w_pt].left;
		temp_edge.midpt = index_in_vertices_with_add(V, temp_vertex);
		temp_edge.right = edges[right_edge_w_pt].right;
		
		
		for (std::vector<int>::iterator vec_iter = edges[left_edge_w_pt].removed_points.begin(); vec_iter!=edges[left_edge_w_pt].removed_points.end(); vec_iter++)
		{
			temp_edge.removed_points.push_back( *vec_iter );
		}
		
		for (std::vector<int>::iterator vec_iter = edges[right_edge_w_pt].removed_points.begin(); vec_iter!=edges[right_edge_w_pt].removed_points.end(); vec_iter++)
		{
			temp_edge.removed_points.push_back( *vec_iter );
		}
		
		temp_edge.removed_points.push_back(edges[left_edge_w_pt].midpt);
		temp_edge.removed_points.push_back(edges[left_edge_w_pt].right);
		temp_edge.removed_points.push_back(edges[right_edge_w_pt].midpt);
		
		
		add_edge(temp_edge);
        //		std::cout << "adding edge " << temp_edge.left << " " << temp_edge.midpt << " " << temp_edge.right << " " << std::endl;
		//add the new_edge
		
		
		
		V.vertices[edges[left_edge_w_pt].midpt].removed = 1;
		V.vertices[edges[left_edge_w_pt].right].removed = 1;
		V.vertices[edges[right_edge_w_pt].midpt].removed = 1;
		
		// delete the two old edges // can't do this in the loop because we were using indexes into the edge set
		std::vector< edge > post_merge_edges;
		// this should be changed.
		for (int ii = 0; ii<this->num_edges; ii++) {
			if ((ii!=left_edge_w_pt) && (ii!=right_edge_w_pt)) { // if don't want to remove the edge
				post_merge_edges.push_back( this->edges[ii]);
			}
			else{
                //				std::cout << "removing edge " << ii << std::endl;
                //				std::cout << this->edges[ii].left << " " << this->edges[ii].midpt << " " << this->edges[ii].right << " " << std::endl;
			}
			//otherwise skip it.
		}
		
		
		
		//swap this's edge info to post-merge info.
		this->edges.swap(post_merge_edges);
		this->num_edges = int(this->edges.size());
		
		edges_to_merge = this->get_merge_candidate(V);
	}// re: while
	clear_mp(half); clear_mp(temp); clear_mp(temp2);
	
	
	clear_mp(new_proj_val);
	
	
} // re: merge















//subfunctions
int curve_decomposition::compute_critical_points(const witness_set & W_curve,
                                                 mat_mp randomizer_matrix,
                                                 std::vector<int> randomized_degrees,
                                                 vec_mp *projections,
                                                 BR_configuration & program_options,
                                                 solver_configuration & solve_options,
                                                 witness_set & W_crit_real)
{
    
	int *declarations = NULL;
	
	W_crit_real.input_filename = W_curve.input_filename;
	

	
    
    
	
	nullspace_config ns_config;
	compute_crit_nullspace(&W_crit_real, // the returned value
                           W_curve,            // input the original witness set
                           randomizer_matrix,
                           projections,
                           randomized_degrees,//this is computed
                           1,  // dimension of ambient complex object
                           1,   //  target dimension to find
                           1,   // COdimension of the critical set to find.
                           program_options,
                           solve_options,
                           &ns_config);
	ns_config.clear();
	
    
	
	W_crit_real.cp_patches(W_curve);
	W_crit_real.cp_linears(W_curve);
	
	W_crit_real.only_first_vars(W_curve.num_variables); // trim the fat, since we are at the lowest level.
	W_crit_real.sort_for_real(solve_options.T);
	
	W_crit_real.print_to_screen();
	
	if (have_sphere_radius) {
		W_crit_real.sort_for_inside_sphere(sphere_radius, sphere_center);
		
//		print_point_to_screen_matlab(sphere_center,"center");
//		print_comp_matlab(sphere_radius,"radius");
	}
	else
	{
//		// this should only happen if we are decomposing a curve at the highest level.  if anything higher dimension, this curve is a subcomponent, and we should NOT compute bounds here.
		std::cout << color::red() << "computing sphere bounds..." << color::console_default() << std::endl;
		compute_sphere_bounds(W_crit_real);
	}
	
	
    witness_set W_additional;
	// now get the sphere intersection critical points and ends of the interval
	get_additional_critpts(&W_additional,  // the returned value
                           W_curve,       // all else here is input
                           randomizer_matrix,  //
                           program_options,
                           solve_options);
	
    W_additional.sort_for_real(solve_options.T);
	W_additional.sort_for_unique(solve_options.T);
	

    
    W_crit_real.merge(W_additional);
	
	
	return SUCCESSFUL;
}




int curve_decomposition::get_additional_critpts(witness_set *W_additional,
                                                const witness_set & W_curve,
                                                mat_mp randomizer_matrix,
                                                BR_configuration & program_options,
                                                solver_configuration & solve_options)
{
	
    
    if (W_curve.num_linears!=1) {
        std::cout << color::red() << "the input witness set to get_additional_critpts had an incorrect number of linears: " << W_curve.num_linears << color::console_default() << std::endl;
        br_exit(-518);
    }
    
    
	//build up the start system
	solve_options.robust = true;
	
	
	int blabla;
	
	
	
	parse_input_file(W_curve.input_filename, &blabla);
	preproc_data_clear(&solve_options.PPD); // ugh this sucks
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	std::vector< int > local_degrees;
	make_randomization_matrix_based_on_degrees(randomizer_matrix, local_degrees,
                                               W_curve.num_variables-W_curve.num_patches-W_curve.num_linears,
                                               solve_options.PPD.num_funcs);
	
	multilin_config ml_config(solve_options); // copies in the randomizer matrix and sets up the SLP & globals.
	
	ml_config.set_randomizer(randomizer_matrix);
	
	vec_mp *multilin_linears = (vec_mp *) br_malloc(1*sizeof(vec_mp));
	init_vec_mp2(multilin_linears[0],W_curve.num_variables,solve_options.T.AMP_max_prec);
	multilin_linears[0]->size = W_curve.num_variables;
    
	
	witness_set W_sphere = W_curve;
    //grab just the shell of the input witness set
    W_sphere.reset_points();
    W_sphere.reset_linears();
    W_sphere.reset_patches();
    
    
	sphere_config sp_config(randomizer_matrix);
	for (int jj=0; jj<W_curve.num_variables; jj++) {
        set_zero_mp(&multilin_linears[0]->coord[jj]);
    }
	
	for (int ii=0; ii<2; ii++) {
		
		for (int jj=0; jj<W_curve.num_natural_vars(); jj++) {
			get_comp_rand_mp(&multilin_linears[0]->coord[jj]);
		}
		
		vec_cp_mp(sp_config.starting_linear[ii], multilin_linears[0]);
		
		witness_set W_temp;
		
		
		
		solve_options.allow_singular = 0;
		solve_options.complete_witness_set = 0;
		solve_options.allow_multiplicity = 0;
		solve_options.allow_unsuccess = 0;
		
		
		
		multilin_solver_master_entry_point(W_curve,         // witness_set
                                           &W_temp, // the new data is put here!
                                           multilin_linears,
                                           ml_config,
                                           solve_options);
		
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
	
	
	
	
	solve_options.allow_singular = 1;
	solve_options.complete_witness_set = 1;
	solve_options.allow_multiplicity = 1;
	solve_options.allow_unsuccess = 0;
	
	
	
	

	sphere_solver_master_entry_point(W_sphere,
                                     W_additional, // returned value
                                     sp_config,
                                     solve_options);
	
	

	return 0;
}

















int get_sum_degrees(char filename[], int num_funcs){
	int degsum = 0, tmpdeg, ii;
	
	FILE *IN;
	
	IN =  safe_fopen_read(filename);
	
	for (ii = 0; ii<num_funcs; ii++) {
		fscanf(IN,"%d",&tmpdeg);
		degsum += tmpdeg;
	}
	
	fclose(IN);
	
	
	return degsum;
}






// will compute a randomizer matrix since you don't provide one. must have current PPD in solve_options for this to work correctly
int verify_projection_ok(const witness_set & W,
                         vec_mp * projection,
                         solver_configuration & solve_options)
{
	
	
	//create a matrix
	mat_mp randomizer_matrix;
	init_mat_mp(randomizer_matrix,W.num_variables-W.num_patches-W.dim,solve_options.PPD.num_funcs); // <--- the PPD had better be current here
	
	//create the array of integers
	std::vector<int> randomized_degrees;
	
	//get the matrix and the degrees of the resulting randomized functions.
	make_randomization_matrix_based_on_degrees(randomizer_matrix, randomized_degrees, W.num_variables-W.num_patches-W.dim, solve_options.PPD.num_funcs);
	
	int invalid_flag = verify_projection_ok(W, randomizer_matrix, projection, solve_options);
	
	clear_mat_mp(randomizer_matrix);
	
	return invalid_flag;
}








int verify_projection_ok(const witness_set & W,
                         mat_mp randomizer_matrix,
                         vec_mp * projection,
                         solver_configuration & solve_options)
{
	int ii,jj;
	
	
	int invalid_flag;
	
	
	vec_mp temp_rand_point;  init_vec_mp(temp_rand_point,W.num_variables); temp_rand_point->size = W.num_variables;
	set_one_mp(&temp_rand_point->coord[0]); // first coordinate must be 1
	for (ii=1; ii<W.num_variables; ++ii) {
		get_comp_rand_mp(&temp_rand_point->coord[ii]);
	}
	
	prog_t SLP;
	setupProg(&SLP, solve_options.T.Precision, 2);
	
	
	comp_mp zerotime; init_mp(zerotime);
	set_zero_mp(zerotime);
	
	
	
	eval_struct_mp ED; init_eval_struct_mp(ED, 0, 0, 0);
	evalProg_mp(ED.funcVals, ED.parVals, ED.parDer, ED.Jv, ED.Jp, temp_rand_point, zerotime, &SLP);
	
	mat_mp AtimesJ; init_mat_mp(AtimesJ, 1, 1); AtimesJ->rows = AtimesJ->cols = 1;
	
	mat_mul_mp(AtimesJ, randomizer_matrix, ED.Jv);
	
	
	mat_mp detme;  init_mat_mp(detme, W.num_variables-1, W.num_variables-1);
	detme->cols = detme->rows= W.num_variables-1;
	
	
	//inflate the matrix
	int dim = W.num_variables - randomizer_matrix->rows - 1;
	
	for (ii=0; ii< AtimesJ->rows; ++ii) {
		for (jj=0; jj<AtimesJ->cols - 1; ++jj) {
			set_mp(&detme->entry[ii][jj],&AtimesJ->entry[ii][jj+1]); // omit the homogeneous coordinate's columns
		}
	}
	
	int offset = W.num_variables - 1 - dim;
	for (jj=0; jj < dim; jj++){
		for (ii=0; ii<W.num_variables-1; ii++) {
			set_mp(&detme->entry[offset+jj][ii],&projection[jj]->coord[ii+1]);
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
	
	return invalid_flag;
	
}










/***************************************************************\
 * USAGE: compute the isolated points for Non-Self-Conjugate case            *
 * ARGUMENTS: witness set, projection, # of variables, name of input file*
 * RETURN VALUES: Curve decomposition*
 * NOTES:                                                        *
 \***************************************************************/
void curve_decomposition::computeCurveNotSelfConj(const witness_set & W_in,
                                                  vec_mp         pi_mp,
                                                  vertex_set			&V,
                                                  int           num_vars,
                                                  BR_configuration & program_options,
                                                  solver_configuration & solve_options)

{
	
	// num_vars includes the homogeneous variable
	
	
    int ii,jj,num_sols,*declarations = NULL;
    
	std::string bertini_system_command = program_options.bertini_command;
	bertini_system_command.append(" input_NSC");
	
    FILE *IN = NULL;
    vec_mp cur_sol,cur_sol_bar;
    init_vec_mp(cur_sol,num_vars); cur_sol->size = num_vars;
	set_one_mp(&cur_sol->coord[0]);
	
    init_vec_mp(cur_sol_bar,num_vars); cur_sol_bar->size = num_vars;
	set_one_mp(&cur_sol_bar->coord[0]);
	
    partition_parse(&declarations, W_in.input_filename, "func_input_nsc", "config_nsc",1);
	
	
    //generate input file
    diag_homotopy_input_file("input_NSC", "func_input_nsc","func_inputbar","config_nsc",W_in.L_mp[0],num_vars-1);
    //generate start file
	diag_homotopy_start_file("start",  W_in);
	
	
    //run bertini
	
	copyfile("witness_data","witness_data_0");
	
    system(bertini_system_command.c_str());
	
	rename("witness_data_0","witness_data");
	
	
    //read the real solutions
    IN = safe_fopen_read("real_solutions");
    
    fscanf(IN, "%d\n\n", &num_sols);
	
	
	vertex temp_vertex;
	//	init_vertex(&temp_vertex);
	change_size_vec_mp(temp_vertex.pt_mp,num_vars); temp_vertex.pt_mp->size = num_vars;
	
	temp_vertex.type = ISOLATED;
	
	comp_mp projection_value;  init_mp(projection_value);
	
	for(ii=0;ii<num_sols;ii++) {
		for(jj=0;jj<num_vars-1;jj++){
			mpf_inp_str(cur_sol->coord[jj+1].r, IN, 10);
			mpf_inp_str(cur_sol->coord[jj+1].i, IN, 10);
			
			mpf_inp_str(cur_sol_bar->coord[jj+1].r, IN, 10);
			mpf_inp_str(cur_sol_bar->coord[jj+1].i, IN, 10);
		}
        
        //check if x=x_bar
		
        if (isSamePoint_homogeneous_input(cur_sol,cur_sol_bar)) { // x=x_bar
			
			vec_cp_mp(temp_vertex.pt_mp,cur_sol);
			
			dot_product_mp(projection_value, temp_vertex.pt_mp, pi_mp);// i think this is wrong, because pi_mp probably has a 0 at the front of it.
			
			index_in_vertices_with_add(V, temp_vertex);
			
		}
	}
	
	clear_mp(projection_value);
    //set the number of vertices
    fclose(IN);
    //clear
    clear_vec_d(cur_sol);
    clear_vec_d(cur_sol_bar);
    free(declarations);
    
	// delete temporary files
    remove("func_input_nsc");
    remove("config_nsc");
    remove("func_inputbar");
    remove("var_names");
	
	
}








int curve_decomposition::setup_edges(boost::filesystem::path INfile)
//setup the vertex structure
{
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
	
    //	std::cout << "printing curve decomposition to folder " << base << std::endl;
	
	decomposition::print(base);
	
	boost::filesystem::path edgefile = base / "E.edge";
	
	curve_decomposition::print_edges(edgefile);
}






/**Output edge structure as follows:
 //this needs to be updated...
 # variables
 # edges
 name of input file
 edge 1
 
 edge 2
 
 
 for each edge, output the following information:
 index to left vertex in vertices
 index to right vertex in vertices
 index to midpoint vertex in vertices
 
 **/
void curve_decomposition::print_edges(boost::filesystem::path outputfile)
{
	int ii;
	FILE *OUT = safe_fopen_write(outputfile);
	
	// output the number of vertices
	fprintf(OUT,"%d\n\n",num_edges);
	
	for(ii=0;ii<num_edges;ii++)
		fprintf(OUT,"%d %d %d \n",
                edges[ii].left,
                edges[ii].midpt,
                edges[ii].right);
	fclose(OUT);
}











// MISC. FUNCTIONS


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
    
	FILE *OUT = safe_fopen_write(outputFile.c_str());
	
    // setup configurations in OUT
    fprintf(OUT, "CONFIG\n");
    IN = safe_fopen_read(configInput.c_str());
    fclose(IN);
	
    fprintf(OUT, "USERHOMOTOPY: 1;\nDeleteTempFiles: 1;\nEND;\nINPUT\n");
	
    // setup variables in OUT
    IN = safe_fopen_read(funcInputx.c_str());
    while ((ch = fgetc(IN)) != EOF )
        fprintf(OUT, "%c", ch);
    fclose(IN);
	
	//setup the function name in OUT
	IN = safe_fopen_read(funcInputy.c_str());
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
	
	comp_mp temp;  init_mp(temp);
    for (ii = 0; ii < L->size; ii++)
    {
        fprintf(OUT, "L%d = ",ii);
        // print output
		print_mp(OUT, 0, &L->coord[ii]); fprintf(OUT, "\n");
		//    fprintf(OUT, fmt, L->coord[ii].r, L->coord[ii].i);
        fprintf(OUT, ";\n");
		
        fprintf(OUT, "Lbar%d = ",ii);
		conjugate_mp(temp, &L->coord[ii]);
        // print output
		print_mp(OUT, 0, temp); fprintf(OUT, "\n");
		//    fprintf(OUT, fmt, L->coord[ii].r, -L->coord[ii].i);
        fprintf(OUT, ";\n");
    }
    fprintf(OUT, "\n");
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
    fprintf(OUT, "\nL=t*(");
    //(Lx-1)*t+(1-t)*A[0]*(x-x_bar)
    for(ii=0;ii<num_vars;ii++)
    {
        fprintf(OUT, "L%d*", ii);
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
    fprintf(OUT, "0);\nLbar=t*(");
    for(ii=0;ii<num_vars;ii++)
    {
        fprintf(OUT, "Lbar%d*", ii);
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
    
    char *fmt = NULL;
    int ii,jj,kk;
	int size,digits=15;
	
	
	FILE *OUT = safe_fopen_write(startFile.c_str());
    
    size = 1 + snprintf(NULL, 0, "%%.%de %%.%de\n", digits, digits);
    // allocate size
    fmt = (char *)br_malloc(size * sizeof(char));
    // setup fmt & fmtb
    sprintf(fmt, "%%.%de %%.%de\n", digits, digits);
	// output the number of start points
    fprintf(OUT,"%d\n\n",W.num_pts*W.num_pts);
	
	
	comp_mp temp; init_mp(temp);
	
	
	for (kk=0; kk<W.num_pts; kk++){
		for (ii=0; ii<W.num_pts; ii++) { // output {w \bar{w}}'
			
			vec_mp result; init_vec_mp(result,0);
			vec_mp result2; init_vec_mp(result2,0);
			
			dehomogenize(&result,W.pts_mp[ii]);
			dehomogenize(&result2,W.pts_mp[kk]);
			
			
			for(jj=0; jj<W.num_variables-1;jj++) {
				print_mp(OUT, 15, &result->coord[jj]); fprintf(OUT, "\n");
				
				conjugate_mp(temp, &result2->coord[jj] )
				print_mp(OUT, 15, temp); fprintf(OUT, "\n");
				//      fprintf(OUT, fmt, result->coord[jj].r,  result->coord[jj].i);
			}
			
			
			
			fprintf(OUT,"\n");
		}// re: ii
	}// re: kk
	
    free(fmt);
    fclose(OUT);
}































