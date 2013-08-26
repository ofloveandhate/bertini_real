#include "surface.hpp"











void surface_decomposition::print(boost::filesystem::path base)
{
	
	std::cout << "printing surface decomposition to folder " << base << std::endl;
	decomposition::print(base);
	
	
	
	boost::filesystem::path summaryname = base;
	summaryname /= "S.surf";
	FILE *OUT = safe_fopen_write(summaryname);
	fprintf(OUT,"%d %d %ld %ld\n\n", num_faces, num_edges, mid_slices.size(), mid_slices.size());
	// what more to print here?
	fclose(OUT);
	
	summaryname = base;
	summaryname /= "F.faces";
	print_faces(summaryname);
	
	
	boost::filesystem::path curve_location = base;
	curve_location /= "curve";
	
	for (int ii=0; ii<mid_slices.size(); ii++) {
		
		std::stringstream converter;
		converter << ii;
		
		boost::filesystem::path specific_loc = curve_location;
		specific_loc += "_midslice_";
		specific_loc += converter.str();
		converter.clear(); converter.str("");
		
		mid_slices[ii].print(specific_loc);
	}
	
	for (int ii=0; ii<crit_slices.size(); ii++) {
		
		std::stringstream converter;
		converter << ii;
		
		boost::filesystem::path specific_loc = curve_location;
		specific_loc += "_critslice_";
		specific_loc += converter.str();
		converter.clear(); converter.str("");
		
		crit_slices[ii].print(specific_loc);
	}
	
	boost::filesystem::path specific_loc = curve_location;
	specific_loc += "_crit";
	crit_curve.print(specific_loc);
}












void surface_decomposition::main(vertex_set & V,
																 witness_set & W_surf,
																 vec_mp *pi,
																 BR_configuration & program_options,
																 solver_configuration & solve_options)
{

#ifdef usetempfolders
	program_options.move_to_temp();
#endif
	
	input_filename = W_surf.input_filename;
	int ambient_dim = 2;
	
	num_variables = W_surf.num_variables;
	component_num = W_surf.comp_num;
	add_projection(pi[0]);
	add_projection(pi[1]);
	
	for (int ii=0; ii<W_surf.num_patches; ii++)
		add_patch(W_surf.patch_mp[ii]);
	
	
	// perform an isosingular deflation
	// note: you do not need witness_data to perform isosingular deflation
	if (program_options.verbose_level>=2)
		printf("performing isosingular deflation\n");
	
	W_surf.write_dehomogenized_coordinates("witness_points_dehomogenized"); // write the points to file
	int num_deflations, *deflation_sequence = NULL;
	
	isosingular_deflation(&num_deflations, &deflation_sequence, program_options,
																 program_options.current_working_filename,
																 "witness_points_dehomogenized",
																 program_options.max_deflations, W_surf.dim, W_surf.comp_num);
	free(deflation_sequence);
	
	
	program_options.input_deflated_filename = program_options.current_working_filename;
	
	std::stringstream converter;
	converter << "_dim_" << W_surf.dim << "_comp_" << W_surf.comp_num << "_deflated";
	program_options.input_deflated_filename += converter.str();
	converter.clear(); converter.str("");
	
	// this wraps around a bertini routine
	parse_input_file(program_options.input_deflated_filename);
	
	preproc_data_clear(&solve_options.PPD);
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	
	if (program_options.verbose_level>=2) {
		printf("checking if component is self-conjugate\n");
	}
	
	
	checkSelfConjugate(W_surf,num_variables,program_options, program_options.current_working_filename);  //later:  could be passed in from user, if we want
	
	//regenerate the various files, since we ran bertini since then.
	parse_input_file(program_options.input_deflated_filename);
	
	
	
	


	//create the matrix
	init_mat_mp2(this->randomizer_matrix,W_surf.num_variables-1-ambient_dim,solve_options.PPD.num_funcs,solve_options.T.AMP_max_prec);
	
	//create the array of integers
	int *randomized_degrees = (int *)br_malloc((W_surf.num_variables-1-ambient_dim)*sizeof(int));
	
	//get the matrix and the degrees of the resulting randomized functions.
	make_randomization_matrix_based_on_degrees(this->randomizer_matrix, &randomized_degrees, W_surf.num_variables-1-ambient_dim, solve_options.PPD.num_funcs);
	
	
	

	
	// find witness points on the critical curve.
	std::cout << "computing witness points for the critical curve" << std::endl;
	
	witness_set W_curve_crit;
	nullspace_config ns_config;
	compute_crit_nullspace(&W_curve_crit, // the returned value
												 W_surf,            // input the original witness set
												 randomizer_matrix,
												 pi,
												 randomized_degrees,
												 ambient_dim,  // dimension of ambient complex object
												 2,   //  target dimension to find
												 2,   // COdimension of the critical set to find.
												 program_options,
												 solve_options,
												 &ns_config);
	
	W_curve_crit.only_first_vars(num_variables);
	W_curve_crit.sort_for_real(solve_options.T);
	W_curve_crit.write_dehomogenized_coordinates("W_curve_crit"); // write the points to file
	
	ns_config.clear();
	
	
	
	
	
	// compute witness points for the critical curve.  
	witness_set W_curve;
	
	std::cout << "computing critical points of the critical curve" << std::endl;
	compute_crit_nullspace(&W_curve, // the returned value
												 W_surf,            // input the original witness set
												 randomizer_matrix,
												 pi,
												 randomized_degrees,
												 ambient_dim,  // dimension of ambient complex object
												 2,   //  target dimension to find
												 1,   // COdimension of the critical set to find.
												 program_options,
												 solve_options,
												 &ns_config);
	
//	W_curve.only_first_vars(num_vars);
	W_curve.write_dehomogenized_coordinates("W_crit"); // write the points to file

	W_curve.input_filename = "input_critical_curve";
	
	// this system describes the system for the critical curve
	create_nullspace_system("input_critical_curve",
													boost::filesystem::path(program_options.called_dir) / program_options.input_filename,
													program_options, &ns_config);
	ns_config.clear();
	
	// this should be deprecated
	program_options.current_working_filename = boost::filesystem::absolute("input_critical_curve");
	
	

	
	vec_mp temp_proj; init_vec_mp2(temp_proj,0,1024);
	vec_cp_mp(temp_proj, pi[0]);
	std::cout << temp_proj->size << std::endl;
	increase_size_vec_mp(temp_proj, W_curve.num_variables); temp_proj->size = W_curve.num_variables;
	for (int ii=W_surf.num_variables; ii<W_curve.num_variables; ii++) {
		set_zero_mp(&temp_proj->coord[ii]);
	}
	
	for (int ii=0; ii<W_curve.num_linears; ii++) {
		increase_size_vec_mp(W_curve.L_mp[ii], W_curve.num_variables); W_curve.L_mp[ii]->size = W_curve.num_variables;
		
		for (int jj=W_surf.num_variables; jj<W_curve.num_variables; jj++) {
			set_zero_mp(&W_curve.L_mp[ii]->coord[jj]);
		}
	}
	
	
	crit_curve.interslice(W_curve,
								 W_curve_crit,
								 randomizer_matrix,
								 &temp_proj,
								 program_options,
								 solve_options,
								 V);
	
	clear_vec_mp(temp_proj);
	
	crit_curve.add_projection(pi[0]);
	
	crit_curve.num_variables = num_variables + num_variables-1;
	
	// DONE COMPUTING THE CRITICAL CURVE NOW.
	std::cout << "done decomposing critical curve" << std::endl;
	
	
	
	vec_mp *multilin_linears = (vec_mp *) br_malloc(2*sizeof(vec_mp));
	
	
	init_vec_mp2(multilin_linears[0], W_surf.num_variables,1024);
	init_vec_mp2(multilin_linears[1], W_surf.num_variables,1024);
	vec_cp_mp(multilin_linears[0],pi[0]);
	vec_cp_mp(multilin_linears[1],W_surf.L_mp[1]);
	
	vec_mp midpoints_downstairs, crit_downstairs;
	init_vec_mp(midpoints_downstairs,0); init_vec_mp(crit_downstairs,0);
	std::vector< int > index_tracker;
	
	W_curve_crit.compute_downstairs_crit_midpts(crit_downstairs, midpoints_downstairs, index_tracker, pi[0]);
	
	print_point_to_screen_matlab(crit_downstairs,"crit_down");
	print_point_to_screen_matlab(midpoints_downstairs,"middown");
	std::vector< witness_set > midpoint_witness_sets;
	midpoint_witness_sets.resize(midpoints_downstairs->size);
	
	mid_slices.resize(midpoints_downstairs->size);
	
	

	
	
	for (int ii=0; ii<midpoints_downstairs->size; ii++){
		std::cout << "decomposing the " << ii << "th midpoint slice" << std::endl;
		std::stringstream converter;
		converter << ii;
		
		int blabla, *declarations;
		parse_input_file(W_surf.input_filename, &blabla);
		partition_parse(&declarations, W_surf.input_filename, "func_input", "config", 0); // the 0 means not self conjugate.
		free(declarations);																													 // i would like to move this.
		preproc_data_clear(&solve_options.PPD);
		parse_preproc_data("preproc_data", &solve_options.PPD);
		
		neg_mp(&multilin_linears[0]->coord[0], &midpoints_downstairs->coord[ii]);
		
		multilintolin_solver_main(solve_options.T.MPType,
															W_surf,         // witness_set
															randomizer_matrix,
															multilin_linears, //  the set of linears we will solve at.
															&midpoint_witness_sets[ii], // the new data is put here!
															solve_options); // already a pointer
		
		boost::filesystem::path slicename = W_surf.input_filename;
		slicename += "_midslice_";
		slicename += converter.str();
		create_sliced_system(W_surf.input_filename, slicename, &multilin_linears[0], 1, W_surf);
		
		
		parse_input_file(slicename, &blabla);
		partition_parse(&declarations, slicename, "func_input", "config", 0); // the 0 means not self conjugate.
		free(declarations);																													 // i would like to move this.
		preproc_data_clear(&solve_options.PPD);
		parse_preproc_data("preproc_data", &solve_options.PPD);
		
		
		midpoint_witness_sets[ii].input_filename = slicename;
		midpoint_witness_sets[ii].add_linear(multilin_linears[1]);
		midpoint_witness_sets[ii].add_patch(W_surf.patch_mp[0]);
		midpoint_witness_sets[ii].variable_names = W_surf.variable_names;
		midpoint_witness_sets[ii].dim = 1;
		
		// we already know the component is self-conjugate (by entry condition), so we are free to call this function
		mid_slices[ii].computeCurveSelfConj(midpoint_witness_sets[ii],
												 &pi[1],
												 V,
												 num_variables,
												 program_options, solve_options);
		
		mid_slices[ii].add_projection(pi[1]);
		mid_slices[ii].num_variables = num_variables;
		
		std::cout << "DONE decomposing the " << ii << "th midpoint slice" << std::endl;

	}

	
	
	
	std::vector< witness_set > critpoint_witness_sets;
	critpoint_witness_sets.resize(crit_downstairs->size);
	
	crit_slices.resize(crit_downstairs->size);
	
	
	for (int ii=0; ii<crit_downstairs->size; ii++){
		std::cout << "decomposing the " << ii << "th critpoint slice" << std::endl;
		
		std::stringstream converter;
		converter << ii;
		
		int blabla, *declarations;
		parse_input_file(W_surf.input_filename, &blabla);
		partition_parse(&declarations, W_surf.input_filename, "func_input", "config", 0); // the 0 means not self conjugate.
		free(declarations);																													 // i would like to move this.
		preproc_data_clear(&solve_options.PPD);
		parse_preproc_data("preproc_data", &solve_options.PPD);
		
		neg_mp(&multilin_linears[0]->coord[0], &crit_downstairs->coord[ii]);
		
		
		
		multilintolin_solver_main(solve_options.T.MPType,
														W_surf,         // witness_set
														randomizer_matrix,
														multilin_linears, //  the set of linears we will solve at.
														&critpoint_witness_sets[ii], // the new data is put here!
														solve_options); // already a pointer
		
		boost::filesystem::path slicename = W_surf.input_filename;
		slicename += "_critslice_";
		slicename += converter.str();
		create_sliced_system(W_surf.input_filename, slicename, &multilin_linears[0], 1, W_surf);
		
		
		parse_input_file(slicename, &blabla);
		partition_parse(&declarations, slicename, "func_input", "config", 0); // the 0 means not self conjugate.
		free(declarations);																													 // i would like to move this.
		preproc_data_clear(&solve_options.PPD);
		parse_preproc_data("preproc_data", &solve_options.PPD);
			
		critpoint_witness_sets[ii].input_filename = slicename;
		critpoint_witness_sets[ii].add_linear(multilin_linears[1]);
		critpoint_witness_sets[ii].add_patch(W_surf.patch_mp[0]);
		critpoint_witness_sets[ii].variable_names = W_surf.variable_names;
		critpoint_witness_sets[ii].dim = 1;
		
		solve_options.backup_tracker_config();

		
		while (crit_slices[ii].num_edges==0) {
			crit_slices[ii].clear();
			
			print_comp_matlab(&crit_downstairs->coord[ii], "target_proj");
			
//			std::cout << "critpoint_witness_sets[ii]" << std::endl;
//			critpoint_witness_sets[ii].print_to_screen();
			// we already know the component is self-conjugate (by entry condition), so we are free to call this function
			crit_slices[ii].computeCurveSelfConj(critpoint_witness_sets[ii],
																					 &pi[1], V, num_variables,
																					 program_options, solve_options);
			
			crit_slices[ii].add_projection(pi[1]);
			crit_slices[ii].num_variables = num_variables;
			solve_options.T.final_tolerance = 0.1*solve_options.T.final_tolerance;
			solve_options.T.securityMaxNorm = 10*solve_options.T.securityMaxNorm;
			
			if (crit_slices[ii].num_edges==0) {
				std::cout << "crit_slice[" << ii << "].num_edges == 0" << std::endl;
				mypause();
			}
		}
		solve_options.reset_tracker_config();
		std::cout << "done decomposing the " << ii << "th critpoint slice" << std::endl;
		
	}
	
	
	//connect the dots
	connect_the_dots(V,
									 pi,
									 program_options,
									 solve_options);
	
	
	

#ifdef usetempfolders
	program_options.move_to_called();
#endif
	
	
	
}







void surface_decomposition::connect_the_dots(vertex_set & V,
																						 vec_mp *pi,
																						 BR_configuration & program_options,
																						 solver_configuration & solve_options)
{
	
	midpoint_config md_config;
	
	
	md_config.setup(*this, solve_options);
	
	
	comp_mp temp, temp2, temp3;
	init_mp(temp);
	init_mp(temp2);
	init_mp(temp3);
	
	comp_mp numer, denom;
	init_mp(numer);
	init_mp(denom);
	
	comp_mp proj_top, proj_bottom, proj_mid;
	init_mp(proj_mid);
	init_mp(proj_bottom);
	init_mp(proj_top);
	
	int num_crit_vars = V.vertices[crit_curve.edges[0].right].pt_mp->size; // is this guaranteed to exist?
	
	witness_set W_midtrack;
	W_midtrack.num_variables = 2*num_crit_vars + this->num_variables;
	
	vec_mp blank_point;  init_vec_mp2(blank_point, this->num_variables + 2*num_crit_vars,1024);
	blank_point->size = this->num_variables + 2*num_crit_vars;
	
	
	W_midtrack.add_point(blank_point);
	
	
	
	for (int ii=0; ii<this->num_patches; ii++) {
		W_midtrack.add_patch(this->patch[ii]);
	}
	
	for (int mm=0;mm<2;mm++){
		for (int ii=0; ii<crit_curve.num_patches; ii++) {
			W_midtrack.add_patch(crit_curve.patch[ii]);
		}
	}
	
	
	for (int ii=0; ii<mid_slices.size(); ii++) {
		for (int jj=0; jj<mid_slices[ii].num_edges; jj++) {
			std::cout << "face " << num_faces << ", edge " << jj << std::endl;
			face F;
			
			//create the face here
			
			F.index = ii; // the index of which midslice this face came from.
			F.midpt = mid_slices[ii].edges[jj].midpt; // index the point
			
			// perform a search to find the top and bottom edges in crit_curve
			F.top = crit_curve.edge_w_midpt(mid_slices[ii].edges[jj].right); // index the *edge*
			F.bottom = crit_curve.edge_w_midpt(mid_slices[ii].edges[jj].left); // index the *edge*
			
			
			
			
			int var_counter = 0;
			for (int kk=0; kk<this->num_variables; kk++) {
				set_mp(&W_midtrack.pts_mp[0]->coord[kk], &V.vertices[mid_slices[ii].edges[jj].midpt].pt_mp->coord[kk]);
				var_counter++;
			}
			
			int offset = var_counter;
			for (int kk=0; kk<num_crit_vars; kk++) {
				set_mp(&W_midtrack.pts_mp[0]->coord[kk+offset], &V.vertices[mid_slices[ii].edges[jj].left].pt_mp->coord[kk]); // y0
				var_counter++;
			}
			
			offset = var_counter;
			for (int kk=0; kk<num_crit_vars; kk++) {
				set_mp(&W_midtrack.pts_mp[0]->coord[kk+offset], &V.vertices[mid_slices[ii].edges[jj].right].pt_mp->coord[kk]); // y2
				var_counter++;
			}
			
			

			
			// the v direction corresponds to pi[1].
			
			
			
			// make u, v target values.
			
			//track
			projection_value_homogeneous_input(md_config.crit_val_left, V.vertices[ crit_slices[ii].edges[0].midpt ].pt_mp,pi[0]);
			projection_value_homogeneous_input(md_config.crit_val_right, V.vertices[crit_slices[ii+1].edges[0].midpt].pt_mp,pi[0]);
			
			// the u direction corresponds to pi[0].
			for (int zz=0; zz<2; zz++) { // go left and right (or up and down, whatever you want, just choose two directions)
				std::cout << zz << "\n";
				if (zz==0){
					set_zero_mp(md_config.u_target);}
				else{
					set_one_mp(md_config.u_target);}
				
				if (crit_slices[ii+zz].is_degenerate()) {
					// can simply set the top or bottom edge to be this one.  know it goes there.  only have to 
					std::cout << "crit_slice[" << ii+zz << "] is degenerate" << std::endl;
					if (zz==0){
						F.left.push_back(0);
						F.num_left++;
					}
					else{
						F.right.push_back(0);
						F.num_right++;
					}
					
				}
				else
				{

					
					if (zz==0) {
						projection_value_homogeneous_input(proj_top,    V.vertices[ crit_curve.edges[F.top].left ].pt_mp,   pi[1]); //w2
						projection_value_homogeneous_input(proj_bottom, V.vertices[ crit_curve.edges[F.bottom].left ].pt_mp,pi[1]); //w0
						
						print_point_to_screen_matlab(V.vertices[ crit_curve.edges[F.top].left ].pt_mp,"top_target");
						print_point_to_screen_matlab(V.vertices[ crit_curve.edges[F.bottom].left ].pt_mp,"bottom_target");
						
						
					}
					else{ // zz==1, and going right
						projection_value_homogeneous_input(proj_top, V.vertices[ crit_curve.edges[F.top].right ].pt_mp,pi[1]);
						projection_value_homogeneous_input(proj_bottom, V.vertices[ crit_curve.edges[F.bottom].right ].pt_mp,pi[1]);
						
						print_point_to_screen_matlab(V.vertices[ crit_curve.edges[F.top].right ].pt_mp,"top_target");
						print_point_to_screen_matlab(V.vertices[ crit_curve.edges[F.bottom].right ].pt_mp,"bottom_target");
						
					}
					
					

					for (int yy=0; yy<crit_slices[ii+zz].num_edges; yy++)
					{
						
						std::cout << "face " << ii << " zz " << zz << " yy " << yy << std::endl;
						
						int curr_left_index = 0;
						int current_edge = yy;
						
						
						//target midpoint e.w from paper.
						projection_value_homogeneous_input(proj_mid, V.vertices[ crit_slices[ii+zz].edges[current_edge].midpt ].pt_mp,pi[1]);

						
						sub_mp(denom, proj_top, proj_bottom); // p2(w2) - p2(w0);
						sub_mp(numer, proj_mid, proj_bottom); // p2(e.w) - p2(w0);
						div_mp(md_config.v_target, numer, denom); // [p2(e.w) - p2(w0)] / [p2(w2) - p2(w0)]
						
						//print some display to screen
						{
							print_point_to_screen_matlab(V.vertices[mid_slices[ii].edges[jj].right].pt_mp,"top_start");
							print_point_to_screen_matlab(V.vertices[mid_slices[ii].edges[jj].left].pt_mp,"bottom_start");
							

							print_point_to_screen_matlab(V.vertices[ crit_slices[ii+zz].edges[yy].midpt ].pt_mp,"midpoint_target");
							print_comp_matlab(proj_top,"upper");
							print_comp_matlab(proj_bottom,"lower");
							print_comp_matlab(proj_mid,"mid");
							
							print_comp_matlab(numer,"numer");
							print_comp_matlab(denom,"denom");
							
							print_comp_matlab(md_config.u_target,"u_target");
							print_comp_matlab(md_config.v_target,"v_target");
							
							print_comp_matlab(md_config.crit_val_left,"crit_val_left");
							print_comp_matlab(md_config.crit_val_right,"crit_val_right");
							
							
							set_one_mp(temp);
							sub_mp(temp2, temp, md_config.v_target);
							mul_mp(temp, temp2, proj_bottom);
							
							mul_mp(temp2, md_config.v_target, proj_top);
							
							add_mp(temp3, temp, temp2);
							
							print_comp_matlab(temp3,"proj_1_target_mid");
						}
						
						
						witness_set W_new;
						midpoint_solver_master_entry_point(W_midtrack, // carries with it the start points, and the linears.
																							 &W_new, // new data goes in here
																							 *this,
																							 md_config,
																							 solve_options);
						
						
						//find the edge with crit_slices[ii+zz].edges[yy].midpt index as midpoint
					}
					
				}

			}
			
			
			std::cout << "F.top " << F.top << std::endl;
			std::cout << "F.bottom " << F.bottom << std::endl;
			std::cout << "F.num_left " << F.num_left << std::endl;
			std::cout << "F.num_right " << F.num_right << std::endl;
			add_face(F);
		}
	}
	
	clear_vec_mp(blank_point);
	clear_mp(temp);
	clear_mp(temp2);
	clear_mp(numer);
	clear_mp(denom);
	
	return;
}



void surface_decomposition::print_faces(boost::filesystem::path outputfile)
{
	std::cout << "printing faces to file " << outputfile << std::endl;
	
	FILE *OUT = safe_fopen_write(outputfile);
	
	// output the number of vertices
	fprintf(OUT,"%d\n\n",num_faces);
	
	for(int ii=0;ii<num_faces;ii++){
		fprintf(OUT,"%d %d\n%d %d\n", faces[ii].midpt, faces[ii].index, faces[ii].top, faces[ii].bottom);
		fprintf(OUT,"%ld\n",faces[ii].left.size());
		for (int jj=0; jj<faces[ii].left.size(); jj++) {
			fprintf(OUT,"%d ",faces[ii].left[jj]);
		}
		fprintf(OUT,"\n");
		fprintf(OUT,"%ld\n",faces[ii].right.size());
		for (int jj=0; jj<faces[ii].right.size(); jj++) {
			fprintf(OUT,"%d ",faces[ii].right[jj]);
		}
		fprintf(OUT,"\n");
	}
	fprintf(OUT,"\n");
	
	fclose(OUT);
}













void create_sliced_system(boost::filesystem::path input_file, boost::filesystem::path output_file,
													vec_mp * linears, int num_to_add,
													const witness_set & W)
{
	
	if (W.variable_names.size()==0) {
		std::cout << "trying to create a sliced system, but witness set does not have the variable names." << std::endl;
		deliberate_segfault();
	}
	int *declarations = NULL;
	
	partition_parse(&declarations, input_file, "func_input", "config", 0); // the 0 means not self conjugate.
	free(declarations);
	
	FILE *OUT = safe_fopen_write(output_file);
	FILE *IN = safe_fopen_read("func_input");
	
	
	
	
	fprintf(OUT,"INPUT\n\n");
	copyfile(IN,OUT);
	fclose(IN);
	
	std::vector< int > indicators;
	for (int ii=0; ii<num_to_add; ii++) {
		indicators.push_back(rand());
	}
	
	for (int ii=0; ii<num_to_add; ii++) {
		std::stringstream linname;
		linname << "supp_lin_" << indicators[ii];
		write_vector_as_constants(linears[ii], linname.str(), OUT);
		
		linname.clear();  linname.str("");
	}
	
	for (int ii=0; ii<num_to_add; ii++) {
		fprintf(OUT,"function supp_lin_%d;\n",indicators[ii]);
	}
	for (int ii=0; ii<num_to_add; ii++) {
		fprintf(OUT,"supp_lin_%d = supp_lin_%d_1",indicators[ii],indicators[ii]);
		for (int jj=1; jj<W.num_variables; jj++) {
			fprintf(OUT," + %s*supp_lin_%d_%d",W.variable_names[jj].c_str(), indicators[ii],jj+1);
		}
		fprintf(OUT, ";\n\n");
	}
	fprintf(OUT,"END\n");
	
	fclose(OUT);
	
	
	
	
}

