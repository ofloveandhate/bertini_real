#include "sampler.hpp"








void sampler_configuration::SetDefaults()
{
	no_duplicates = true;
	use_distance_condition = false;

	target_num_samples = 10;

	stifle_membership_screen = 1;
	stifle_text = " > /dev/null ";

	max_num_ribs = 20;
	min_num_ribs = 3;

	minimum_num_iterations = 2;
	maximum_num_iterations = 10;

	mpf_init(TOL);
	mpf_set_d(TOL, 1e-1); // this should be made adaptive to the span of the projection values or the endpoints

	use_gamma_trick = 0;

	mode = Mode::AdaptivePredMovement;

	save_ribs = false;
}



void sampler_configuration::splash_screen()
{
	printf("\n Sampler module for Bertini_real(TM) v%s\n\n", VERSION);
	printf(" s. amethyst, \n with D.J. Bates, W. Hao, \n J.D. Hauenstein, A.J. Sommese, and C.W. Wampler\n\n");
	printf("(using GMP v%d.%d.%d, MPFR v%s)\n\n",
				 __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
	printf("See the website at www.bertinireal.com\n\n");
	printf("Send email to %s for assistance.\n\n",PACKAGE_BUGREPORT);



	bertini_splash_screen();


}


void sampler_configuration::print_usage()
{
	int opt_w = 20;
	int type_w = 10;
	int def_w = 8;
	int note_w = 30;
	auto line = [=](std::string name, std::string type, std::string def_val = " ", std::string note = " ")
		{std::cout << std::setw(opt_w) << std::left << name << std::setw(type_w) << type << std::setw(def_w) << def_val <<  std::setw(note_w) << note << "\n";};

	std::cout << "Sampler has the following options:\n\n";
	line("option","type","default","note");
	std::cout << '\n';
	line("-v -version",  " --  ");
	line("-h -help",  " -- ");
	line("-t -tol -tolerance" ,  "<double>", "0.1", " > 0");
	line("-verb",  "<int>", "0", "how much stuff to print to screen");
	line("-minits",  "<int>", "2","minimum number of passes for adaptive curve or surface refining");
	line("-maxits",  "<int>", "10","maximum number of passes for adaptive curve or surface refining");
	line("-maxribs ",  "<int>", "3","maximum number of ribs for adaptive surface refining");
	line("-minribs",  "<int>", "20", "minimum number of ribs for adaptive surface refining");
	line("-numsamples ",  "<int>", "10", "target number samples per edge");
	line("-mode -m ",  "<char>", "a", "sampling mode.  'a' adaptive by movement, 'd' adaptive by distance, 'f' fixed, ");
	line("-cyclenum",  "<int>", "2", "cycle number to use for rib spacing in face sampling");
	line("-nouniformcyclenum",  " -- ", " ", "turn OFF uniform cycle number usage in surface sampling.  buggy.");
	line("-uniformcyclenum",  " -- ", " ", "turn ON uniform cycle number usage in surface sampling.  works well.");
	line("-saveribs",  " -- ", " ", "turn ON saving of ribs for each face.  off by default.");
	line("-stitchmethod -S", "<char>", "t (trailingangle)", "choose between methods for stitching together triangles. t (trailingangle), s (sumsquaresangles), p (projectionbinning)");
	std::cout << "\n\n\n";
	std::cout.flush();
	return;
}

int  sampler_configuration::parse_commandline(int argc, char **argv)
{
	// this code created based on gnu.org's description of getopt_long
	int choice;

	while (1)
	{
		static struct option long_options[] =
		{
			/* These options set a flag. */
			{"nostifle", no_argument,       0, 's'},
			{"ns", no_argument,       0, 's'},
			{"help",		no_argument,			 0, 'h'},
			{"h",		no_argument,			 0, 'h'},
			{"version",		no_argument,			 0, 'v'},
			{"v",		no_argument,			 0, 'v'},
			{"verb",		required_argument,			 0, 'V'},
			{"tolerance",		required_argument,			 0, 't'},
			{"tol",		required_argument,			 0, 't'},
			{"t",		required_argument,			 0, 't'},
			{"minits",		required_argument,			 0, 'l'},
			{"maxits",		required_argument,			 0, 'm'},
			{"maxribs",		required_argument,			 0, 'R'},
			{"minribs",		required_argument,			 0, 'r'},
			{"numsamples",		required_argument,			 0, 'n'},
			{"nd", no_argument,0,'d'},
			{"m",		required_argument,			 0, 'M'},
			{"mode",		required_argument,			 0, 'M'},
			{"uniformcyclenum", no_argument, 0, 'u'},
			{"nouniformcyclenum", no_argument, 0, 'U'},
			{"cyclenum", required_argument, 0, 'c'},
			{"saveribs", no_argument, 0, 'I'},
			{"stitchmethod", required_argument, 0, 'S'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;

		choice = getopt_long_only (argc, argv, "bdf:svt:V:l:m:R:r:hM:uUc:IS:", // colon requires option, two is optional
															 long_options, &option_index);

		/* Detect the end of the options. */
		if (choice == -1)
			break;

		switch (choice)
		{
			case 'd':
				no_duplicates = false;
				break;


			case 'n':
				target_num_samples = atoi(optarg);

				if (target_num_samples <= 3) {
					std::cout << "The number of desired samples must be larger than 3, but you provided " << target_num_samples << std::endl;
					exit(0);
				}
				break;

			case 's':
				this->stifle_text = "\0";
				break;

			case 'v':
				printf("\n Sampler module for Bertini_real(TM) version %s\n\n", VERSION);
				std::cout << "for help, use option '-h'\n\n";
				exit(0);
				break;

			case 't':

				mpf_set_str(this->TOL,optarg,10);
				break;

			case 'V':
				this->verbose_level(atoi(optarg));
				break;

			case 'l':
				this->minimum_num_iterations = atoi(optarg);
				break;

			case 'm':
				this->maximum_num_iterations = atoi(optarg);
				break;

			case 'M':
			{
				std::string curr_opt{optarg};
				if (curr_opt.size()!=1){
					std::cout << "option to 'mode' must be a single character.  valid modes are 'a' and 'f'\n";
					exit(1);
				}

				switch (curr_opt[0])
				{
					case 'd':
						mode = Mode::AdaptiveConsecDistance;
						break;

					case 'a':
						mode = Mode::AdaptivePredMovement;
						break;

					case 'f':
						mode = Mode::Fixed;
						break;
				}
				break;
			}

			case 'R':
				this->max_num_ribs = atoi(optarg);
				break;

			case 'r':
				this->min_num_ribs = atoi(optarg);
				break;

			case 'h':

				sampler_configuration::print_usage();
				exit(0);

			case 'u':
				this->use_uniform_cycle_num = true;
				break;
			case 'U':
				this->use_uniform_cycle_num = false;
				break;

			case 'c':
				this->cycle_num = atoi(optarg);
				break;

			case 'I':
				this->save_ribs = true;
				break;

			case 'S':{
				std::string curr_opt{optarg};
				switch (curr_opt[0]){
					case 't':
						this->stitch_method = StitchMethod::TrailingAngle;
						std::cout << "using Morgan's triangulate method" << std::endl;
						break;
					case 'p':
						this->stitch_method = StitchMethod::ProjectionBinning;
						break;
					case 's':
						this->stitch_method = StitchMethod::SumOfSquaresAnglesFrom60;
						break;
					default:
						sampler_configuration::print_usage();
						std::cout << "option to stitchmethod or S invalid.  See printed help." << std::endl;
						exit(1);
				}
			}
				
			case '?':
				/* getopt_long already printed an error message. */
				break;

			default:
				sampler_configuration::print_usage();
				exit(1);
		}
	}

	/* Instead of reporting ‘--verbose’
	 and ‘--brief’ as they are encountered,
	 we report the final status resulting from them. */


	/* Print any remaining command line arguments (not options). */
	if (optind < argc)
	{
		printf ("non-option ARGV-elements: ");
		while (optind < argc)
			printf ("%s ", argv[optind++]);
		putchar ('\n');
	}

	return 0;
}







int main(int argC, char *args[])
{
	MPI_Init(&argC,&args);

	sampler_configuration sampler_options;
	sampler_options.parse_commandline(argC, args);


	if (false) { // sampler_options.debugwait()

		if (sampler_options.is_head()) {
			std::cout << "in debug mode, waiting to start so you can attach a debugger to this process" << std::endl;
			std::cout << "master PID: " << getpid() << std::endl;
			int duration = 15;
			for (int ii=0; ii<duration; ii++) {
				std::cout << "starting program in " << duration-ii << " seconds" << std::endl;
				sleep(1);
				std::cout << "\033[F"; // that's a line up movement.  only works with some terminals
			}
		}


		if (sampler_options.use_parallel()) {
			MPI_Barrier(MPI_COMM_WORLD);
		}

	}



	if (sampler_options.is_head())
	{
		boost::timer::auto_cpu_timer t;
		SamplerMaster(sampler_options);
	}
	else
	{
		SamplerWorker(sampler_options);
	}

	clearMP();
	MPI_Finalize();

	return 0;
}



void SamplerMaster(sampler_configuration & sampler_options)
{
	sampler_options.splash_screen();

	boost::filesystem::path inputName, witnessSetName, samplingNamenew;
	int MPType, dimension;

	boost::filesystem::path directoryName;

	get_dir_mptype_dimen( directoryName, MPType, dimension); // i really do hate this.
	sampler_options.output_dir(directoryName);
	witnessSetName = directoryName / "WitnessSet";
	samplingNamenew = directoryName;




	//only one will be used.  i don't know how to avoid this duplicity.  sorry.
	// no, wait, i do.  it's by using polymorphism, i think.  that's not a pressing matter to me.
	Curve curve;
	Surface surf;

	Decomposition * decom_pointy; // this feels unnecessary
	switch (dimension) {
		case 1:
		{
			curve.setup(directoryName);
			decom_pointy = &curve;
			break;
		}

		case 2:
		{
			surf.setup(directoryName);
			decom_pointy = &surf;
			break;
		}
		default:
		{
			throw std::runtime_error("invalid dimension of object to sample.  you sure you have a valid decomposition?");
		}
	}



	SolverConfiguration solve_options;

	common_sampler_startup(*decom_pointy,
						   sampler_options,
						   solve_options);


	if (solve_options.use_parallel())
	{
		int routine = TRACKER_CONFIG;
		MPI_Bcast(&routine, 1, MPI_INT, solve_options.head(), MPI_COMM_WORLD);
		bcast_tracker_config_t(&solve_options.T, solve_options.id(), solve_options.head() );
	}


	VertexSet V(decom_pointy->num_variables());

	V.set_tracker_config(&solve_options.T);
	V.setup_vertices(directoryName / "V.vertex"); //setup V structure from V.vertex
	V.set_same_point_tolerance(1e1*solve_options.T.real_threshold);


	/////////
	////////
	//////
	////
	//
	//  Generate new sampling data
	//


	switch (dimension) {
		case 1:
		{
			switch (sampler_options.mode){
				case sampler_configuration::Mode::Fixed:
					curve.FixedSampler(V,
									sampler_options,
									solve_options,
									sampler_options.target_num_samples);
					break;
				case sampler_configuration::Mode::AdaptiveConsecDistance:

					curve.AdaptiveDistanceSampler(V,
												sampler_options,
												solve_options);
					break;

				case sampler_configuration::Mode::AdaptivePredMovement:
					curve.AdaptiveMovementSampler(V,
												sampler_options,
												solve_options);
					break;

				case sampler_configuration::Mode::SemiFixed:
					throw std::runtime_error("semi fixed not available for root-level curves.");
					break;
			} // switch
			curve.output_sampling_data(directoryName);
			V.print(directoryName / "V_samp.vertex");

			break;
		}


		case 2:
		{
			switch (sampler_options.mode){
				case sampler_configuration::Mode::Fixed:
				{
					surf.fixed_sampler(V,
											 sampler_options,
											 solve_options);

					break;
				}
				case sampler_configuration::Mode::AdaptivePredMovement:
					std::cout << color:: magenta() << "adaptive by movement not implemented for surfaces, using adaptive by distance" << color::console_default() << "\n\n";
				case sampler_configuration::Mode::AdaptiveConsecDistance:
				{
					surf.AdaptiveSampler(V,
											 sampler_options,
											 solve_options);
					break;
				}

				case sampler_configuration::Mode::SemiFixed:
					throw std::runtime_error("semi fixed not available for root-level curves.");
					break;

			} // switch

			surf.output_sampling_data(directoryName);
			V.print(directoryName / "V_samp.vertex");

			break;
		}
		default:
			break;
	}

	// dismiss the workers
	int sendme = TERMINATE;
	MPI_Bcast(&sendme, 1, MPI_INT, 0, MPI_COMM_WORLD);
}


void SamplerWorker(sampler_configuration & sampler_options)
{

	int routine = INITIAL_STATE;

	SolverConfiguration solve_options;

	solve_options.robust = sampler_options.robustness()>=1;
	
	bool have_tracker_config = false;



	while (routine != TERMINATE) {
		// get the task to perform.  yes, everyone has to participate, because of B1 things.  ugh.
		MPI_Bcast(&routine, 1, MPI_INT, solve_options.head(), MPI_COMM_WORLD);

		if ( (solve_options.id()==1) && (sampler_options.verbose_level()>=2)) { //(routine!=0) &&
			std::cout << "sampler worker received call for help for solver " << enum_lookup(routine) << " (code " << routine << ")" << std::endl;
		}

		switch (routine) {
			case SAMPLE_CURVE:
			{
				WorkerSampleCurve(sampler_options, solve_options);
				break;
			}
			case SAMPLE_SURFACE:
			{
				WorkerSampleSurface(sampler_options, solve_options);
				break;
			}
			case PARSING:
			{
				int single_int_buffer = 0;
				MPI_Bcast(&single_int_buffer, 1, MPI_INT, solve_options.head(), solve_options.comm()); // this catches a broadcast from Bertini's parser...
				break;
			}
			case TRACKER_CONFIG:
			{
				if (have_tracker_config)
					throw std::runtime_error("worker getting tracker config AGAIN");

				bcast_tracker_config_t(&solve_options.T, solve_options.id(), solve_options.head() );
				have_tracker_config = true;
				initMP(solve_options.T.Precision);
				break;
			}
			case TERMINATE:
				break;
			default:
				std::cout << "received unknown routine: " << enum_lookup(routine) << std::endl;
				break;
		}
	}
}




void common_sampler_startup(const Decomposition & D,
							sampler_configuration & sampler_options,
							SolverConfiguration & solve_options)
{

	parse_input_file(D.input_filename()); // restores all the temp files generated by the parser, to this folder.  i think this can be removed?  but the PPD and tracker config both depend on it...  so maybe not.



	// set up the solver configuration
	get_tracker_config(solve_options, solve_options.T.MPType);

	initMP(solve_options.T.Precision);


	parse_preproc_data("preproc_data", &solve_options.PPD);





	solve_options.verbose_level(sampler_options.verbose_level());



	solve_options.T.ratioTol = 1; // manually assert to be more permissive.  i don't really like this.

	solve_options.use_midpoint_checker = 0;
}








//dehomogenizes, takes the average, computes the projection.
//takes in the full projection \pi, including the homogenizing coordinate.
void estimate_new_projection_value(comp_mp result, vec_mp left, vec_mp right, vec_mp pi){
	int ii;

	if (left->size < pi->size) {
		printf("left point too short in estimate new projection value\n");
		deliberate_segfault();
	}

	if (right->size < pi->size) {
		printf("left point too short in estimate new projection value\n");
		deliberate_segfault();
	}
    //	print_point_to_screen_matlab(left,"left");
    //	print_point_to_screen_matlab(right,"right");
    //	print_point_to_screen_matlab(pi,"pi");

	vec_mp dehom_left, dehom_right;
	init_vec_mp(dehom_left,pi->size-1);   dehom_left->size  = pi->size-1;
	init_vec_mp(dehom_right,pi->size-1);  dehom_right->size = pi->size-1;

	dehomogenize(&dehom_left,left,pi->size);
	dehomogenize(&dehom_right,right,pi->size);


	comp_mp temp, temp2, half; init_mp(temp); init_mp(temp2);  init_mp(half);

	mpf_set_d(half->r, 0.5); mpf_set_d(half->i, 0.0);


	set_zero_mp(result);                                           // result = 0;  initialize

	for (ii = 0; ii<pi->size-1; ii++) {
		add_mp(temp,&dehom_left->coord[ii],&dehom_right->coord[ii]); //  a = (x+y)
		mul_mp(temp2, temp, half);                                   //  b = a/2
		mul_mp(temp,&pi->coord[ii+1],temp2);                          //  a = b.pi
		set_mp(temp2,result);                                        //  b = result
		add_mp(result, temp, temp2);                                  //  result = a+b
	}
	// in other words, result += (x+y)/2 \cdot pi


    real_threshold(result,1e-9);


	clear_mp(temp); clear_mp(temp2); clear_mp(half);
	clear_vec_mp(dehom_right);clear_vec_mp(dehom_left);

	return;
}



//dehomogenizes, takes the average, computes the projection.
//takes in the full projection \pi, including the homogenizing coordinate.
// this version returns the estimated point as well.
void estimate_new_projection_value(comp_mp result, vec_mp estimated_point, vec_mp left, vec_mp right, vec_mp pi){


	if (left->size != right->size) {
		printf("attempting to estimate_new_projection_value on vectors of different size\n");
		br_exit(98128);
	}

	comp_mp temp;  init_mp(temp);

	comp_mp half; init_mp(half);
	mpf_set_str(half->r, "0.5", 10); mpf_set_str(half->i, "0.0", 10);



	vec_mp dehom_left, dehom_right;
	init_vec_mp(dehom_left,left->size-1);   dehom_left->size  = left->size-1;
	init_vec_mp(dehom_right,right->size-1); dehom_right->size = right->size-1;

	dehomogenize(&dehom_left,left,pi->size);
	dehomogenize(&dehom_right,right,pi->size);



	vec_add_mp(estimated_point, dehom_left, dehom_right);
	vec_mulcomp_mp(estimated_point, estimated_point, half);

	set_zero_mp(result);                                           // result = 0;  initialize

	for (int ii = 0; ii<pi->size-1; ii++) {                                 //  b = a/2
		mul_mp(temp,&pi->coord[ii+1],&estimated_point->coord[ii]);                          //  a = b.pi
		add_mp(result, result, temp);                                  //  result = a+b
	}

	// in other words, result += (x+y)/2 \cdot pi


	//  i think this thresholding should be moved to outside this call.
	mpf_t zerothresh; mpf_init(zerothresh);
	mpf_set_d(zerothresh, 1e-9);
	if (mpf_cmp(result->i, zerothresh) < 0){
		mpf_set_str(result->i, "0.0", 10);
	}

	mpf_clear(zerothresh);

	clear_mp(temp); clear_mp(half);

	clear_vec_mp(dehom_right);clear_vec_mp(dehom_left);

	return;
}





void UnpackProjvals(vec_mp projvals, const std::vector< int > & rib, int projind, const VertexSet & V)
{
	for (int ii=0; ii< rib.size(); ++ii)
		set_mp(&(projvals->coord[ii]), &(V[rib[ii]].projection_values()->coord[projind]));
}

void ScaleToUnitInverval(vec_mp projvals)
{
	for (int ii=projvals->size-1; ii>=0; --ii)
		sub_mp(&(projvals->coord[ii]), &(projvals->coord[ii]), &(projvals->coord[0]));

	auto ind = projvals->size-1;
	for (int ii=1; ii< projvals->size; ++ii)
		div_mp(&(projvals->coord[ii]), &(projvals->coord[ii]), &(projvals->coord[ind]));

}


void TailEndOfRibs(const std::vector< int > & rib1, const std::vector< int > & rib2, int curr_index_rib1, int curr_index_rib2, std::vector< Triangle> & current_samples)
{
	const std::vector< int > *exhausted_rib, *rib_still_going;
	unsigned int index_still_going, terminal_index;
	bool flip;

	if (curr_index_rib1==rib1.size()-1) {
		exhausted_rib = &rib1;
		rib_still_going = &rib2;
		index_still_going = curr_index_rib2;
		terminal_index = curr_index_rib1;
		flip = false;
	}
	else
	{
		exhausted_rib = &rib2;
		rib_still_going = &rib1;
		index_still_going = curr_index_rib1;
		terminal_index = curr_index_rib2;
		flip = true;
	}


	for (; index_still_going < rib_still_going->size()-1; index_still_going++) { // initializer deliberately empty

		long long v1, v2, v3;
		if (flip) {
			v1 = rib_still_going->at(index_still_going);
			v2 = rib_still_going->at(index_still_going+1);
		}
		else
		{
			v1 = rib_still_going->at(index_still_going+1);
			v2 = rib_still_going->at(index_still_going);
		}
		v3 = exhausted_rib->at(terminal_index);

		current_samples.push_back( Triangle(v1,v2,v3) );
	}
}


void triangulate_two_ribs_by_projection_binning(const std::vector< int > & rib1, const std::vector< int > & rib2,
											  VertexSet & V, double real_thresh,
											  std::vector< Triangle> & current_samples)
{
#ifdef functionentry_output
	std::cout << "triangulate_two_ribs_by_projection_binning" << std::endl;
#endif

	bool bail_out = false;

	if (rib1.size()==0) {
		std::cout << "rib1 had 0 size!" << std::endl;
		bail_out = true;
	}
	if (rib2.size()==0) {
		std::cout << "rib2 had 0 size!" << std::endl;
		bail_out = true;
	}

	if (rib1.size()==1 && rib2.size()==1) {
		std::cout << "both ribs have size 1!" << std::endl;
		bail_out = true;
	}

	if (bail_out) {
		return;
	}

	int num_vars = V.num_natural_variables();

	vec_mp projvals1, projvals2;
	init_vec_mp(projvals1, rib1.size()); projvals1->size = rib1.size();
	init_vec_mp(projvals2, rib2.size()); projvals2->size = rib2.size();

	UnpackProjvals(projvals1, rib1, 1, V);
	ScaleToUnitInverval(projvals1);

	UnpackProjvals(projvals2, rib2, 1, V);
	ScaleToUnitInverval(projvals2);

	vec_mp *pi_long, *pi_short;
	const std::vector<int> *rib_long, *rib_short;
	if (projvals1->size >= projvals2->size)
	{
		pi_long = &projvals1; rib_long = &rib1;
		pi_short = &projvals2;	rib_short = &rib2;
	}
	else
	{
		pi_long = &projvals2; rib_long = &rib2;
		pi_short = &projvals1;	rib_short = &rib1;
	}
	// ok now we have the projection values.  we're going to work from the outside in, left and right simultaneously, to connect the triangles.


	// print_point_to_screen_matlab(*pi_short, "pi_short");
	// print_point_to_screen_matlab(*pi_long, "pi_long");

	// std::cout << "t = [...\n";
	int Q = 1;


	for (int ii=1; ii<(*pi_short)->size; ii++)
	{
		int I = ii-1;

		while (Q < (*pi_long)->size && mpf_cmp( (*pi_long)->coord[Q].r,(*pi_short)->coord[ii].r)<0)
		{
			current_samples.push_back(
									  Triangle(
											   (*rib_short)[I], 
											   (*rib_long)[Q], 
											   (*rib_long)[Q-1]
											   ) 
									  );
			Q++;
		}

		// then finish up the last one to advance
		current_samples.push_back(
								  Triangle(
										   (*rib_short)[I], 
										   (*rib_short)[ii],
										   (*rib_long)[Q-1]
										   ) 
								  );
	}

	TailEndOfRibs((*rib_short), (*rib_long), rib_short->size()-1, Q-1, current_samples);

	clear_vec_mp(projvals1); clear_vec_mp(projvals2);
}



void triangulate_two_ribs_by_trailing_angle(const std::vector< int > & rib1, const std::vector< int > & rib2,
											  VertexSet & V, double real_thresh,
											  std::vector< Triangle> & current_samples)
{
#ifdef functionentry_output
	std::cout << "triangulate_two_ribs_by_trailing_angle" << std::endl;
#endif

	

	bool bail_out = false;

	if (rib1.size()==0) {
		std::cout << "rib1 had 0 size!" << std::endl;
		bail_out = true;
	}
	if (rib2.size()==0) {
		std::cout << "rib2 had 0 size!" << std::endl;
		bail_out = true;
	}

	if (rib1.size()==1 && rib2.size()==1) {
		std::cout << "both ribs have size 1!" << std::endl;
		bail_out = true;
	}

	if (bail_out) {
		return;
	}


	int num_vars = V.num_natural_variables();

	unsigned int num_vectors_needed = 9;
	vec_mp * bulk_vectors = (vec_mp *) br_malloc(num_vectors_needed* sizeof(vec_mp));

	for (unsigned int ii=0; ii<num_vectors_needed; ii++) {
		init_vec_mp(bulk_vectors[ii],num_vars);  bulk_vectors[ii]->size = num_vars;

	}

	vec_mp *A = &bulk_vectors[0], *B = &bulk_vectors[1], *C = &bulk_vectors[2], *D = &bulk_vectors[3];


	vec_mp *AB = &bulk_vectors[4];
	vec_mp *AC = &bulk_vectors[5];
	vec_mp *BC = &bulk_vectors[6];
	vec_mp *DA = &bulk_vectors[7];
	vec_mp *DC = &bulk_vectors[8];



	// seed the advancing loop.
	dehomogenize(B, V[rib1[0]].point(), num_vars);
	(*B)->size = num_vars-1;
	real_threshold(*B,real_thresh);
	dehomogenize(D, V[rib2[0]].point(), num_vars);
	(*D)->size = num_vars-1;
	real_threshold(*D,real_thresh);


	comp_mp cos_angle_CAB, cos_angle_BCA, cos_angle_ABC;
	init_mp(cos_angle_CAB); init_mp(cos_angle_BCA); init_mp(cos_angle_ABC);

	comp_mp cos_angle_DCA, cos_angle_ADC, cos_angle_CAD;
	init_mp(cos_angle_DCA); init_mp(cos_angle_ADC); init_mp(cos_angle_CAD);


	comp_mp length_AB, length_AC, length_BC, length_DA, length_DC;
	init_mp(length_AB); init_mp(length_AC); init_mp(length_BC); init_mp(length_DA); init_mp(length_DC);

	comp_mp dot_CAB, dot_BCA, dot_ABC;
	init_mp(dot_CAB); init_mp(dot_BCA); init_mp(dot_ABC);

	comp_mp dot_DCA, dot_ADC, dot_CAD;
	init_mp(dot_DCA); init_mp(dot_ADC); init_mp(dot_CAD);

	comp_mp temp;  init_mp(temp);
	comp_d temp_d;

	unsigned int curr_index_rib1 = 0, curr_index_rib2 = 0;
	bool moved_1 = true, moved_2 = true;  //this is an intial condition to get them both set properly.  all subsequent iterations have only one as moved==true, and the other is always false.
	while (curr_index_rib1 < rib1.size()-1 // neither rib size is 0, so this -1 is ok, won't underflow
		   &&
		   curr_index_rib2 < rib2.size()-1)
	{

#ifdef debug_compile
		int a =rib1[curr_index_rib1];
		int b =rib1[curr_index_rib1+1];
		int c =rib2[curr_index_rib2];
		int d =rib2[curr_index_rib2+1];

		std::cout << rib1[curr_index_rib1] << " " << rib1[curr_index_rib1+1] << std::endl;
		std::cout << rib2[curr_index_rib2] << " " << rib2[curr_index_rib2+1] << std::endl;
#endif


		dehomogenize(A, V[rib1[curr_index_rib1]].point(), num_vars); (*A)->size = num_vars-1;
		dehomogenize(B, V[rib1[curr_index_rib1+1]].point(), num_vars); (*B)->size = num_vars-1;

		real_threshold(*A,real_thresh);
		real_threshold(*B,real_thresh);


		vec_sub_mp(*AB, *B,*A);

		twoNormVec_mp(*AB, length_AB);




		dehomogenize(C, V[rib2[curr_index_rib2]].point(), num_vars); (*C)->size = num_vars-1;
		dehomogenize(D, V[rib2[curr_index_rib2+1]].point(), num_vars); (*D)->size = num_vars-1;

		real_threshold(*C,real_thresh);
		real_threshold(*D,real_thresh);

		vec_sub_mp(*DC, *C,*D);
		twoNormVec_mp(*DC, length_DC);




//		A           --->           B
//		  ***********************
//		  *- <--              *
//		  * -  --           *
//		  *  -   --       *    --
//		| *   -         *   ---
//		| *    -      *  <--
//		| *     -   *
//		\/*      -*
//		  *     * -
//		  *   *    -
//		  * *       -
//		  *-----------
//		C    <--       D



//        AC, BC, DA,

		vec_sub_mp(*AC, *C,*A);
		vec_sub_mp(*BC, *C,*B);
		vec_sub_mp(*DA, *A,*D);


#ifdef debug_compile
		print_point_to_screen_matlab(*A,"A");
		print_point_to_screen_matlab(*B,"B");
		print_point_to_screen_matlab(*C,"C");
		print_point_to_screen_matlab(*D,"D");

		print_point_to_screen_matlab(*AB,"AB");
		print_point_to_screen_matlab(*AC,"AC");
		print_point_to_screen_matlab(*BC,"BC");
		print_point_to_screen_matlab(*DA,"DA");
		print_point_to_screen_matlab(*DC,"DC");
#endif


		// now have the 5 vectors for the test computed.  (5 because the two triangles share a common leg)

		dot_product_mp(dot_CAB, *AC, *AB);

		dot_product_mp(dot_BCA, *AC, *BC);

		dot_product_mp(dot_ABC, *AB, *BC);
		neg_mp(dot_ABC,dot_ABC);

		dot_product_mp(dot_DCA, *DC, *AC);

		dot_product_mp(dot_ADC, *DA, *DC);


		dot_product_mp(dot_CAD, *DA, *AC);
		neg_mp(dot_CAD,dot_CAD);

		twoNormVec_mp(*AC, length_AC);
		twoNormVec_mp(*BC, length_BC);
		twoNormVec_mp(*DA, length_DA);


		

		//CAB
		double error_CAB = compute_abs_of_difference_from_sixtydegrees(temp, length_AC, length_AB, dot_CAB);

		//BCA
		double error_BCA = compute_abs_of_difference_from_sixtydegrees(temp, length_BC, length_AC, dot_BCA);

		//ABC
		double error_ABC = compute_abs_of_difference_from_sixtydegrees(temp, length_BC, length_AB, dot_ABC);


		//DCA
		double error_DCA = compute_abs_of_difference_from_sixtydegrees(temp, length_DC, length_AC, dot_DCA);

		//ADC
		double error_ADC = compute_abs_of_difference_from_sixtydegrees(temp, length_DA, length_DC, dot_ADC);

		//CAD
		double error_CAD = compute_abs_of_difference_from_sixtydegrees(temp, length_AC, length_DA, dot_CAD);




#ifdef debug_compile
		print_comp_matlab(dot_CAB,"CAB");
		std::cout << c << " " << a << " " << b << " " << total_error_rib1 << std::endl;
#endif

#ifdef debug_compile
		print_comp_matlab(dot_BCA,"BCA");
		std::cout << b << " " << c << " " << a << " " << angle_BCA << std::endl;
#endif
#ifdef debug_compile
		print_comp_matlab(dot_ABC,"ABC");
		std::cout << a << " " << b << " " << c << " " << angle_ABC << std::endl;
#endif
#ifdef debug_compile
		print_comp_matlab(dot_DCA,"DCA");
		std::cout << d << " " << c << " " << a << " " << total_error_rib2 << std::endl;
#endif
#ifdef debug_compile
		print_comp_matlab(dot_ADC,"ADC");
		std::cout << a << " " << d << " " << c << " " << angle_ADC << std::endl;
#endif
#ifdef debug_compile
		print_comp_matlab(dot_CAD,"CAD");
		std::cout << c << " " << a << " " << d << " " << angle_CAD << std::endl;
#endif




		// logic for which one to advance
		int advance = 0;
		if (error_CAB < error_DCA){
			advance = 1;
		}
		else
		{	advance = 2;}




		if (advance==1) { // if the 1 Triangle is more equilateral than the 2 Triangle.
			current_samples.push_back(
									  Triangle(// Triangle A B C
											   rib1[curr_index_rib1], //A
											   rib1[curr_index_rib1+1], //B
											   rib2[curr_index_rib2]) //C
									  );
			moved_1 = true;  curr_index_rib1++;
			moved_2 = false;

		}
		else
		{
			current_samples.push_back(
									  Triangle(// Triangle C A D
											   rib2[curr_index_rib2], //C
											   rib1[curr_index_rib1], //A
											   rib2[curr_index_rib2+1]) //D
									  );
			moved_1 = false;
			moved_2 = true; curr_index_rib2++;
		}

	} // re: while loop.



	// now down here, we have triangulated until one of the ribs is on its last point, so there is no more testing that can be done.  you simply have to connect the rest into triangles.

	TailEndOfRibs(rib1, rib2, curr_index_rib1, curr_index_rib2, current_samples);




	clear_mp(cos_angle_CAB); clear_mp(cos_angle_BCA); clear_mp(cos_angle_ABC);
	clear_mp(cos_angle_DCA); clear_mp(cos_angle_ADC); clear_mp(cos_angle_CAD);


	clear_mp(length_AB); clear_mp(length_AC); clear_mp(length_BC); clear_mp(length_DA); clear_mp(length_DC);
	clear_mp(dot_CAB); clear_mp(dot_BCA); clear_mp(dot_ABC);
	clear_mp(dot_DCA); clear_mp(dot_ADC); clear_mp(dot_CAD);

	clear_mp(temp);
	// clean up at the end.  i wish scope was deletion!

	for (unsigned int ii=0; ii<num_vectors_needed; ii++) {
		clear_vec_mp(bulk_vectors[ii]);
	}
	free(bulk_vectors);



	return;
}




void triangulate_two_ribs_by_angle_optimization(const std::vector< int > & rib1, const std::vector< int > & rib2,
											  VertexSet & V, double real_thresh,
											  std::vector< Triangle> & current_samples)
{
#ifdef functionentry_output
	std::cout << "triangulate_two_ribs_by_angle_optimization" << std::endl;
#endif


	bool bail_out = false;

	if (rib1.size()==0) {
		std::cout << "rib1 had 0 size!" << std::endl;
		bail_out = true;
	}
	if (rib2.size()==0) {
		std::cout << "rib2 had 0 size!" << std::endl;
		bail_out = true;
	}

	if (rib1.size()==1 && rib2.size()==1) {
		std::cout << "both ribs have size 1!" << std::endl;
		bail_out = true;
	}

	if (bail_out) {
		return;
	}


	int num_vars = V.num_natural_variables();

	unsigned int num_vectors_needed = 9;
	vec_mp * bulk_vectors = (vec_mp *) br_malloc(num_vectors_needed* sizeof(vec_mp));

	for (unsigned int ii=0; ii<num_vectors_needed; ii++) {
		init_vec_mp(bulk_vectors[ii],num_vars);  bulk_vectors[ii]->size = num_vars;

	}

	vec_mp *A = &bulk_vectors[0], *B = &bulk_vectors[1], *C = &bulk_vectors[2], *D = &bulk_vectors[3];


	vec_mp *AB = &bulk_vectors[4];
	vec_mp *AC = &bulk_vectors[5];
	vec_mp *BC = &bulk_vectors[6];
	vec_mp *DA = &bulk_vectors[7];
	vec_mp *DC = &bulk_vectors[8];



	// seed the advancing loop.
	dehomogenize(B, V[rib1[0]].point(), num_vars);
	(*B)->size = num_vars-1;
	real_threshold(*B,real_thresh);
	dehomogenize(D, V[rib2[0]].point(), num_vars);
	(*D)->size = num_vars-1;
	real_threshold(*D,real_thresh);


	comp_mp cos_angle_CAB, cos_angle_BCA, cos_angle_ABC;
	init_mp(cos_angle_CAB); init_mp(cos_angle_BCA); init_mp(cos_angle_ABC);

	comp_mp cos_angle_DCA, cos_angle_ADC, cos_angle_CAD;
	init_mp(cos_angle_DCA); init_mp(cos_angle_ADC); init_mp(cos_angle_CAD);


	comp_mp length_AB, length_AC, length_BC, length_DA, length_DC;
	init_mp(length_AB); init_mp(length_AC); init_mp(length_BC); init_mp(length_DA); init_mp(length_DC);

	comp_mp dot_CAB, dot_BCA, dot_ABC;
	init_mp(dot_CAB); init_mp(dot_BCA); init_mp(dot_ABC);

	comp_mp dot_DCA, dot_ADC, dot_CAD;
	init_mp(dot_DCA); init_mp(dot_ADC); init_mp(dot_CAD);

	comp_mp temp;  init_mp(temp);
	comp_d temp_d;

	unsigned int curr_index_rib1 = 0, curr_index_rib2 = 0;
	bool moved_1 = true, moved_2 = true;  //this is an intial condition to get them both set properly.  all subsequent iterations have only one as moved==true, and the other is always false.
	while (curr_index_rib1 < rib1.size()-1 // neither rib size is 0, so this -1 is ok, won't underflow
		   &&
		   curr_index_rib2 < rib2.size()-1)
	{

#ifdef debug_compile
		int a =rib1[curr_index_rib1];
		int b =rib1[curr_index_rib1+1];
		int c =rib2[curr_index_rib2];
		int d =rib2[curr_index_rib2+1];

		std::cout << rib1[curr_index_rib1] << " " << rib1[curr_index_rib1+1] << std::endl;
		std::cout << rib2[curr_index_rib2] << " " << rib2[curr_index_rib2+1] << std::endl;
#endif

		if (moved_1) {
			vec_mp * temp_vec = A; // swap
			A = B;
			B = temp_vec;
			dehomogenize(B, V[rib1[curr_index_rib1+1]].point(), num_vars);
			(*B)->size = num_vars-1;
			real_threshold(*B,real_thresh);


			vec_sub_mp(*AB, *B,*A);



			twoNormVec_mp(*AB, length_AB);
		}

		if (moved_2){ // moved_smaller
			vec_mp * temp_vec = C; // swap
			C = D;
			D = temp_vec;
			dehomogenize(D, V[rib2[curr_index_rib2+1]].point(), num_vars);
			(*D)->size = num_vars-1;

			real_threshold(*D,real_thresh);

			vec_sub_mp(*DC, *C,*D);
			twoNormVec_mp(*DC, length_DC);
		}



//		A           --->           B
//		  ***********************
//		  *- <--              *
//		  * -  --           *
//		  *  -   --       *    --
//		| *   -         *   ---
//		| *    -      *  <--
//		| *     -   *
//		\/*      -*
//		  *     * -
//		  *   *    -
//		  * *       -
//		  *-----------
//		C    <--       D



//        AC, BC, DA,

		vec_sub_mp(*AC, *C,*A);
		vec_sub_mp(*BC, *C,*B);
		vec_sub_mp(*DA, *A,*D);


#ifdef debug_compile
		print_point_to_screen_matlab(*A,"A");
		print_point_to_screen_matlab(*B,"B");
		print_point_to_screen_matlab(*C,"C");
		print_point_to_screen_matlab(*D,"D");

		print_point_to_screen_matlab(*AB,"AB");
		print_point_to_screen_matlab(*AC,"AC");
		print_point_to_screen_matlab(*BC,"BC");
		print_point_to_screen_matlab(*DA,"DA");
		print_point_to_screen_matlab(*DC,"DC");
#endif


		// now have the 5 vectors for the test computed.  (5 because the two triangles share a common leg)

		dot_product_mp(dot_CAB, *AC, *AB);

		dot_product_mp(dot_BCA, *AC, *BC);

		dot_product_mp(dot_ABC, *AB, *BC);
		neg_mp(dot_ABC,dot_ABC);

		dot_product_mp(dot_DCA, *DC, *AC);

		dot_product_mp(dot_ADC, *DA, *DC);


		dot_product_mp(dot_CAD, *DA, *AC);
		neg_mp(dot_CAD,dot_CAD);

		twoNormVec_mp(*AC, length_AC);
		twoNormVec_mp(*BC, length_BC);
		twoNormVec_mp(*DA, length_DA);


		double thresh = 3;
		int advance = 0;
		bool aspect_ok_rib1 = true, aspect_ok_rib2 = true;

		div_mp(temp,length_BC,length_AB);
		mp_to_d(temp_d, temp);
//		double aspect11 = temp_d->r;
		if ( temp_d->r  > thresh) {
			aspect_ok_rib1 = false;
		}
		div_mp(temp,length_BC,length_AC);
		mp_to_d(temp_d, temp);
//		double aspect12 = temp_d->r;
		if ( temp_d->r  > thresh) {
			aspect_ok_rib1 = false;
		}
//


		div_mp(temp,length_DA,length_DC);
		mp_to_d(temp_d, temp);
//		double aspect21 = temp_d->r;
		if (temp_d->r > thresh) {
			aspect_ok_rib2 = false;
		}
		div_mp(temp,length_DA,length_AC);
		mp_to_d(temp_d, temp);
//		double aspect22 = temp_d->r;
		if (temp_d->r > thresh) {
			aspect_ok_rib2 = false;
		}






		// the below lines computes the sum of the squares of the differences of the absolute values of the cosines of two of the angles in one of the triangles.

		//CAB
		double total_error_rib1 = compute_square_of_difference_from_sixtydegrees(temp, length_AC, length_AB, dot_CAB);
#ifdef debug_compile
		print_comp_matlab(dot_CAB,"CAB");
		std::cout << c << " " << a << " " << b << " " << total_error_rib1 << std::endl;
#endif


		//BCA
		double angle_BCA = compute_square_of_difference_from_sixtydegrees(temp, length_BC, length_AC, dot_BCA);
#ifdef debug_compile
		print_comp_matlab(dot_BCA,"BCA");
		std::cout << b << " " << c << " " << a << " " << angle_BCA << std::endl;
#endif
		total_error_rib1 += angle_BCA;

		//ABC
		double angle_ABC = compute_square_of_difference_from_sixtydegrees(temp, length_BC, length_AB, dot_ABC);
		total_error_rib1 += angle_ABC;
#ifdef debug_compile
		print_comp_matlab(dot_ABC,"ABC");
		std::cout << a << " " << b << " " << c << " " << angle_ABC << std::endl;
#endif






		//DCA
		double total_error_rib2 = compute_square_of_difference_from_sixtydegrees(temp, length_DC, length_AC, dot_DCA);
#ifdef debug_compile
		print_comp_matlab(dot_DCA,"DCA");
		std::cout << d << " " << c << " " << a << " " << total_error_rib2 << std::endl;
#endif

		//ADC

		double angle_ADC = compute_square_of_difference_from_sixtydegrees(temp, length_DA, length_DC, dot_ADC);

		total_error_rib2+= angle_ADC;
#ifdef debug_compile
		print_comp_matlab(dot_ADC,"ADC");
		std::cout << a << " " << d << " " << c << " " << angle_ADC << std::endl;
#endif
		//CAD
		double angle_CAD = compute_square_of_difference_from_sixtydegrees(temp, length_AC, length_DA, dot_CAD);
		total_error_rib2 += angle_CAD;

#ifdef debug_compile
		print_comp_matlab(dot_CAD,"CAD");
		std::cout << c << " " << a << " " << d << " " << angle_CAD << std::endl;
#endif





		if ((total_error_rib1>0.25) && (total_error_rib2>0.25)){
			// both triangles more or less equilateral, or both are bad
			if ( mpf_cmp(length_DA->r,length_BC->r)<0)
				advance = 2;
			else
				advance = 1;
		}
		else{
			if (total_error_rib1 < total_error_rib2) {
				advance = 1;
			}
			else
				advance = 2;
		}


//		if(!aspect_ok_rib1 && aspect_ok_rib2) {
//			advance = 2;
//		}
//		else if(!aspect_ok_rib2 && aspect_ok_rib1) {
//			advance = 1;
//		}
//		else if ( (!aspect_ok_rib2 && !aspect_ok_rib1) )
//		{
//
//
//
//
//
//		}
//		else
//		{
//
//
//			if (total_error_rib1 < total_error_rib2) {
//				advance = 1;
//			}
//			else{
//				advance = 2;
//			}
//		}




		if (advance==1) { // if the 1 Triangle is more equilateral than the 2 Triangle.
			current_samples.push_back(
									  Triangle(// Triangle A B C
											   rib1[curr_index_rib1], //A
											   rib1[curr_index_rib1+1], //B
											   rib2[curr_index_rib2]) //C
									  );
			moved_1 = true;  curr_index_rib1++;
			moved_2 = false;

		}
		else
		{
			current_samples.push_back(
									  Triangle(// Triangle C A D
											   rib2[curr_index_rib2], //C
											   rib1[curr_index_rib1], //A
											   rib2[curr_index_rib2+1]) //D
									  );
			moved_1 = false;
			moved_2 = true; curr_index_rib2++;
		}

	} // re: while loop.



	// now down here, we have triangulated until one of the ribs is on its last point, so there is no more testing that can be done.  you simply have to connect the rest into triangles.

	TailEndOfRibs(rib1, rib2, curr_index_rib1, curr_index_rib2, current_samples);




	clear_mp(cos_angle_CAB); clear_mp(cos_angle_BCA); clear_mp(cos_angle_ABC);
	clear_mp(cos_angle_DCA); clear_mp(cos_angle_ADC); clear_mp(cos_angle_CAD);


	clear_mp(length_AB); clear_mp(length_AC); clear_mp(length_BC); clear_mp(length_DA); clear_mp(length_DC);
	clear_mp(dot_CAB); clear_mp(dot_BCA); clear_mp(dot_ABC);
	clear_mp(dot_DCA); clear_mp(dot_ADC); clear_mp(dot_CAD);

	clear_mp(temp);
	// clean up at the end.  i wish scope was deletion!

	for (unsigned int ii=0; ii<num_vectors_needed; ii++) {
		clear_vec_mp(bulk_vectors[ii]);
	}
	free(bulk_vectors);



	return;
}




double compute_abs_of_difference_from_sixtydegrees(comp_mp temp, comp_mp length1, comp_mp length2, comp_mp dot_prod)
{
	comp_d cos;

	mul_mp(temp, length1, length2);
	div_mp(temp, dot_prod, temp);
	mp_to_d(cos, temp);
	double error = fabs(cos->r - 0.5);


	return error;
}




double compute_square_of_difference_from_sixtydegrees(comp_mp temp, comp_mp length1, comp_mp length2, comp_mp dot_prod)
{
	comp_d cos;

	mul_mp(temp, length1, length2);
	div_mp(temp, dot_prod, temp);
	mp_to_d(cos, temp);
	double error = fabs(cos->r - 0.5);


//	double angle = acos(cos->r);
//	double pi_over_three = acos(-1)/3;
//	double error = fabs(angle-pi_over_three);


//	div_mp(temp, length1,length2)
//	mp_to_d(cos,temp);
//	double error2 = fabs(cos->r - 1);
//
//	div_mp(temp, length2,length1)
//	mp_to_d(cos,temp);
//	double error3 = fabs(cos->r - 1);
//
//	double error = error1 + error2 + error3;

	return error*error;//*error*error
}





// if x> 0.5
// 	cycle_num = c2;
// 	pi_out = 1;
// 	pi_mid = 0.5;
// 	p = (x-0.5)*2;
// 	r = pi_out + (pi_mid - pi_out) * (1-p)^cycle_num;
// else
// 	cycle_num = c1;
// 	pi_out = 0;
// 	pi_mid = 0.5;
// 	p = x*2;
// 	r = pi_out + (pi_mid - pi_out) * (1-p)^cycle_num;
// end
void ScaleByCycleNum(comp_mp result, comp_mp input, int cycle_num_l, int cycle_num_r)
{
	comp_mp one;  init_mp(one); set_one_mp(one);
	comp_mp two;  init_mp(two); mpf_set_str(two->r, "2.0", 10); mpf_set_str(two->i, "0.0", 10);
	comp_mp half;  init_mp(half); mpf_set_str(half->r, "0.5", 10); mpf_set_str(half->i, "0.0", 10);

	comp_mp temp1, temp2; init_mp(temp1); init_mp(temp2);
	comp_mp p; init_mp(p);

	// pi_out + (pi_mid - pi_out) * (1-p)^cycle_num;
	if ( mpf_get_d(input->r)>=0.5 )
	{
		// p = (x-0.5)*2;
		sub_mp(temp1, input, half);
		mul_mp(p, temp1 ,two);

		sub_mp(temp1, one, p);
		exp_mp_int(temp2, temp1, cycle_num_r);
		mul_mp(temp1, half, temp2);
		// 1 - (0.5) * (1-p)^cycle_num_r;


		sub_mp(result, one, temp1);
	}
	else
	{
		// p = x*2
		sub_mp(temp1, half, input);
		mul_mp(p, temp1, two);
		sub_mp(temp1, one, p);
		exp_mp_int(temp2, temp1, cycle_num_l);
		mul_mp(result, half, temp2);
		// r = 0.5 * (1-p)^cycle_num
	}

	clear_mp(one); clear_mp(half); clear_mp(temp1); clear_mp(temp2); clear_mp(p); clear_mp(two);
}







void set_witness_set_mp(WitnessSet & W, vec_mp new_linear, vec_mp new_point)
{
	W.reset_points();
	W.add_point(new_point);

	W.reset_linears();
	W.add_linear(new_linear);
}








int get_dir_mptype_dimen(boost::filesystem::path & Dir_Name, int & MPType, int & dimension){

	
	std::ifstream fin("Dir_Name");

	if (!fin.is_open())
	{
		std::cout << color::red() << "did not find a decomposition in this directory.  currently uses file `Dir_Name` to record the name of the directory of the decomposition.  please ensure you have a completed decomposition in this directory, and `Dir_Name` is intact, specifying the correct location.\n" << color::console_default();
		dimension = -1;
		MPType = -1;
		Dir_Name = "missing";
		return MPType;
	}

	std::string line;
	std::getline(fin, line);

	fin >> MPType;
	fin >> dimension;

	Dir_Name = line;
	Dir_Name = Dir_Name.filename();
	return MPType;
}
