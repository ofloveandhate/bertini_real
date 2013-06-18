#include "programConfiguration.h"





void get_projection(vec_mp *pi,
										program_configuration program_options,
										solver_configuration solve_options,
										int num_vars,
										int num_projections)
{
	
	int ii,jj;
	for (ii=0; ii<num_projections; ii++) {
		change_size_vec_mp(pi[ii], num_vars);  pi[ii]->size = num_vars;
	}
	
	
	
	//assumes the vector pi is already initialized
	if (program_options.user_projection==1) {
		FILE *IN = safe_fopen_read(program_options.projection_filename); // we are already assured this file exists, but safe fopen anyway.
		int tmp_num_vars;
		fscanf(IN,"%d",&tmp_num_vars); scanRestOfLine(IN);
		if (tmp_num_vars!=num_vars-1) {
			printf("the number of variables appearing in the projection\nis not equal to the number of non-homogeneous variables in the problem\n");
			printf("please modify file to have %d coordinate pairs.\n",num_vars-1);
			abort();
		}
		
		for (ii=0; ii<num_projections; ii++) {
			set_zero_mp(&pi[ii]->coord[0]);
			for (jj=1; jj<num_vars; jj++) {
				mpf_inp_str(pi[ii]->coord[jj].r, IN, 10);
				mpf_inp_str(pi[ii]->coord[jj].i, IN, 10);
				scanRestOfLine(IN);
			}
		}
		fclose(IN);
	}
	else{
		for (ii=0; ii<num_projections; ii++) {
			set_zero_mp(&pi[ii]->coord[0]);
			for (jj=1; jj<num_vars; jj++) 
				get_comp_rand_real_mp(&pi[ii]->coord[jj]);
			
		}
	}
	
	return;
}









void parse_input_file(char filename[], int *MPType)
{
	
	unsigned int currentSeed;
	int trackType, genType = 0,  sharpenOnly, needToDiff, remove_temp, useParallelDiff = 0,userHom = 0;
  int my_id = 0, num_processes = 1, headnode = 0; // headnode is always 0
	
	
	
	//end parser-bertini essentials
	
	
	
	parse_input(filename, &trackType, MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);
	
}

void get_tracker_config(solver_configuration *solve_options,int MPType)
{

	//necessary for the setupConfig call
	double intrinsicCutoffMultiplier;
	int userHom = 0, useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, pathMod = 0, reducedOnly = 0, supersetOnly = 0, paramHom = 0;
	//end necessaries for the setupConfig call.
	
	
  setupConfig(&solve_options->T, &solve_options->midpoint_tol, &userHom, &useRegen, &regenStartLevel, &maxCodim, &specificCodim, &pathMod, &intrinsicCutoffMultiplier, &reducedOnly, &supersetOnly, &paramHom, MPType);

	return;
}







int BR_startup(program_configuration options)
/***************************************************************\
 * USAGE:    prepares the variables inputname and startname
 *      for use later in the program
 * ARGUMENTS:                                                    *
 * RETURN VALUES:           an integer                           *
 * NOTES:                                                        *
 \***************************************************************/
{
	
	
	
	// check for write privilege
  if (checkWritePrivilege())
  {
    printf("ERROR: BertiniReal does not have write privileges!\n");
    bexit(ERROR_WRITE_PRIVILEGE);
  }
	
	
	//test for presence of necessary files
	FILE *IN;
	IN = safe_fopen_read(options.input_filename);
	fclose(IN);
	
	IN = safe_fopen_read(options.witness_set_filename);
	fclose(IN);
	
	if (options.user_randomization) {
		IN = safe_fopen_read(options.randomization_filename);
		fclose(IN);
	}
	
	if (options.user_projection) {
		IN = safe_fopen_read(options.projection_filename);
		fclose(IN);
	}
	
	return 0;
}

void BR_splash_screen(){
	printf("\n BertiniReal(TM) v%s\n\n", BERTINI_REAL_VERSION_STRING);
  printf(" D.J. Bates, D. Brake,\n W. Hao, J.D. Hauenstein,\n A.J. Sommese, C.W. Wampler\n\n");
  printf("(using GMP v%d.%d.%d, MPFR v%s)\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
	
}



void BR_display_current_options(program_configuration options){
	printf("current options:\n\n");
	
	printf("user_projection: %d",options.user_projection);
	if (options.user_projection)
		printf(", %s\n",options.projection_filename);
	else
		printf("\n");
	
	
	printf("user_randomization: %d",options.user_randomization);
	if (options.user_randomization)
		printf(", %s\n",options.randomization_filename);
	else
		printf("\n");
	
	printf("input_filename: %s\n",options.input_filename);
	printf("witness_set_filename: %s\n",options.witness_set_filename);
	
	printf("stifle_text: %s\n",options.stifle_text);
	mypause();
}


void BR_print_usage(){
	printf("bertini_real has the following options:\n");
	printf("option name(s)\t\t\targument\n\n");
	printf("-p -pi -projection \t\t\t'filename'\n");
	printf("-r -randomization \t\t'filename'\n");
	printf("-w -witness\t\t\t'filename'\n");
	printf("-i -input\t\t\t'filename'\n");
	printf("-ns -nostifle\t\t\t   --\n");
	printf("-v -version\t\t\t   -- \n");
	printf("-h -help\t\t\t   --");
	printf("-box -b\t\t\t   bool");
	
	printf("\n\n\n");
	return;
}



int  BR_parse_commandline(int argc, char **argv, program_configuration *options){
	// this code created based on gnu.org's description of getopt_long
	int choice;
	
	while (1)
	{
		static struct option long_options[] =
		{
			/* These options set a flag. */
			{"nostifle", no_argument,       0, 's'}, {"ns", no_argument,					0, 's'},
			{"projection",		required_argument,			 0, 'p'}, {"p",		required_argument,			 0, 'p'}, {"pi",		required_argument,			 0, 'p'},
			{"randomization",		required_argument,			 0, 'r'}, {"r",		required_argument,			 0, 'r'},
			{"input",		required_argument,			 0, 'i'}, {"i",		required_argument,			 0, 'i'},
			{"witness",		required_argument,			 0, 'w'}, {"w",		required_argument,			 0, 'w'},
			{"help",		no_argument,			 0, 'h'}, {"h",		no_argument,			 0, 'h'},
			{"version",		no_argument,			 0, 'v'}, {"v",		no_argument,			 0, 'v'},
			{"output",		required_argument,			 0, 'o'}, {"out",		required_argument,			 0, 'o'}, {"o",		required_argument,			 0, 'o'},
			{"verb",		required_argument,			 0, 'V'},
			{"box",		required_argument,			 0, 'b'}, {"b",		required_argument,			 0, 'b'},
			{"gammatrick",		required_argument,			 0, 'g'}, {"g",		required_argument,			 0, 'g'},
			{"detjac",		no_argument,			 0, 'd'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		
		choice = getopt_long_only (argc, argv, "r:p:v:i:w:h:v:s",
															 long_options, &option_index);
		
		/* Detect the end of the options. */
		if (choice == -1)
			break;
		
		switch (choice)
		{
			case 'd':
				options->crit_solver = LINPRODTODETJAC;
				break;
			case 'g':
				options->use_gamma_trick = atoi(optarg);
				if (! (options->use_gamma_trick==0 || options->use_gamma_trick==1) ) {
					printf("value for 'gammatrick' or 'g' must be 1 or 0\n");
					exit(689);
				}
				break;
				
				
			case 'b':
				options->use_bounding_box = atoi(optarg);
				if (! (options->use_bounding_box==0 || options->use_bounding_box==1) ) {
					printf("value for 'box' must be 1 or 0\n");
					exit(689);
				}
				
			case 'V':
				options->verbose_level = atoi(optarg);
				break;
				
			case 'o':
				options->output_basename = optarg;
				break;
				
			case 's':
				options->stifle_text = "\0";
				break;
				
			case 'r':
				options->user_randomization=1;
				options->randomization_filename = optarg;
				break;
				
				
			case 'p':
				options->user_projection=1;
				options->projection_filename = optarg;
				break;
				
				
			case 'w': // witness_set_filename
				options->witness_set_filename = optarg;
				break;
				
			case 'i': // input filename
				options->input_filename = optarg;
				break;
				
			case 'v':
				printf("\n BertiniReal(TM) v %s\n\n", BERTINI_REAL_VERSION_STRING);
				exit(0);
				break;
				
			case 'h':
				printf("\nThis is BertiniReal v %s, developed by\nDan J. Bates, Daniel Brake,\nWenrui Hao, Jonathan D. Hauenstein,\nAndrew J. Sommmese, and Charles W. Wampler.\n\n", BERTINI_REAL_VERSION_STRING);
				printf("Send email to brake@math.colostate.edu for details about BertiniReal.\n\n");
				BR_print_usage();
				exit(0);
				break;
				
			case '?':
				/* getopt_long already printed an error message. */
				break;
				
			default:
				BR_print_usage();
				exit(0);
		}
	}
	
	
	
	
	/* Print any remaining command line arguments (not options). */
	if (optind < argc)
	{
		printf ("these options not processed: ");
		while (optind < argc)
			printf ("%s ", argv[optind++]);
		putchar ('\n');
	}
	
	return 0;
}


void BR_init_config(program_configuration *options) {
	
	
	options->user_projection = 0;
	options->projection_filename = (char *)bmalloc(MAX_STRLEN*sizeof(char));
	options->projection_filename = "";
	
	options->user_randomization = 0;
	options->randomization_filename = (char *)bmalloc(MAX_STRLEN*sizeof(char));
	options->randomization_filename = "";
	
	options->input_filename = (char *)bmalloc(MAX_STRLEN*sizeof(char));
	options->input_filename = "input\0";
	
	options->witness_set_filename = (char *)bmalloc(MAX_STRLEN*sizeof(char));
	options->witness_set_filename = "witness_set\0";
	
	options->output_basename = (char *)bmalloc(MAX_STRLEN*sizeof(char));
	options->output_basename = "output";
	
	options->input_deflated_filename = (char *)bmalloc(MAX_STRLEN*sizeof(char));
	//this must be set after deflation is run.
	
	options->crit_solver = NULLSPACE;
	
	options->stifle_membership_screen = 1;
	options->stifle_text = (char *)bmalloc(MAX_STRLEN*sizeof(char));
	options->stifle_text = " > /dev/null ";
	
	options->verbose_level = 0; // default to 0
	
	options->MPType = 2;
	
	options->use_bounding_box = 0;
	options->use_gamma_trick = 0;
	return;
}

void BR_clear_config(program_configuration *options) {
	
	if (options->user_projection) {
		free(options->projection_filename);
	}
	
	if (options->user_randomization) {
		free(options->randomization_filename);
	}
	
	free(options->input_filename);
	free(options->input_deflated_filename);
	
	return;
}






void sampler_splash_screen(){
	printf("\n Sampler module for BertiniReal(TM) v%s\n\n", BERTINI_REAL_VERSION_STRING);
	printf(" D.J. Bates, D. Brake,\n W. Hao, J.D. Hauenstein,\n A.J. Sommese, C.W. Wampler\n\n");
	printf("(using GMP v%d.%d.%d, MPFR v%s)\n\n",
				 __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
	
	printf("Send email to brake@math.colostate.edu for details.\n\n");
}


void sampler_print_usage(){
	printf("bertini_real has the following options:\n");
	printf("option name(s)\t\t\targument\n\n");
	printf("-ns -nostifle\t\t\t   --\n");
	printf("-v -version\t\t\t   -- \n");
	printf("-h -help\t\t\t   --\n");
	printf("-t -tol -tolerance \t\tdouble > 0\n");
	printf("-verb\t\t\t\tint\n");
	printf("-maxits -m \t\t\tint\n");
	printf("-gammatrick -g \t\t\tbool\n");
	printf("\n\n\n");
	return;
}

int  sampler_parse_commandline(int argc, char **argv, sampler_configuration *options){
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
			{"m",		required_argument,			 0, 'm'},
			{"maxits",		required_argument,			 0, 'm'},
			{"gammatrick",		required_argument,			 0, 'g'},
			{"g",		required_argument,			 0, 'g'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		
		choice = getopt_long_only (argc, argv, "h:v:s:ve",
															 long_options, &option_index);
		
		/* Detect the end of the options. */
		if (choice == -1)
			break;
		
		switch (choice)
		{
			case 's':
				options->stifle_text = "\0";
				break;
				
			case 'v':
				printf("\n Sampler module for BertiniReal(TM) v %s\n\n", BERTINI_REAL_VERSION_STRING);
				exit(0);
				break;
				
			case 't':
				
				mpf_set_str(options->TOL,optarg,10);
				break;
				
			case 'g':
				options->use_gamma_trick = atoi(optarg);
				if (! (options->use_gamma_trick==0 || options->use_gamma_trick==1) ) {
					printf("value for 'gammatrick' or 'g' must be 1 or 0\n");
					exit(689);
				}
				break;
				
			case 'V':
				options->verbose_level = atoi(optarg);
				break;
				
			case 'm':
				options->maximum_num_iterations = atoi(optarg);
				break;
				
			case 'h':
				
				sampler_print_usage();
				exit(0);
				break;
				
			case '?':
				/* getopt_long already printed an error message. */
				break;
				
			default:
				sampler_print_usage();
				exit(0);
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



void sampler_init_config(sampler_configuration *options) {
	options->stifle_membership_screen = 1;
	options->stifle_text = (char *)bmalloc(MAX_STRLEN*sizeof(char));
	options->stifle_text = " > /dev/null ";
	
	options->verbose_level = 0; // default to 0
	
	options->maximum_num_iterations = 10;
	
	mpf_init(options->TOL);
	mpf_set_d(options->TOL, 1e-1); // this should be made adaptive to the span of the projection values or the endpoints

	options->use_gamma_trick = 0;
	
	return;
}

void sampler_clear_config(sampler_configuration *options) {
	
	mpf_clear(options->TOL);
	return;
}








void solver_init_config(solver_configuration *options){
	options->allow_multiplicity = 0;
	options->allow_singular = 0;
	options->allow_infinite = 0;
	options->allow_unsuccess = 0;
	options->use_midpoint_checker = 1;
	options->show_status_summary = 0;
	options->verbose_level = 0; // default to 0.  higher is more verbose
	options->use_gamma_trick = 0;
	
	options->complete_witness_set = 1;
}


void solver_clear_config(solver_configuration *options){
	//has no fields which require clearing.
	return;
}




