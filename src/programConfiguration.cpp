#include "programConfiguration.hpp"






void parse_preproc_data(boost::filesystem::path filename,  preproc_data *PPD)
{

	setupPreProcData(const_cast< char*> (filename.c_str()), PPD);

}



void parse_input_file(boost::filesystem::path filename)
{
	
	unsigned int currentSeed;
	int trackType, genType = 0,  sharpenOnly, needToDiff, remove_temp, useParallelDiff = 0,userHom = 0;
	int my_id = 0, num_processes = 1, headnode = 0; // headnode is always 0
	int MPType;
	
	int bcastme = PARSING;
	MPI_Bcast(&bcastme, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//end parser-bertini essentials
	

	parse_input(const_cast< char*> (filename.c_str()), &trackType, &MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);
	
}




void parse_input_file(boost::filesystem::path filename, int * MPType)
{
	
	unsigned int currentSeed;
	int trackType, genType = 0,  sharpenOnly, needToDiff, remove_temp, useParallelDiff = 0,userHom = 0;
	int my_id = 0, num_processes = 1, headnode = 0; // headnode is always 0
	
	
	//end parser-bertini essentials
	
	
	parse_input(const_cast< char*> (filename.c_str()), &trackType, MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);
	
}





void prog_config::move_to_temp()
{
	if (!boost::filesystem::exists(this->working_dir)) {
		boost::filesystem::create_directory(this->working_dir);
	}
	
	if (!boost::filesystem::is_directory(this->working_dir)) {
		std::cerr << "trying to move into a directory which is a regular file!" << std::endl;
		// add error code here
	}
	
	chdir(this->working_dir.c_str());
	
	if (this->verbose_level>=3)
		std::cout << "moved to working_dir '" << this->working_dir.string() << "'" << std::endl;
}

void prog_config::move_to_called()
{
	
	chdir(this->called_dir.c_str());
	
	if (this->verbose_level>=3)
		std::cout << "moved to called_dir '" << this->called_dir.string() << "'" << std::endl;
}



int BR_configuration::startup()
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
	IN = safe_fopen_read(this->input_filename.c_str());
	fclose(IN);
	
	IN = safe_fopen_read(this->witness_set_filename.c_str());
	fclose(IN);
	
	
	if (this->user_projection) {
		IN = safe_fopen_read(this->projection_filename.c_str());
		fclose(IN);
	}
	
	if (this->user_sphere) {
		IN = safe_fopen_read(this->bounding_sphere_filename.c_str());
		fclose(IN);
	}
	
	return 0;
}

void BR_configuration::splash_screen()
{
	printf("\n BertiniReal(TM) v%s\n\n", BERTINI_REAL_VERSION_STRING);
	printf(" Daniel A Brake with\n D.J. Bates,\n W. Hao, J.D. Hauenstein,\n A.J. Sommese, C.W. Wampler\n\n");
	printf("(using GMP v%d.%d.%d, MPFR v%s)\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
	
}



void BR_configuration::display_current_options()
{
	printf("current options:\n\n");
	
	printf("user_projection: %d",this->user_projection);
	if (this->user_projection)
		printf(", %s\n",this->projection_filename.c_str());
	else
		printf("\n");
	
	

	
	
	printf("user_sphere: %d",user_sphere);
	if (user_sphere)
		printf(", %s\n",bounding_sphere_filename.c_str());
	else
		printf("\n");
	
	
	printf("input_filename: %s\n",this->input_filename.c_str());
	printf("witness_set_filename: %s\n",this->witness_set_filename.c_str());
	
	std::cout << "stifle_text: " << this->stifle_text << std::endl;
	std::cout << "bertini_command: " << this->bertini_command << std::endl;
	std::cout << "output_directory base name: " << this->output_dir << std::endl;
	
}




int  BR_configuration::parse_commandline(int argc, char **argv)
{
	// this code created based on gnu.org's description of getopt_long
	int choice;
	
	while (1)
	{
		static struct option long_options[] =
		{
			/* These options set a flag. */
			{"nostifle", no_argument,       0, 's'}, {"ns", no_argument,					0, 's'},
			{"projection",		required_argument,			 0, 'p'}, {"p",		required_argument,			 0, 'p'}, {"pi",		required_argument,			 0, 'p'},
			{"input",		required_argument,			 0, 'i'}, {"i",		required_argument,			 0, 'i'},
			{"witness",		required_argument,			 0, 'w'}, {"w",		required_argument,			 0, 'w'},
			{"help",		no_argument,			 0, 'h'}, {"h",		no_argument,			 0, 'h'},
			{"version",		no_argument,			 0, 'v'}, {"v",		no_argument,			 0, 'v'},
			{"output",		required_argument,			 0, 'o'}, {"out",		required_argument,			 0, 'o'}, {"o",		required_argument,			 0, 'o'},
			{"verb",		required_argument,			 0, 'V'},
			{"sphere",		required_argument,			 0, 'S'}, {"s",		required_argument,			 0, 'S'},
			{"gammatrick",		required_argument,			 0, 'g'}, {"g",		required_argument,			 0, 'g'},
			{"detjac",		no_argument,			 0, 'd'},
			{"debug", no_argument, 0, 'D'},
			{"quick",no_argument,0,'q'},{"q",no_argument,0,'q'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		
		choice = getopt_long_only (argc, argv, "D:d:g:b:V:o:s:r:p:w:i:v:h",
								   long_options, &option_index);
		
		/* Detect the end of the options. */
		if (choice == -1)
			break;
		
		switch (choice)
		{
			case 'q':
				this->quick_run = true;
				break;

				
			case 'D':
				this->debugwait = 1;
				break;
				
			case 'd':
				this->crit_solver = LINPRODTODETJAC;
				break;
			case 'g':
				this->use_gamma_trick = atoi(optarg);
				if (! (this->use_gamma_trick==0 || this->use_gamma_trick==1) ) {
					printf("value for 'gammatrick' or 'g' must be 1 or 0\n");
					exit(689);
				}
				break;
				
				
			case 'S':
				user_sphere = true;
				this->bounding_sphere_filename = boost::filesystem::absolute(optarg);
				break;
				
			case 'V':
				this->verbose_level = atoi(optarg);
				break;
				
			case 'o':
				this->output_dir = boost::filesystem::absolute(optarg);
				break;
				
			case 's':
				this->stifle_text = "\0";
				break;
				
				
				
			case 'p':
				this->user_projection=1;
				this->projection_filename = optarg;
				break;
				
				
			case 'w': // witness_set_filename
				this->witness_set_filename = optarg;
				break;
				
			case 'i': // input filename
				this->input_filename = optarg;
				break;
				
			case 'v':
				printf("\n BertiniReal(TM) v %s\n\n", BERTINI_REAL_VERSION_STRING);
				exit(0);
				break;
				
			case 'h':
				
				splash_screen();
				
				printf("\nThis is BertiniReal v %s, developed by\nDaniel A. Brake with Dan J. Bates,\nWenrui Hao, Jonathan D. Hauenstein,\nAndrew J. Sommmese, and Charles W. Wampler.\n\n", BERTINI_REAL_VERSION_STRING);
				printf("Send email to brake@math.colostate.edu for details about BertiniReal.\n\n");
				BR_configuration::print_usage();
				exit(0);
				break;
				
			case '?':
				/* getopt_long already printed an error message. */
				break;
				
			default:
				BR_configuration::print_usage();
				exit(0); // this indices a memory leak.
		}
	}
	
	
    
	this->called_dir = boost::filesystem::absolute(boost::filesystem::current_path());
	this->output_dir = boost::filesystem::absolute(this->output_dir);
	this->working_dir = this->called_dir;
	this->working_dir/="temp";
	
	/* Print any remaining command line arguments (not options). */
	if (optind < argc)
	{
		printf ("these options were not processed: ");
		while (optind < argc)
			printf ("%s ", argv[optind++]);
		putchar ('\n');
	}
	
	return 0;
}


void BR_configuration::print_usage()
{
//	
//	{"nostifle", no_argument,       0, 's'}, {"ns", no_argument,					0, 's'},
//	{"projection",		required_argument,			 0, 'p'},
//	{"p",		required_argument,			 0, 'p'}, {"pi",		required_argument,			 0, 'p'},
//	{"r",		required_argument,			 0, 'r'},
//	{"input",		required_argument,			 0, 'i'}, {"i",		required_argument,			 0, 'i'},
//	{"witness",		required_argument,			 0, 'w'}, {"w",		required_argument,			 0, 'w'},
//	{"help",		no_argument,			 0, 'h'}, {"h",		no_argument,			 0, 'h'},
//	{"version",		no_argument,			 0, 'v'}, {"v",		no_argument,			 0, 'v'},
//	{"output",		required_argument,			 0, 'o'}, {"out",		required_argument,			 0, 'o'}, {"o",		required_argument,			 0, 'o'},
//	{"verb",		required_argument,			 0, 'V'},
//	{"box",		required_argument,			 0, 'b'}, {"b",		required_argument,			 0, 'b'},
//	{"gammatrick",		required_argument,			 0, 'g'}, {"g",		required_argument,			 0, 'g'},
//	{"detjac",		no_argument,			 0, 'd'},
//	{"debug", no_argument, 0, 'D'},
//	{"nomerge", no_argument, 0, 'm'},
//	{"quick",no_argument,0,'q'},{"q",no_argument,0,'q'},
	
	printf("bertini_real has the following options:\n----------------------\n");
	printf("option name(s)\t\t\targument\n\n");
	printf("-p -pi -projection \t\t\t'filename'\n");
	printf("-w -witness\t\t\t'filename'\n");
	printf("-i -input\t\t\t'filename'\n");
	printf("-ns -nostifle\t\t\t   --\n");
	printf("-v -version\t\t\t   -- \n");
	printf("-h -help\t\t\t   --\n");
	printf("-sphere -b\t\t\t   'filename'\n");
	printf("-q -quick\t\t\t --\n");
	printf("-debug\t\t\t --\n");
	printf("-gammatrick\t\t\t bool\n");
	printf("\n\n\n");
	return;
}

void BR_configuration::init()
{
	

	quick_run = false;
	this->debugwait = 0;
	this->max_deflations = 10;
	
	this->user_projection = 0;
	this->projection_filename = "";
	
	
	
	user_sphere = false;
	bounding_sphere_filename = "";
	
	this->input_filename = "input";
	
	this->witness_set_filename = "witness_set";
	
	this->output_dir = boost::filesystem::absolute("output");
	
	
	this->crit_solver = NULLSPACE;
	
	this->stifle_membership_screen = 1;
	this->stifle_text = " > /dev/null ";
	
	this->bertini_command = "~/bin/bertini_serial";
	this->matlab_command = "matlab -nosplash -nodesktop -nojvm";
	this->verbose_level = 0; // default to 0
	
	this->MPType = 2;
	
	this->use_gamma_trick = 0;
	
	merge_edges = true;
	
	return;
}















void sampler_configuration::splash_screen()
{
	printf("\n Sampler module for BertiniReal(TM) v%s\n\n", BERTINI_REAL_VERSION_STRING);
	printf(" D.J. Bates, D. Brake,\n W. Hao, J.D. Hauenstein,\n A.J. Sommese, C.W. Wampler\n\n");
	printf("(using GMP v%d.%d.%d, MPFR v%s)\n\n",
				 __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
	
	printf("Send email to brake@math.colostate.edu for details.\n\n");
}


void sampler_configuration::print_usage()
{
	printf("bertini_real has the following options:\n");
	printf("option name(s)\t\t\targument\n\n");
	printf("-ns -nostifle\t\t\t   --\n");
	printf("-v -version\t\t\t   -- \n");
	printf("-h -help\t\t\t   --\n");
	printf("-t -tol -tolerance \t\tdouble > 0\n");
	printf("-verb\t\t\t\tint\n");
	printf("-maxits -m \t\t\tint\n");
	printf("-gammatrick -g \t\t\tbool\n");
	printf("-fixed \t\t\tint number samples per edge\n");
	printf("\n\n\n");
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
			{"m",		required_argument,			 0, 'm'},
			{"maxits",		required_argument,			 0, 'm'},
			{"gammatrick",		required_argument,			 0, 'g'},
			{"g",		required_argument,			 0, 'g'},
			{"fixed",		required_argument,			 0, 'f'},
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
			case 'f':
				use_fixed_sampler = true;
				target_num_samples = atoi(optarg);
				
				break;
				
			case 's':
				this->stifle_text = "\0";
				break;
				
			case 'v':
				printf("\n Sampler module for BertiniReal(TM) v %s\n\n", BERTINI_REAL_VERSION_STRING);
				exit(0);
				break;
				
			case 't':
				
				mpf_set_str(this->TOL,optarg,10);
				break;
				
			case 'g':
				this->use_gamma_trick = atoi(optarg);
				if (! (this->use_gamma_trick==0 || this->use_gamma_trick==1) ) {
					printf("value for 'gammatrick' or 'g' must be 1 or 0\n");
					exit(689);
				}
				break;
				
			case 'V':
				this->verbose_level = atoi(optarg);
				break;
				
			case 'm':
				this->maximum_num_iterations = atoi(optarg);
				break;
				
			case 'h':
				
				sampler_configuration::print_usage();
				exit(0);
				break;
				
			case '?':
				/* getopt_long already printed an error message. */
				break;
				
			default:
				sampler_configuration::print_usage();
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
















