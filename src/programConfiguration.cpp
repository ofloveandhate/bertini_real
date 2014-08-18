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
	printf(" Daniel A. Brake with\n Dan J. Bates, Wenrui Hao, Jonathon D. Hauenstein,\n Andrew J. Sommese, Charles W. Wampler\n\n");
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
			{"debug", no_argument, 0, 'D'},
			{"dim",required_argument,0,'d'}, {"d",required_argument,0,'d'},
			{"component",required_argument,0,'c'}, {"comp",required_argument,0,'c'}, {"c",required_argument,0,'c'},
			{"gammatrick",required_argument,0, 'g'}, {"g",required_argument, 0, 'g'},
			{"verb",	required_argument,0, 'V'},
			{"output",	required_argument,0, 'o'}, {"out",required_argument, 0, 'o'}, {"o",	required_argument, 0, 'o'},
			{"nostifle", no_argument,       0, 's'}, {"ns", no_argument, 0, 's'},
			{"nomerge",no_argument,0,'m'}, {"nm",no_argument,0,'m'},
			{"projection",required_argument,0, 'p'}, {"p",required_argument,0, 'p'}, {"pi",	required_argument,0,'p'},
			{"sphere",required_argument, 0, 'S'}, {"s",required_argument, 0, 'S'},
			{"input",required_argument,	0, 'i'}, {"i",required_argument, 0, 'i'},
			{"quick",no_argument,0,'q'}, {"q",no_argument,0,'Q'},
			{"veryquick",no_argument,0,'q'}, {"vq",no_argument,0,'Q'},
			{"version",		no_argument,			 0, 'v'}, {"v",		no_argument,			 0, 'v'},
			{"help",		no_argument,			 0, 'h'}, {"h",		no_argument,			 0, 'h'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		
		choice = getopt_long_only (argc, argv, "d:c:Dg:V:o:smp:S:i:qvh", // if followed by colon, requires option.  two colons is optional
								   long_options, &option_index);
		
		/* Detect the end of the options. */
		if (choice == -1)
			break;
		
		switch (choice)
		{
			case 'd':
				target_dimension = atoi(optarg);
				break;
				
			case 'c':
				target_component = atoi(optarg);
				break;
				
			case 'D':
				this->debugwait = 1;
				break;
		
			case 'g':
				this->use_gamma_trick = atoi(optarg);
				if (! (this->use_gamma_trick==0 || this->use_gamma_trick==1) ) {
					printf("value for 'gammatrick' or 'g' must be 1 or 0\n");
					br_exit(689);
				}
				break;
				
			case 'V':
				verbose_level = atoi(optarg);
				break;

				
			case 'o':
				output_dir = boost::filesystem::absolute(optarg);
				break;
				
			case 's':
				stifle_text = "\0";
				break;
				
			case 'm':
				merge_edges = false;
				break;
				
			case 'p':
				user_projection=1;
				projection_filename = optarg;
				break;
				
				
				
			case 'S':
				user_sphere = true;
				this->bounding_sphere_filename = boost::filesystem::absolute(optarg);
				break;
				
				
			case 'i': // input filename
				input_filename = optarg;
				break;
				
			case 'q':
				quick_run = 1;
				break;
				
			case 'Q':
				quick_run = 2;
				break;

			case 'v':
				printf("\n BertiniReal(TM) v %s\n\n", BERTINI_REAL_VERSION_STRING);
				exit(0);
				break;
				
			case 'h':
				
				printf("\nBertiniReal(TM) v %s.\n\n", BERTINI_REAL_VERSION_STRING);
				printf("Online at bertinireal.com\n\n");
				printf("For immediate support, send email to danielthebrake@gmail.com\n\n");
				BR_configuration::print_usage();
				exit(0);
				break;
				
			case '?':
				/* getopt_long already printed an error message. */
				break;
				
			default:
				BR_configuration::print_usage();
				exit(0); //
		}
	}
	
	
    /* Print any remaining command line arguments (not options). */
	if (optind < argc)
	{
		printf ("these options were not processed: ");
		while (optind < argc)
			printf ("%s ", argv[optind++]);
		printf("\n");
		br_exit(-2363);
	}
	
	
	
	
	this->called_dir = boost::filesystem::absolute(boost::filesystem::current_path());
	this->output_dir = boost::filesystem::absolute(this->output_dir);
	this->working_dir = this->called_dir;
	this->working_dir/="temp";
	
	
	
	return 0;
}


void BR_configuration::print_usage()
{
	printf("bertini_real has the following options:\n----------------------\n");
	printf("option name(s)\t\t\targument\n\n");
	printf("-p -pi -projection \t\t\t'filename'\n");
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
	target_component = -2;
	target_dimension = -1;
	
	quick_run = 0;
	debugwait = 0;
	max_deflations = 10;
	
	user_projection = 0;
	projection_filename = "";
	
	
	
	user_sphere = false;
	bounding_sphere_filename = "";
	
	input_filename = "input";
	
	
	output_dir = boost::filesystem::absolute("output");
	
	
	stifle_membership_screen = 1;
	stifle_text = " > /dev/null ";
	
	bertini_command = "~/bin/bertini_serial";
	matlab_command = "matlab -nosplash -nodesktop -nojvm";
	verbose_level = 0; // default to 0
	
	MPType = 2;
	
	use_gamma_trick = 0;
	
	merge_edges = true;
	
	return;
}















void sampler_configuration::splash_screen()
{
	printf("\n Sampler module for Bertini_real(TM) v%s\n\n", SAMPLER_VERSION_STRING);
	printf(" D. Brake, \n with D.J. Bates, W. Hao, \n J.D. Hauenstein, A.J. Sommese, and C.W. Wampler\n\n");
	printf("(using GMP v%d.%d.%d, MPFR v%s)\n\n",
				 __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
	printf("See the website at www.bertinireal.com\n\n");
	printf("Send email to danielthebrake@gmail.com for assistance.\n\n");
}


void sampler_configuration::print_usage()
{
	printf("Bertini_real has the following options:\n");
	printf("option name(s)\t\t\targument\n\n");
	printf("-ns -nostifle\t\t\t   --\n");
	printf("-v -version\t\t\t   -- \n");
	printf("-h -help\t\t\t   --\n");
	printf("-t -tol -tolerance \t\tdouble > 0\n");
	printf("-verb\t\t\t\tint\n");
	printf("-maxits -m \t\t\tint\n");
	printf("-gammatrick -g \t\t\tbool\n");
	printf("-fixed \t\t\tint number samples per edge\n");
	printf("-projbin -pb \t\t\t -- turn on projection-based binning rather than distance based.\n");
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
			{"nd", no_argument,0,'d'},
			{"projbin",		no_argument,			 0, 'b'},
			{"pb", no_argument,0,'b'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		
		choice = getopt_long_only (argc, argv, "bdf:svt:g:V:m:h",
															 long_options, &option_index);
		
		/* Detect the end of the options. */
		if (choice == -1)
			break;
		
		switch (choice)
		{
			case 'b':
				
				use_projection_binning = true;
				
				break;
				
			case 'd':
				no_duplicates = false;
				break;
				
				
			case 'f':
				use_fixed_sampler = true;
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
				printf("\n Sampler module for Bertini_real(TM) version %s\n\n", SAMPLER_VERSION_STRING);
				exit(0);
				break;
				
			case 't':
				
				mpf_set_str(this->TOL,optarg,10);
				break;
				
			case 'g':
				this->use_gamma_trick = atoi(optarg);
				if (! (this->use_gamma_trick==0 || this->use_gamma_trick==1) ) {
					printf("value for 'gammatrick' or 'g' must be 1 or 0\n");
					exit(0);
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
















