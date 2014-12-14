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
	
	
	int bcastme = PARSING;
	MPI_Bcast(&bcastme, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//end parser-bertini essentials
	
	
	parse_input(const_cast< char*> (filename.c_str()), &trackType, MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);
	
}





void prog_config::move_to_temp()
{
	if (!boost::filesystem::exists(this->working_dir())) {
		boost::filesystem::create_directory(this->working_dir());
	}
	
	if (!boost::filesystem::is_directory(this->working_dir())) {
		std::cerr << "trying to move into a directory which is a regular file!" << std::endl;
		// add error code here
	}
	
	chdir(this->working_dir().c_str());
	
	if (this->verbose_level()>=3)
		std::cout << "moved to working_dir '" << this->working_dir().string() << "'" << std::endl;
}

void prog_config::move_to_called()
{
	
	chdir(this->called_dir().c_str());
	
	if (this->verbose_level()>=3)
		std::cout << "moved to called_dir '" << this->called_dir().string() << "'" << std::endl;
}



int BertiniRealConfig::startup()
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
	IN = safe_fopen_read(this->input_filename());
	fclose(IN);
	
	
	
	if (this->user_projection()) {
		IN = safe_fopen_read(this->projection_filename());
		fclose(IN);
	}
	
	if (this->user_sphere()) {
		IN = safe_fopen_read(this->bounding_sphere_filename());
		fclose(IN);
	}
	
	return 0;
}

void BertiniRealConfig::splash_screen()
{
	printf("\n BertiniReal(TM) v%s\n\n", BERTINI_REAL_VERSION_STRING);
	printf(" D.A. Brake with\n D.J. Bates, W. Hao, J.D. Hauenstein,\n A.J. Sommese, C.W. Wampler\n\n");
	printf("(using GMP v%d.%d.%d, MPFR v%s)\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
	
	
	
	
	
	bertini_splash_screen();
	
	
	
}



void BertiniRealConfig::display_current_options()
{
	std::cout << "current options:\n\n";
	
	std::cout << "user_projection: " << user_projection();
	if (this->user_projection())
		std::cout << ", " << projection_filename_.string() << ",\n";
	else
		std::cout << "\n";
	
	

	
	
	std::cout << "user_sphere: " << user_sphere();
	if (user_sphere())
		std::cout << ", " << bounding_sphere_filename_.string() << "\n";
	else
		std::cout << "\n";
	
	
	std::cout << "input_filename: " << input_filename_.string() << "\n";

	
	std::cout << "stifle_text: " << stifle_text() << "\n";
	std::cout << "matlab_command: " << matlab_command() << "\n";
	std::cout << "output_directory base name: " << output_dir() << std::endl;
	
}




int  BertiniRealConfig::parse_commandline(int argc, char **argv)
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
			{"mode",required_argument,0,'M'}, {"m",required_argument,0,'M'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		
		choice = getopt_long_only (argc, argv, "d:c:Dg:V:o:smp:S:i:qvhM:", // if followed by colon, requires option.  two colons is optional
								   long_options, &option_index);
		
		/* Detect the end of the options. */
		if (choice == -1)
			break;
		
		switch (choice)
		{
			case 'd':
				target_dimension_ = atoi(optarg);
				break;
				
			case 'c':
				target_component_ = atoi(optarg);
				break;
				
			case 'D':
				this->debugwait(true);
				break;
		
			case 'g':
				this->use_gamma_trick_ = atoi(optarg);
				break;
				
			case 'V':
				verbose_level(atoi(optarg));
				break;

				
			case 'o':
				output_dir(boost::filesystem::absolute(optarg));
				break;
				
			case 's':
				stifle_text("\0");
				break;
				
			case 'm':
				merge_edges(false);
				break;
				
			case 'p':
				user_projection(true);
				projection_filename_ = boost::filesystem::absolute(optarg);
				break;
				
				
				
			case 'S':
				user_sphere(true);
				this->bounding_sphere_filename_ = boost::filesystem::absolute(optarg);
				break;
				
				
			case 'i': // input filename
				input_filename_ = optarg;
				break;
				
			case 'q':
				quick_run(1);
				break;
				
			case 'Q':
				quick_run(2);
				break;

			case 'v':
				printf("\n BertiniReal(TM) v %s\n\n", BERTINI_REAL_VERSION_STRING);
				exit(0);
				break;
				
			case 'h':
				
				printf("\nBertiniReal(TM) v %s.\n\n", BERTINI_REAL_VERSION_STRING);
				printf("Online at bertinireal.com\n\n");
				printf("For immediate support, send email to danielthebrake@gmail.com\n\n");
				BertiniRealConfig::print_usage();
				exit(0);
				break;
				
				
			case 'M':
			{
				std::string usermode = optarg;
				
				std::cout << usermode << std::endl;
				
				if (usermode=="bertini_real") {
					this->primary_mode_ = BERTINIREAL;
				}
				else if (usermode=="crit") {
					this->primary_mode_ = CRIT;
				}
				else
				{
					std::cout << "bad mode of operation.  acceptable options are [bertini_real] and crit." << std::endl;
					exit(0);
				}
				
				break;
			}
			case '?':
				/* getopt_long already printed an error message. */
				break;
				
			default:
				BertiniRealConfig::print_usage();
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
	
	
	
	
	this->set_called_dir(boost::filesystem::absolute(boost::filesystem::current_path()));
	this->output_dir(boost::filesystem::absolute(this->output_dir()));
	
	boost::filesystem::path new_name = this->called_dir();
	new_name/="temp";
	this->working_dir(new_name);

	
	
	
	return 0;
}


void BertiniRealConfig::print_usage()
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

void BertiniRealConfig::init()
{
	target_component_ = -2;
	target_dimension_ = -1;
	
	quick_run_ = 0;
	debugwait_ = false;
	max_deflations_ = 10;
	
	user_projection_ = false;
	projection_filename_ = "";
	
	orthogonal_projection_ = true;
	
	user_sphere_ = false;
	bounding_sphere_filename_ = "";
	
	input_filename_ = "input";
	
	
	output_dir(boost::filesystem::absolute("./output"));
	
	
	stifle_membership_screen_ = true;
	stifle_text_ = " > /dev/null ";
	
	matlab_command_ = "matlab -nosplash -nodesktop -nojvm";
	verbose_level(0); // default to 0
	
	
	use_gamma_trick_ = false;
	
	merge_edges_ = true;
	
	primary_mode_ = BERTINIREAL;
	
	return;
}















void sampler_configuration::splash_screen()
{
	printf("\n Sampler module for Bertini_real(TM) v%s\n\n", SAMPLER_VERSION_STRING);
	printf(" D.A. Brake, \n with D.J. Bates, W. Hao, \n J.D. Hauenstein, A.J. Sommese, and C.W. Wampler\n\n");
	printf("(using GMP v%d.%d.%d, MPFR v%s)\n\n",
				 __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
	printf("See the website at www.bertinireal.com\n\n");
	printf("Send email to danielthebrake@gmail.com for assistance.\n\n");
	
	
	
	bertini_splash_screen();
	
	
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
				this->verbose_level(atoi(optarg));
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
















