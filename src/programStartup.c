#include "programStartup.h"

void splash_screen(){
	printf("\n BertiniReal(TM) v%s\n\n", BERTINI_REAL_VERSION_STRING);
  printf(" D.J. Bates, D. Brake,\n W. Hao, J.D. Hauenstein,\n A.J. Sommese, C.W. Wampler\n\n");
  printf("(using GMP v%d.%d.%d, MPFR v%s)\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
	
}



void display_current_options(program_configuration options){
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
	
	mypause();
}


void print_usage(){
	printf("bertini_real has the following options:\n");
	printf("option name(s)\t\t\targument\n\n");
	printf("-p -projection \t\t\t'filename'\n");
	printf("-r -randomization \t\t'filename'\n");
	printf("-w -witness\t\t\t'filename'\n");
	printf("-i -input\t\t\t'filename'\n");
	printf("-v -version\t\t\t   -- \n");
	printf("-h -help\t\t\t   --");
	
	
	printf("\n\n\n");
	return;
}


/* Flag set by ‘--verbose’. */
static int verbose_flag;

int parse_options(int argc, char **argv, program_configuration *options){
	 // this code created based on gnu.org's description of getopt_long
	int choice;
	
	while (1)
	{
		static struct option long_options[] =
		{
			/* These options set a flag. */
			{"verbose", no_argument,       &verbose_flag, 1},
			{"brief",   no_argument,       &verbose_flag, 0},
			{"projection",		required_argument,			 0, 'p'},
			{"p",		required_argument,			 0, 'p'},
			{"randomization",		required_argument,			 0, 'r'},
			{"r",		required_argument,			 0, 'r'},
			{"input",		required_argument,			 0, 'i'},
			{"i",		required_argument,			 0, 'i'},
			{"witness",		required_argument,			 0, 'w'},
			{"w",		required_argument,			 0, 'w'},
			{"help",		no_argument,			 0, 'h'},
			{"h",		no_argument,			 0, 'h'},
			{"version",		no_argument,			 0, 'v'},
			{"v",		no_argument,			 0, 'v'},
//			{"allnew",	no_argument,			 &reuseflag, 0},
			/* These options don't set a flag.
			 We distinguish them by their indices. */
//			{"add",     no_argument,       0, 'a'},
//			{"append",  no_argument,       0, 'b'},
//			{"delete",  required_argument, 0, 'd'},
//			{"create",  required_argument, 0, 'choice'},
//			{"file",    required_argument, 0, 'f'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		
		choice = getopt_long_only (argc, argv, "r:p:v:i:w:h:v",
										 long_options, &option_index);
		
		/* Detect the end of the options. */
		if (choice == -1)
			break;
		
		switch (choice)
		{
			case 0:
				/* If this option set a flag, do nothing else now. */
				if (long_options[option_index].flag != 0)
					break;
				printf ("option %s", long_options[option_index].name);
				if (optarg)
					printf (" with arg %s", optarg);
				printf ("\n");
				break;
				
			case 'r':
				options->user_projection=1;
				options->projection_filename = optarg;
				break;
				
			case 'p':
				options->user_randomization=1;
				options->randomization_filename = optarg;
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
				print_usage();
				exit(0);
				break;
				
			case '?':
				/* getopt_long already printed an error message. */
				break;
				
			default:
				print_usage();
				exit(0);
		}
	}
	
	/* Instead of reporting ‘--verbose’
	 and ‘--brief’ as they are encountered,
	 we report the final status resulting from them. */
	if (verbose_flag)
		puts ("verbose flag is set");
	
	
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





int startup(program_configuration options)
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


void get_tracker_config(tracker_config_t *T,int MPType)
{

	//necessary for the setupConfig call
	double midpoint_tol, intrinsicCutoffMultiplier;
	int userHom = 0, useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, pathMod = 0, reducedOnly = 0, supersetOnly = 0, paramHom = 0;
	//end necessaries for the setupConfig call.
	
	
  setupConfig(T, &midpoint_tol, &userHom, &useRegen, &regenStartLevel, &maxCodim, &specificCodim, &pathMod, &intrinsicCutoffMultiplier, &reducedOnly, &supersetOnly, &paramHom, MPType);

	return;
}
