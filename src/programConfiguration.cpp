#include "programConfiguration.hpp"
#include <iomanip>





void parse_preproc_data(boost::filesystem::path filename,  preproc_data *PPD)
{

	setupPreProcData(const_cast< char*> (filename.c_str()), PPD);

}



unsigned int parse_input_file(boost::filesystem::path filename)
{
	int MPType;
	return parse_input_file(filename, &MPType);
}




unsigned int parse_input_file(boost::filesystem::path filename, int * MPType)
{

	unsigned int currentSeed;
	int trackType, genType = 0,  sharpenOnly, needToDiff, remove_temp, useParallelDiff = 0,userHom = 0;
	int my_id = 0, num_processes = 1, headnode = 0; // headnode is always 0


	int bcastme = PARSING;
	MPI_Bcast(&bcastme, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//end parser-bertini essentials

	parse_input(const_cast< char*> (filename.c_str()), &trackType, MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);
	return currentSeed;
}



void ParallelismConfig::init()
{

	force_no_parallel_ = false;
	numprocs_ = 1;
	headnode_ = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &this->numprocs_);
	MPI_Comm_rank(MPI_COMM_WORLD, &this->my_id_);


	if (is_head())
		worker_level_ = 0;
	else
		worker_level_ = 1;



	my_communicator_ = MPI_COMM_WORLD; // default communicator is MPI_COMM_WORLD




	//		MPI_Group orig_group, new_group;
	//
	//		/* Extract the original group handle */
	//
	//		MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
	//
	//		/* Create new communicator and then perform collective communications */
	//
	//		MPI_Comm_create(MPI_COMM_WORLD, orig_group, &my_communicator);
	//
	//		MPI_Comm_size(my_communicator, &this->numprocs);
	//		MPI_Comm_rank(my_communicator, &this->my_id);
	return;
}


void ParallelismConfig::init_active_workers()
{
	if (!available_workers_.empty()) {
		throw std::logic_error("set of available workers is not empty at call of init_active_workers...  it must be empty.  some previous process did not finish properly, dismissing all workers at the end.");
	}

	while (!available_workers_.empty())
		available_workers_.pop();


	for (int ii=1; ii<this->numprocs_; ii++) {
		available_workers_.push(ii);
		worker_status_[ii] = INACTIVE;
	}


}


int ParallelismConfig::activate_next_worker()
{
	int worker_id = available_workers_.front();
	available_workers_.pop();

	if (worker_status_[worker_id] == ACTIVE) {
		std::cout << "master tried making worker" << worker_id << " active when it was already active" << std::endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	worker_status_[worker_id] = ACTIVE;


	return worker_id;
}


void ParallelismConfig::deactivate(int worker_id)
{
	if (worker_status_[worker_id] == INACTIVE) {
		std::cout << "master tried decativating worker" << worker_id << " when it was already inactive" << std::endl;
		MPI_Abort(MPI_COMM_WORLD,2);
	}
	worker_status_[worker_id] = INACTIVE;


	available_workers_.push(worker_id);
}



void ParallelismConfig::send_all_available(int numtosend)
{
	while (available_workers_.size()>0)  {
		int sendtome = available_workers_.front();
		MPI_Send(&numtosend, 1, MPI_INT, sendtome, NUMPACKETS, MPI_COMM_WORLD);
		available_workers_.pop();
	}
}


void ParallelismConfig::call_for_help(int solver_type)
{

	MPI_Bcast(&solver_type, 1, MPI_INT, id(), MPI_COMM_WORLD);
	init_active_workers();

}

bool ParallelismConfig::have_available()
{
	if (available_workers_.size()==0) {
		return false;
	}
	else
	{
		return true;
	}

}

bool ParallelismConfig::have_active()
{
	bool yep = false;
	for (int ii=1; ii<this->numprocs_; ii++) {
		if (this->worker_status_[ii]==ACTIVE) {
			yep = true;
			break;
		}
	}
	return yep;
}

int ParallelismConfig::num_active()
{
	int num = 0;
	for (int ii=1; ii<this->numprocs_; ii++) {
		if (this->worker_status_[ii]==ACTIVE) {
			num++;
		}
	}
	return num;
}













void ProgramConfigBase::move_to_temp()
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

void ProgramConfigBase::move_to_called()
{

	chdir(this->called_dir().c_str());

	if (this->verbose_level()>=3)
		std::cout << "moved to called_dir '" << this->called_dir().string() << "'" << std::endl;
}


void ProgramConfigBase::PrintMetadata(boost::filesystem::path const& filename) const
{
	FILE *OUT = safe_fopen_write(filename);

	fprintf(OUT, "%s\n", VERSION);
	fprintf(OUT, "%s\n", called_dir_.c_str());
	fprintf(OUT, "%s\n", timer_.format().c_str());
	fprintf(OUT, "%d\n", this->num_procs());


	fclose(OUT);
}


void PrintPointTypeMapping(boost::filesystem::path const& filename)
{
	//The following lets us use words instead of numbers to indicate vertex type.
//enum {UNSET= 100, CRITICAL, SEMICRITICAL, MIDPOINT, ISOLATED, NEW, CURVE_SAMPLE_POINT, SURFACE_SAMPLE_POINT, REMOVED, PROBLEMATIC};

	FILE *OUT = safe_fopen_write(filename);

	fprintf(OUT,"11\n\n");
	fprintf(OUT,"Unset %d\n",Unset);
	fprintf(OUT,"Critical %d\n",Critical);
	fprintf(OUT,"Semicritical %d\n",Semicritical);
	fprintf(OUT,"Midpoint %d\n",Midpoint);
	fprintf(OUT,"Isolated %d\n",Isolated);
	fprintf(OUT,"New %d\n",New);
	fprintf(OUT,"Curve_sample_point %d\n",Curve_sample_point);
	fprintf(OUT,"Surface_sample_point %d\n",Surface_sample_point);
	fprintf(OUT,"Removed %d\n",Removed);
	fprintf(OUT,"Problematic %d\n",Problematic);
	fprintf(OUT,"Singular %d\n",Singular);

	fclose(OUT);

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

void BertiniRealConfig::splash_screen() const
{
	printf("\n BertiniReal(TM) v%s\n\n", VERSION);
	printf(" S. Amethyst with\n D.J. Bates, W. Hao, J.D. Hauenstein,\n A.J. Sommese, C.W. Wampler\n\n");
	printf("(using GMP v%d.%d.%d, MPFR v%s)\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());





	bertini_splash_screen();



}



void BertiniRealConfig::display_current_options() const
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
	std::cout << "output_directory base name: " << output_dir() << '\n';

	// Which symbolic Engine
	switch(symbolic_engine())
	  {
	  case SymEngine::Matlab:
	    std::cout << "Using Matlab as the symbolic engine.\n";
	    break;
	  case SymEngine::Python:
	    std::cout << "Using Python as the symbolic engine.\n";
	    break;
	  }

	 std::cout << "same_point_tol: " << same_point_tol() << std::endl;
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
			{"patch",required_argument, 0, 'A'},
			{"robustness",required_argument, 0, 'r'},{"r",required_argument, 0, 'r'},
			{"input",required_argument,	0, 'i'}, {"i",required_argument, 0, 'i'},
			{"version",		no_argument,			 0, 'v'}, {"v",		no_argument,			 0, 'v'},
			{"help",		no_argument,			 0, 'h'}, {"h",		no_argument,			 0, 'h'},
			{"mode",required_argument,0,'M'}, {"m",required_argument,0,'M'},
			{"symengine",required_argument,0,'E'},{"E",required_argument,0, 'E'},
			{"pycommand",required_argument,0,'P'},{"P",required_argument,0, 'P'},
			{"symnosubst",	no_argument,			 0, 't'},
			{"symallowsubst",	no_argument,			 0, 'T'},
			{"samepointtol",	required_argument,		 0, 'e'},
			{"ignoresing", no_argument, 0, 'w'},
			{"realify", no_argument, 0, 'R'},

			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;

		choice = getopt_long_only (argc, argv, "d:c:Dg:V:o:smp:S:i:rvhM:E:P:tTe:wA:R", // if followed by colon, requires option.  two colons is optional
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

			case 'A':
				user_patch(true);
				patch_filename_ = boost::filesystem::absolute(optarg);
				break;

			case 'S':
				user_sphere(true);
				this->bounding_sphere_filename_ = boost::filesystem::absolute(optarg);
				break;


			case 'i': // input filename
				input_filename_ = optarg;
				break;

			case 'r':
				robustness(atoi(optarg));
				break;

			case 'R':
				realify_ = true;
				break;

			case 'v':
				printf("\n BertiniReal (TM) v %s\n\n", VERSION);
				exit(0);
				break;

			case 'h':

				printf("\nBertiniReal (TM) v %s.\n\n", VERSION);
				printf("Online at bertinireal.com\n\n");
				printf("For support, \nfile an issue on Github at github.com/ofloveandhate/bertini_real/issues\n\n");
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

			case 'E': // symbolic Engine
			{
			  std::string use_engine = optarg;

			  // converting input_engine to lowercase
			  std::transform(use_engine.begin(),use_engine.end(),use_engine.begin(), ::tolower);
			  std::string check1="matlab";
			  std::string check2="python";

			  if (use_engine.compare(check1)==0)
				{
					this->engine_ = SymEngine::Matlab;
				}
			  else if (use_engine.compare(check2)==0)
				{
					this->engine_ = SymEngine::Python;
				}
				else
				{
					std::cout << "bad mode of symbolic engine.  acceptable options are Matlab and Python." << std::endl;
					exit(0);
				}

				break;
			}

			case 'P':
			{
				this->python_command_ = std::string(optarg) + " ";
			}

			case 't':
			{
				this->sym_prevent_subst(true);
				break;
			}

			case 'T':
			{
				this->sym_prevent_subst(false);
				break;
			}

			case 'e':
			{
				std::stringstream converter;
				converter << optarg;
				double d; converter >> d;
				this->same_point_tol(d);
				break;
			}

			case 'w':
			{
				this->ignore_singular(true);
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


void BertiniRealConfig::print_usage() const
{
	int opt_w = 20;
	int type_w = 11;
	int def_w = 9;
	int note_w = 30;
	auto line = [=](std::string name, std::string type, std::string def_val = " ", std::string note = " ")
		{std::cout << std::setw(opt_w) << std::left << name << std::setw(type_w) << type << std::setw(def_w) << def_val <<  std::setw(note_w) << note << "\n";};

	std::cout << "Bertini_real has the following options:\n\n";
	line("option","type","default","note");
	std::cout << '\n';

	line("-p -pi -projection", 	"string", 	" -- ", "projection file name");
	line("-i -input", 			"string", 	"input", "bertini input file name");
	line("-ns -nostifle", 		" -- ", 	" ", "turn off stifling of screen out when running system calls");

	line("-v -version", 		" -- ", 	" ", "get version number");
	line("-h -help", 			" --", 		" ", "print this help menu");
	line("-sphere -b", 			"string", 	" -- ", "name of sphere file");
	line("-r -robustness", 			"int", 	" 1 ", "use lower robustness to speed up computation -- but get worse results, probably");
	line("-debug", 				" -- ", 	" ", "make bertini_real wait 30 seconds for you to attach a debugger");
	line("-symengine -E", 		"string", 	"matlab", "select a symbolic engine.  choices are 'matlab' and 'python'");
	line("-pycommand -P", 		"string", 	"python", "indicate how python should be called.  default is 'python'");
	line("-symnosubst", 		" -- ", 	" ", "prevent substitution of subfunctions during deflation and other sym ops.");
	line("-symallowsubst", 		" -- ", 	" ", "allow substitution of subfunctions during deflation and other sym ops.  default");
	line("-samepointtol", 		"<double>", "1e-7" , "(scaled) infinity-norm distance between two points to be considered distinct");
	line("-nomerge", 		" -- ", "off" , "turn off merging for top-dimensional *curve* decompositions (does not affect surface decompositions).  no argument necessary, choosing this turns the nomerge option on.");
	line("-gammatrick -g", 		"bool", "0" , "use the complex gamma trick for all paths.  is this good?  does it even work at all?  does using this option produce complete garbage, or speed things up like racing stripes?  i don't know, but it's implemented and an option.  choose your own adventure.  enjoy.");
	line("-ignoresing", " -- ", "off", "ignore singular curve(s); only use if singular curves are naked.  no argument necessary, choosing this turns it on.");
	line("-realify", " -- ", "off", "change patch and discard imaginary parts where possible throughout decomposition.  no argument necessary, choosing this turns it on.");
	printf("\n\n\n");
	return;
}

void BertiniRealConfig::init()
{
	target_component_ = -2;
	target_dimension_ = -1;

	debugwait_ = false;
	max_deflations_ = 10;

	user_projection_ = false;
	projection_filename_ = "";

	orthogonal_projection_ = true;

	user_patch_ = false;
	patch_filename_ = "";

	user_sphere_ = false;
	bounding_sphere_filename_ = "";

	input_filename_ = "input";


	output_dir(boost::filesystem::absolute("./output"));


	stifle_membership_screen_ = true;
	stifle_text_ = " > /dev/null ";

	matlab_command_ = "matlab -nosplash -nodesktop -nojvm -r ";
	python_command_ = "python ";

	verbose_level(0); // default to 0


	use_gamma_trick_ = false;

	realify_ = false;
	merge_edges_ = true;

	primary_mode_ = BERTINIREAL;
	engine_ = SymEngine::Matlab; // setting default to Matlab symbolic engine
	prevent_sym_substitution_ = false;

	same_point_tol_ = 1e-7;

	ignore_singular_ = false;
	return;
}











void get_projection(vec_mp *pi,
					BertiniRealConfig const& program_options,
					int num_vars,
					int num_projections)
{


	for (int ii=0; ii<num_projections; ii++) {
		change_size_vec_mp(pi[ii], num_vars);  pi[ii]->size = num_vars;
	}



	//assumes the vector pi is already initialized
	if (program_options.user_projection()) {
		FILE *IN = safe_fopen_read(program_options.projection_filename()); // we are already assured this file exists, but safe fopen anyway.
		int tmp_num_vars;
		fscanf(IN,"%d",&tmp_num_vars); scanRestOfLine(IN);
		if (tmp_num_vars!=num_vars-1) {
			printf("the number of variables declared in the projection\nis not equal to the number of non-homogeneous variables in the problem\n");
			printf("please modify file to have %d coordinates, one per line.\n(imaginary part will be ignored if provided).\n",num_vars-1);
			abort();
		}

		for (int ii=0; ii<num_projections; ii++) {
			set_zero_mp(&pi[ii]->coord[0]);
			for (int jj=1; jj<num_vars; jj++) {
				set_zero_mp(&pi[ii]->coord[jj]);
				mpf_inp_str(pi[ii]->coord[jj].r, IN, 10);
				scanRestOfLine(IN);
			}
		}
		fclose(IN);
	}
	else{
		if (program_options.orthogonal_projection()) {
			mat_mp temp_getter;
			init_mat_mp2(temp_getter,0,0,1024);
			make_matrix_random_real_mp(temp_getter,num_projections, num_vars-1, 1024); // this matrix is ~orthogonal

			for (int ii=0; ii<num_projections; ii++) {
				set_zero_mp(&pi[ii]->coord[0]);
				for (int jj=1; jj<num_vars; jj++)
					set_mp(&pi[ii]->coord[jj], &temp_getter->entry[ii][jj-1]);

			}

			clear_mat_mp(temp_getter);

		}
		else
		{
			for (int ii=0; ii<num_projections; ii++) {
				set_zero_mp(&pi[ii]->coord[0]);
				for (int jj=1; jj<num_vars; jj++)
					get_comp_rand_real_mp(&pi[ii]->coord[jj]);//, &temp_getter->entry[ii][jj-1]);

			}

		}

	}


	return;
}



void get_patch(vec_mp *patch,
					BertiniRealConfig const& program_options,
					int num_vars)
{
	int num_patches = 1;

	for (int ii=0; ii<num_patches; ii++) {
		change_size_vec_mp(patch[ii], num_vars);  patch[ii]->size = num_vars;
	}



	//assumes the vector patch is already initialized
	if (program_options.user_patch()) {
		FILE *IN = safe_fopen_read(program_options.patch_filename()); // we are already assured this file exists, but safe fopen anyway.
		int tmp_num_vars;
		fscanf(IN,"%d",&tmp_num_vars); scanRestOfLine(IN);
		if (tmp_num_vars!=num_vars-1) {
			printf("the number of variables declared in the patch\nis not equal to the number of homogenized variables in the problem (#affine +1)\n");
			printf("please modify file to have %d real coordinates, one per line.\n(imaginary part will be ignored if provided).\n",num_vars);
			abort();
		}

		for (int ii=0; ii<num_patches; ii++) {
			for (int jj=0; jj<num_vars; jj++) {
				set_zero_mp(&patch[ii]->coord[jj]);
				mpf_inp_str(patch[ii]->coord[jj].r, IN, 10);
				scanRestOfLine(IN);
			}
		}
		fclose(IN);
	}
	else{
		for (int ii=0; ii<num_patches; ii++) {
			set_one_mp(&patch[ii]->coord[0]);
			for (int jj=1; jj<num_vars; jj++)
				set_zero_mp(&patch[ii]->coord[jj]);

		}

	}
}

