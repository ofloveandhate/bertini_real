bool create_Matlab_Crit_NoSubst(std::string const& filename,vec_mp *pi, boost::filesystem::path inputOutputName, int dim)
{

	FILE* OUT = safe_fopen_write("matlab_nullspace_system.m");

	fprintf(OUT,"filename = '%s';\n\n",filename.c_str());

	int num_vars = pi[1]->size-1;

	fprintf(OUT,"pi_vals = cell(%i,%i);\n\n",dim,num_vars);
	for (int ii = 0; ii < dim; ++ii)
	{
		for (int jj = 1; jj < pi[ii]->size; ++jj)
		{
			fprintf(OUT,"pi_vals{%i,%i} = '",ii,jj);
			mpf_out_str(OUT,10,0,pi[ii]->coord[jj].r);
			fprintf(OUT,"';\n");
		}
	}

	fprintf(OUT,"OutputName = '%s';\n\n",inputOutputName.c_str());

	fprintf(OUT, "crit_no_subst(filename,pi_vals,OutputName);\n\n");

	fprintf(OUT, "exit");

	fclose(OUT);

	return false;
}
