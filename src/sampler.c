#include "sampler.h"


int main(int argC, char *args[])
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES:                                                        *
 \***************************************************************/
{
    int num_vars=0;  //1=self-conjugate; 0=not
    char *inputName = NULL, *witnessSetName = NULL, *samplingName = NULL, *samplingNamenew = NULL;
    curveDecomp_d C;  //new data type; stores vertices, edges, etc.
    witness_set_d W;
    sample_d   S_old,S_new;
    //mat_d n_minusone_randomizer_matrix;

    ////
    //  begin the actual program
    ////
    
    if(setup_curveDecomp(argC, args, &inputName, &witnessSetName,&samplingName,&samplingNamenew,&C,&num_vars))
        return 1;  
    //< prints the welcome message, also gets the inputName, witnessSetName and C 
    srand(time(NULL));
    // essentials for using the bertini parser
    prog_t SLP;
	
    unsigned int currentSeed;
    int trackType, genType = 0, MPType,  sharpenOnly, needToDiff, remove_temp, useParallelDiff = 0;
    int my_id, num_processes, headnode = 0; // headnode is always 0
    int precision = 53;
    num_processes = 1;
    int num_var_gps = 0, userHom = 0;
    //end parser-bertini essentials
	
    parse_input(inputName, &trackType, &MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);
	
    preproc_data PPD;
    setupPreProcData("preproc_data", &PPD);
    
    num_var_gps = PPD.num_var_gp;
    num_vars = setupProg(&SLP, precision, MPType); // num_vars includes the number of homogeneous coordinates.
    // the number of homogeneous coordinates is the num_var_gps.
	
	
    printf("parsing witness set\n");
    witnessSetParse(&W, witnessSetName,num_vars); 
    printf("Loading sampling data\n");
    Load_sampling_data(&S_old,samplingName,C,num_vars);
    
    printf("The number of current sampling points is %d\n please enter the number of new sampling points",S_old.num_pts);
    while(1)
    {
        scanf("%d",&(S_new.num_pts));
        if(S_new.num_pts>S_old.num_pts)
        {
            printf("Sampler will generate %d sampling points for each edge\n",S_new.num_pts);
            break;
        }
        printf("The number you entered (%d) is less than the number of current sampling points (%d).",S_new.num_pts,S_old.num_pts);
        printf("Please re-enter...\n");
    }
    

    //Generate new sampling data
    generate_new_sampling_pts(&S_new,S_old, C,W, num_vars,currentSeed,MPType);

    //output
    output_sampling_data(S_new,samplingNamenew,num_vars);
    
    // clear memory
    free(inputName);
    free(witnessSetName);
	
    clear_witness_set(W);
    clear_sample_d(&S_old);
    clear_sample_d(&S_new);

  //TMP END
  return 0;
}

void generate_new_sampling_pts(sample_d *S_new,sample_d S_old, curveDecomp_d C,witness_set_d W, int num_vars,unsigned int currentSeed,int  MPType)
{
    witness_set_d Wnew;
	vec_mp         new_linears;
	vec_d L,startpt;
    int           i,j,k,cur_old_pt;
    comp_d     sl,se,temp;
    mat_mp n_minusone_randomizer_matrix;

	//should be reading in or creating an appropriately sized matrix right here.  so that we can use MPType2, this matrix should be of maximum precision.
    init_mat_mp(n_minusone_randomizer_matrix,1,1);
    make_matrix_ID_mp(n_minusone_randomizer_matrix, num_vars-1, num_vars-1);//this is likely of incorrect size

    init_witness_set_d(&Wnew);
    cp_patches(&Wnew,W);

    init_vec_mp(new_linears,num_vars);
    new_linears->size = num_vars;

    init_vec_d(L,num_vars); // what is the purpose of this?
    L->size = num_vars; // what is L?

    init_vec_d(startpt,num_vars);
    startpt->size = num_vars;

    S_new->num_E = S_old.num_E;
    S_new->V = (mat_d *)bmalloc(S_new->num_E * sizeof(mat_d));
    S_new->pV = (vec_d *)bmalloc(S_new->num_E * sizeof(vec_d));

    for(i=0;i<S_old.num_E;i++)
    {
        init_mat_d(S_new->V[i],num_vars,S_new->num_pts);//three points: left mid right
        init_vec_d(S_new->pV[i],S_new->num_pts); // what is pV?

        //copy left & right points
        set_d(&(S_new->pV[i]->coord[0]),&(S_old.pV[i]->coord[0]));
        set_d(&(S_new->pV[i]->coord[S_new->num_pts-1]),&(S_old.pV[i]->coord[S_old.num_pts-1]));

        for(j=0;j<num_vars;j++)
        {
            set_d(&(S_new->V[i]->entry[j][0]),&(S_old.V[i]->entry[j][0]));
            set_d(&(S_new->V[i]->entry[j][S_new->num_pts-1]),&(S_old.V[i]->entry[j][S_old.num_pts-1]));
        }
        //setup the length of interval of the edge projection
        sub_d(sl,&(S_old.pV[i]->coord[S_old.num_E-1]),&(S_old.pV[i]->coord[0])); // what is sl?

        for(k=1;k<S_new->num_pts-1;k++)
        {
					
            cur_old_pt=k*(S_old.num_pts-1)/(S_new->num_pts-1); // what does this line do?
            temp->r=(double)(k/(S_new->num_pts-1));temp->i=0.0; //what is happening here?
            mul_d(se,temp,sl); // what is se?
					
					//what is the purpose of this loop?  to set the start point for the solve?
            for(j=0;j<num_vars;j++)
            {
                set_d(&(L->coord[j]),&(C.E[i].pi->coord[j]));
                if(!j)
                    sub_d(&(L->coord[j]),&(L->coord[j]),&(S_old.pV[i]->coord[cur_old_pt]));
                d_to_mp(&(new_linears->coord[j]),&(C.E[i].pi->coord[j]));
                if(!j)
								{
                    sub_d(temp,&(L->coord[j]),se);
									d_to_mp(&(new_linears->coord[j]),temp);
								}
                set_d(&(startpt->coord[j]),&(S_old.V[i]->entry[j][cur_old_pt]));
            }
					
					// copy in the data to the source witness set?
            set_witness_set_d(&W, L,startpt,num_vars);
					
            lin_to_lin_solver_main(MPType,
															W,         // witness_set
															n_minusone_randomizer_matrix,
															&new_linears, //  the set of linears we will solve at.
															1, // the number of new linears.
															&Wnew); // the new data is put here!
					
					
					//Wnew has only the new point and new linear...  no info from W, save patches and num_variables, is copied over.
					
            //get the point coords
            set_zero_d(&(S_new->pV[i]->coord[k]));
            for(j=0;j<num_vars;j++)
            {
                set_d(&(S_new->V[i]->entry[j][k]),&(Wnew.W.pts[0]->coord[j]));
                mul_d(temp,&(S_new->V[i]->entry[j][k]),&(C.E[i].pi->coord[j]));
                add_d(&(S_new->pV[i]->coord[k]),&(S_new->pV[i]->coord[k]),temp);
            }    
        }
    }

    clear_witness_set(Wnew);
    clear_vec_d(L);
    clear_mat_d(n_minusone_randomizer_matrix);
    clear_vec_d(new_linears);

}

void set_witness_set_d(witness_set_d *W, vec_d L,vec_d pts,int num_vars)
{

    int ii;
	
	W->num_variables = num_vars;
	
	
	W->W_mp.num_pts=1;
	W->W_mp.pts=(point_mp *)bmalloc(sizeof(point_mp)); // apparently you can only pass in a single point to copy in.

	
    W->W.num_pts=1;
    W->W.pts=(point_d *)bmalloc(sizeof(point_d)); // apparently you can only pass in a single point to copy in.
    

    //initialize the memory
	init_point_d(W->W.pts[0],num_vars); W->W.pts[0]->size = num_vars;
	init_point_mp(W->W_mp.pts[0],num_vars); W->W_mp.pts[0]->size = num_vars;

	
	point_cp_d(W->W.pts[0],pts);
	point_d_to_mp(W->W_mp.pts[0],pts);
//    //read the witness points into memory
//    for (ii=0; ii < num_vars; ii++) 
//    {
//        W->W.pts[0]->coord[ii].r = pts->coord[ii].r;
//        W->W.pts[0]->coord[ii].i = pts->coord[ii].i;
//    }

    W->num_linears = 1;
    W->L = (vec_d *)bmalloc(sizeof(vec_d));
		 W->L_mp = (vec_mp *)bmalloc(sizeof(vec_mp));
	
   init_vec_d( W->L[0],   num_vars); W->L[0]->size = num_vars;
	 init_vec_mp(W->L_mp[0],num_vars); W->L_mp[0]->size = num_vars;

	vec_cp_d(W->L[0],L);
	vec_d_to_mp(W->L_mp[0],L);
//   for (ii=0; ii < num_vars; ii++)
//   {
//       W->L[0]->coord[ii].r = L->coord[ii].r;
//       W->L[0]->coord[ii].i = L->coord[ii].i;
//   }

}

void  output_sampling_data(sample_d S,char *samplingName,int num_vars)
{
    FILE *OUT =  fopen(samplingName, "w");
    int i,j,k,size,digits=15;
    char *fmt = NULL;
    size = 1 + snprintf(NULL, 0, "%%.%de %%.%de\n", digits, digits);
    // allocate size
    fmt = (char *)bmalloc(size * sizeof(char));
    // setup fmt & fmtb
    sprintf(fmt, "%%.%de %%.%de\n", digits, digits);
    // output the number of vertices
    fprintf(OUT,"%d\n\n",S.num_pts);
    
    for(i=0;i<S.num_E;i++)
    {
        for(j=0;j<S.num_pts;j++)
        {
            fprintf(OUT, fmt, S.pV[i]->coord[j].r, S.pV[i]->coord[j].i);
            for(k=0;k<num_vars;k++)
            {
                fprintf(OUT, fmt, S.V[i]->entry[k][j].r, S.V[i]->entry[k][j].i);
            }
            fprintf(OUT,"\n");
        }
    }
}

int  Load_sampling_data(sample_d *S,char *samlingName, curveDecomp_d C,int num_vars)
{
    int i,j,k,index;
    if(samlingName==NULL)
    {
        comp_d temp;
        S->num_E = C.num_E;
        S->num_pts = 3;
        S->V = (mat_d *)bmalloc(C.num_E * sizeof(mat_d));
        S->pV = (vec_d *)bmalloc(C.num_E * sizeof(vec_d));
        for(i=0;i<S->num_E;i++)
        {
            init_mat_d(S->V[i],num_vars,S->num_pts);//three points: left mid right
            //left point
	    index = C.E[i].left;
            init_vec_d(S->pV[i],S->num_pts);
            set_zero_d(&(S->pV[i]->coord[0]));
            for(j=0;j<num_vars;j++)
            {
                set_d(&(S->V[i]->entry[j][0]),&(C.V1[index].pt->coord[j]));
                mul_d(temp,&(S->V[i]->entry[j][0]),&(C.E[i].pi->coord[j]));
                add_d(&(S->pV[i]->coord[0]),&(S->pV[i]->coord[0]),temp);
            }
            //mid point
            for(j=0;j<num_vars;j++)
            {
                set_d(&(S->V[i]->entry[j][1]),&(C.E[i].midpt->coord[j]));
                mul_d(temp,&(S->V[i]->entry[j][1]),&(C.E[i].pi->coord[j]));
                add_d(&(S->pV[i]->coord[1]),&(S->pV[i]->coord[1]),temp);
            }
            //right point
            index = C.E[i].right;
            for(j=0;j<num_vars;j++)
            {
                set_d(&(S->V[i]->entry[j][2]),&(C.V1[index].pt->coord[j]));
                mul_d(temp,&(S->V[i]->entry[j][2]),&(C.E[i].pi->coord[j]));
                add_d(&(S->pV[i]->coord[2]),&(S->pV[i]->coord[2]),temp);
            }

        }
    }
    else
    {
        FILE *IN=safe_fopen_read(samlingName);
        fscanf(IN, "%d\n", &(S->num_pts));
 	
        S->num_E = C.num_E;
        S->V = (mat_d *)bmalloc(C.num_E * sizeof(mat_d));
        S->pV = (vec_d *)bmalloc(C.num_E * sizeof(vec_d));
       
        for(i=0;i<S->num_E;i++)
        {
            init_mat_d(S->V[i],num_vars,S->num_pts);
            init_vec_d(S->pV[i],S->num_pts);
            for(j=0;j<S->num_pts;j++)
            {
                fscanf(IN, "\n%lf %lf", &(S->pV[i]->coord[j].r), &(S->pV[i]->coord[j].i));
                for(k=0;k<num_vars;k++)
                {
                    fscanf(IN, "%lf %lf", &(S->V[i]->entry[k][j].r),  &(S->V[i]->entry[k][j].i));
                }
            }
        }
         
        fclose(IN);
    }
    return 0;
}


int  setup_curveDecomp(int argC, char *args[], char **inputName, char **witnessSetName, char **samplingName, char **samplingNamenew, curveDecomp_d *C,int *num_vars)
/***************************************************************\
 * USAGE:    setup curveDecomp structure and inputname 
 * ARGUMENTS:                                                    *
 * RETURN VALUES:           an integer                           *
 * NOTES:                                                        *
 \***************************************************************/
{
    FILE *IN;
    int strLength;
    char *directoryName=NULL,tmp_file[1000];
    DIR *dp;
    struct dirent *ep;

    printf("\n Sampler module for BertiniReal(TM) v%s\n\n", BERTINI_REAL_VERSION_STRING);
    printf(" D.J. Bates, D. Brake,\n W. Hao, J.D. Hauenstein,\n A.J. Sommese, C.W. Wampler\n\n");
    printf("(using GMP v%d.%d.%d, MPFR v%s)\n\n", 
           __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());

    // check for write privilege
    if (checkWritePrivilege())
    {
      printf("ERROR: BertiniReal does not have write privileges!\n");
      bexit(ERROR_WRITE_PRIVILEGE);
    }

    if (argC > 1 && args[1] != NULL && (!strcmp(args[1], "--help") || !strcmp(args[1], "-help"))) // help
    { // print information about Bertini
        printf("\nThis is a sampler module for BertiniReal v%s, developed by \n Dan J. Bates, Daniel Brake,\n Wenrui Hao, Jonathan D. Hauenstein, \n Andrew J. Sommmese, and Charles W. Wampler.\n\n", BERTINI_REAL_VERSION_STRING);
        printf("See ??? for details about BertiniReal.\n\n");
    }


    //setup the name of directory
    IN = safe_fopen_read("Dir_Name");
    fscanf(IN, "%d\n", &strLength);
    directoryName = (char *)bmalloc(strLength * sizeof(char));
    fgets(directoryName, strLength, IN);
    fclose(IN);
    
    //setup E structure from E.edge
    sprintf(tmp_file,  "%s/E.edge", directoryName);
    C->num_E = setup_edges(&(C->E),tmp_file,num_vars,inputName,directoryName);
    //setup V0 structure from V0.vert
    sprintf(tmp_file,  "%s/V0.vert", directoryName);
    C->num_V0 = setup_vertices(&(C->V0),tmp_file,*num_vars);

    //setup V1 structure from V1.vert
    sprintf(tmp_file,  "%s/V1.vert", directoryName);
    C->num_V1 = setup_vertices(&(C->V1),tmp_file,*num_vars);

    strLength = 1 + snprintf(NULL, 0, "%s/witness_set", directoryName);
    *witnessSetName = (char *)bmalloc(strLength * sizeof(char));
    sprintf(*witnessSetName, "%s/witness_set", directoryName);

    dp = opendir (directoryName);
    if (dp != NULL)
    {
        while ( (ep = readdir (dp)) )
           if (ep->d_name[0] != '.' && !strcmp(ep->d_name,"samp.dat"))
	   {
               strLength = 1 + snprintf(NULL, 0, "%s/samp.dat", directoryName);
               *samplingName = (char *)bmalloc(strLength * sizeof(char));
               sprintf(*samplingName, "%s/samp.dat", directoryName);
           }                 

    }
    (void) closedir (dp);

    strLength = 1 + snprintf(NULL, 0, "%s/samp.dat", directoryName);
    *samplingNamenew = (char *)bmalloc(strLength * sizeof(char));
    sprintf(*samplingNamenew, "%s/samp.dat", directoryName);
    if(! C->num_E)
    {
         printf("sampler will not generate sampling data since there is no edges\n");
         return 1;
    }

    return 0;
}

int setup_edges(edge_d **E,char *INfile,int *num_vars, char **inputName, char *directoryName)
//setup the vertex structure
{
    FILE *IN = safe_fopen_read(INfile);
    int num,i,j,strLength;
    char *input_deflated_Name=NULL;

    fscanf(IN, "%d\n", num_vars);
    fscanf(IN, "%d\n", &num);
    fscanf(IN, "%d\n", &strLength);
    
    input_deflated_Name = (char *)bmalloc((strLength+1) * sizeof(char));
    fgets(input_deflated_Name, (strLength+1), IN);
    
    strLength = 1 + snprintf(NULL, 0, "%s/%s", directoryName,input_deflated_Name);
    *inputName = (char *)bmalloc(strLength * sizeof(char));
    sprintf(*inputName,  "%s/%s", directoryName,input_deflated_Name);
    
    *E=(edge_d*) bmalloc(num*sizeof(edge_d));
    for(i=0;i<num;i++)
    {
        fscanf(IN,"%d\n",&((*E)[i].left));
        fscanf(IN,"%d\n",&((*E)[i].right));
        init_point_d((*E)[i].midpt,*num_vars);
        for(j=0;j<*num_vars;j++)
            fscanf(IN, "%lf %lf", &((*E)[i].midpt->coord[j].r), &((*E)[i].midpt->coord[j].i));

        init_point_d((*E)[i].pi,*num_vars);
        for(j=0;j<*num_vars;j++)
            fscanf(IN, "%lf %lf", &((*E)[i].pi->coord[j].r), &((*E)[i].pi->coord[j].i));

    }
    fclose(IN);
    return num;
}


int setup_vertices(vertex_d **V,char *INfile,int num_vars)
//setup the vertex structure
{
    FILE *IN = safe_fopen_read(INfile);
    int num,i,j;
    fscanf(IN, "%d\n\n", &num);
    *V=(vertex_d* )bmalloc(num*sizeof(vertex_d));
    for(i=0;i<num;i++)
    {
        init_point_d((*V)[i].pt,num_vars);
        for(j=0;j<num_vars;j++)
        {
            fscanf(IN, "%lf %lf", &((*V)[i].pt->coord[j].r), &((*V)[i].pt->coord[j].i));
        }
    }
    fclose(IN);
    return num;
}
