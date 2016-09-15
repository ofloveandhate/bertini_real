#pragma once


#include "cells/vertex.hpp"
#include "nag/witness_set.hpp"

/**
 \brief the main structure for storing vertices.  
 
 The VertexSet is bertini_real's main method for storing data.  We essentially construct a graph of vertices, consisting of edges and faces.
 
 there are methods in place to add vertices, and perform lookups.
 */
class VertexSet
{

protected:
	
	vec_mp *projections_; ///< a pointer array of projection vectors.
	int num_projections_; ///< the number of projections.  this should match the dimension of the object being decomposed.
	int curr_projection_; ///< the projection currently being used.
	
	int curr_input_index_; ///< the index of the current input file.
	std::vector< boost::filesystem::path > filenames_; ///< the set of filenames from which vertices arise.
	
	std::vector<Vertex> vertices_;  ///< the main storage of points in the real numerical cellular decomposition.
	size_t num_vertices_; ///< the number of vertices found so far.
	
	int num_natural_variables_;  ///< the number of natural variables appearing in the problem to solve.
	
	
	double same_point_tolerance_;
	mpf_t abs_;
	mpf_t zerothresh_;
	comp_mp diff_;
	vec_mp checker_1_;
	vec_mp checker_2_;
	
	tracker_config_t * T_;
public:
	
	/**
	 get the tracker config struct
	 
	 \return a pointer to the tracker config.
	 */
	tracker_config_t * T() const
	{
		return T_;
	}
	
	/**
	 set the tracker config pointer
	 
	 \param new_T a pointer to the tracker config.
	 */
	void set_tracker_config(tracker_config_t * new_T)
	{
		T_ = new_T;
	}
	
	/**
	 get the tolerance for two points being the same
	 \return the L2 tolerance for whether two points are the same.
	 */
	double same_point_tolerance() const
	{
		return same_point_tolerance_;
	}
	
	/**
	 \brief set the tolerance for points being the same
	 \param new_tolerance The new tolerance.  This is for telling whether two points are the same.
	 */
	void set_same_point_tolerance(double new_tolerance)
	{
		same_point_tolerance_ = new_tolerance;
	}
	
	/**
	 \brief get the currently active projection
	 
	 \return the index of the current projection
	 */
	inline int curr_projection() const
	{
		return curr_projection_;
	}
	
	
	
	
	/**
	 \brief get the number of natural variables
	 
	 \return the number of natural variables
	 */
	inline int num_natural_variables() const
	{
		return num_natural_variables_;
	}
	
	
	/**
	 \brief get the number of vertices.
	 
	 \return the number of vertices stored in this Vertex set
	 */
	inline unsigned int num_vertices() const
	{
		return num_vertices_;
	}
	
	/**
	 \return the ith Vertex, or a reference to it.
	 \param index The index of the vertex to get.
	 */
	const Vertex& GetVertex(unsigned int index) const
	{
		if (index >= vertices_.size()) {
			throw std::out_of_range("trying to access Vertex out of range in VertexSet.");
		}
		return vertices_[index];
	}

	/**
	 \return the ith Vertex, or a reference to it.
	 \param index The index of the vertex to get.
	 */
	const Vertex& operator[](unsigned int index) const
	{
		if (index >= vertices_.size()) {
			throw std::out_of_range("trying to access Vertex out of range in VertexSet.");
		}
		return vertices_[index];
	}
	
	/**
	 \return the ith Vertex, or a reference to it.
	 \param index The index of the vertex to get.
	 */
	Vertex & operator[](unsigned int index)
	{
		if (index >= vertices_.size()) {
			throw std::out_of_range("trying to access Vertex out of range in VertexSet.");
		}
		return vertices_[index];
	}
	
	
	boost::filesystem::path filename(unsigned int index) const
	{
		if (index >= filenames_.size()) {
			throw std::out_of_range("trying to access filename out of range");
		}
		
		
		return filenames_[index];
	}
	
	
	void print_to_screen() const; ///< operator for displaying information to screen
	
	
	/**
	 \brief add a new Vertex to the set.
	 
	 \param new_vertex Vertex to add to the set.
	 \return the index of the added Vertex
	 */
	int add_vertex(const Vertex & new_vertex);
	
	
	/**
	 \brief create a VertexSet from a file.
	 
	 Read in a VertexSet from a file.
	 
	 \param INfile the file to parse and store in a VertexSet
	 \return the number of vertices read in.
	 */
	int setup_vertices(boost::filesystem::path INfile);
	
	
	
	
	
	VertexSet(){
		init();
	}
	
	VertexSet(int num_vars){
		init();
		
		set_num_vars(num_vars);
	}
	

	

	
	
	VertexSet & operator=( const VertexSet& other) {
		init();
		copy(other);
		return *this;
	}
	
	VertexSet(const VertexSet &other)
	{
		init();
		copy(other);
	}
	
	~VertexSet()
	{
		VertexSet::clear();
	}
	
	
	
	/**
	 sets the number of variables for the VertexSet.
	 
	 \param num_vars the number of {\em natural} variables
	 */
	void set_num_vars(int num_vars)
	{
		this->num_natural_variables_ = num_vars;
		
		change_size_vec_mp(checker_1_, num_vars-1);
		change_size_vec_mp(checker_2_, num_vars-1);
		checker_1_->size = checker_2_->size = num_vars-1;
	}
	
	
	
	/**
	 \brief write VertexSet to a file, readable by bertini_real again.
	 
	 
	 write the VertexSet to a file.
	 \see setup_vertices
	 
	 Output Vertex structure as follows:
	 
	 [
	 num_vertices num_projections num_natural_variables filenames.size()
	 
	 the projections, as bertini points
	 
	 the names of the files as
	 length_of_name name   pairs
	 
	 then the points as
	 
	 
	 num_variables in point \\
	 coordinates
	 
	 num_projection_coordinates \\
	 projection values 
	 
	 filename_index 
	 
	 type
	]
	 
	 
	 \param outputfile the name of the file to write the VertexSet to.
	 */
	void print(boost::filesystem::path outputfile) const;
	
	
	
	/**
	 set the name of the current input file.  while it is set to this, all added vertices will inherit the index of this name.
	 
	 \return the index of the set filename.
	 \param el_nom the name of the file
	 */
	int set_curr_input(boost::filesystem::path el_nom){
		
		int nom_index = -1;
        
		if ( el_nom.string().compare("unset_filename")==0 )
		{
			throw std::logic_error("trying to set curr_input from unset_filename");
		}
		
		
		int counter = 0;
		for (std::vector<boost::filesystem::path>::iterator ii = filenames_.begin(); ii!= filenames_.end(); ++ii) {
			
			if (*ii == el_nom) {
				nom_index = counter;
				break;
			}
			counter++;
		}
		
		if (nom_index==-1) {
			filenames_.push_back(el_nom);
			nom_index = counter;
		}
		
		curr_input_index_ = nom_index;
		return nom_index;
	}
	
	
    

    /**
	 \brief find the index of a point.
	 
	 find the index of a point.
	 first, we check the active points, then the inactive, then give up.  the method dehomogenizes the points, and checks only the natural variables
	 
	 \return the index of the testpoint, or -1 if it is not found.
	 \param testpoint the mp point to find.
	 */
    int search_for_point(vec_mp testpoint);
	
	
	/**
	 \brief find the index of a point among active (non-removed) points only.
	 
	 find the index of a point.
	 
	 
	 \return the index of the testpoint, or -1 if it is not found.
	 \param testpoint the mp point to find.
	 */
    int search_for_active_point(vec_mp testpoint);
	
	
	
	/**
	 \brief find the index of a point among removed points only.
	 
	 find the index of a point.
	 
	 
	 \return the index of the testpoint, or -1 if it is not found.
	 \param testpoint the mp point to find.
	 */
    int search_for_removed_point(vec_mp testpoint);
    
	
	/**
	 \brief Compute projections values, and midpoints, of a set of points with respect to a projection.
	 
	 This function computes \f$\pi(x)\f$ for each of the points \f$x\f$ in WitnessSet W.  Then it sorts them, and computes averages.
	 
	 The output is stored in crit_downstairs and midpoints_downstairs, both pre-initialized vec_mp's.  This function also produces a std::vector<int> named index_tracker which contains the sorting of W according to \f$\pi(W)\f$.
	 
     \param W witness				set containing points of which we wish to retrieve projections values.
     \param crit_downstairs			the projection values of the input witness set, sorted for uniqueness and increasingness.
     \param midpoints_downstairs	the bisection of each interval in crit_downstairs.
     \param index_tracker			the indices of the points in W.
     \param pi						the projection we are retrieving projection values with respect to.
	 \param T pointer to a Bertini tracker_config_t object, holding the necessary settings for this method.
     \return the integer SUCCESSFUL.
     */
    int compute_downstairs_crit_midpts(const WitnessSet & W,
                                       vec_mp crit_downstairs,
                                       vec_mp midpoints_downstairs,
                                       std::vector< int > & index_tracker,
									   vec_mp pi,
									   tracker_config_t * T);
	
	
	/**
	 
	 sets the value of the [current] projection for each Vertex which has index in the set of relevant indices.
	 
	 \param relevant_indices set of Vertex indices
	 \param new_value the new value you want to set the projection value to
	 \return a vector of indices for which this operation failed, because the new and old values were too far away from each other.
	 */
    std::vector<int>  assert_projection_value(const std::set< int > & relevant_indices, comp_mp new_value);
	
	/**
	 
	 sets the value of the [proj_index] projection for each Vertex which has index in the set of relevant indices.
	 
	 \param relevant_indices set of Vertex indices
	 \param new_value the new value you want to set the projection value to
	 \param proj_index the index of the projection you want to assert.
	 \return a vector of indices for which this operation failed, because the new and old values were too far away from each other.
	 */
    std::vector<int>  assert_projection_value(const std::set< int > & relevant_indices, comp_mp new_value, int proj_index);
	
	
	/**
	 
	 
	 set the current projection.
	 
	 \return the index of the projection
	 \param new_proj the projection to set as current
	 */
	int set_curr_projection(vec_mp new_proj){
        
        int proj_index = get_proj_index(new_proj);
        
        if (proj_index==-1) {
            int init_size = new_proj->size;
            new_proj->size = this->num_natural_variables_;
            
			proj_index = add_projection(new_proj);
            
            new_proj->size = init_size;
		}
		
        curr_projection_ = proj_index;
        
//        std::cout << "curr_projection is now: " << curr_projection << std::endl;
		return proj_index;
	}
	
    
	/**
	 \brief query the index of a projection.
	 
	 Want to find out the index of a projection in this VertexSet? pass it into this function.
	 
	 \param proj the projection to query
	 \return the index of the projection, or -1 if it doesn't exist.
	 */
    int get_proj_index(vec_mp proj) const
    {
        int init_size = proj->size;
        proj->size = this->num_natural_variables_;
        
        
        int proj_index = -1;
        
        for (int ii=0; ii<num_projections_; ii++) {
			if (isSamePoint_inhomogeneous_input(projections_[ii],proj,same_point_tolerance_)) {
				proj_index = ii;
				break;
			}
		}
        
        proj->size = init_size;
        
        return proj_index;
    }
    
    
	/**
	 \brief add a projection to the set, and get its index.
	 
	 Add a projection to the set, and get its index.  Does not test for uniqueness of the projection, assumes it is not in there yet.
	 
	 \return the index of the added projection
	 \param proj the projection to add.
	 
	 */
	int add_projection(vec_mp proj){
		
		if (this->num_projections_==0) {
			projections_ = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else{
			this->projections_ = (vec_mp *)br_realloc(this->projections_, (this->num_projections_+1) * sizeof(vec_mp));
		}
		
		
		init_vec_mp2(this->projections_[num_projections_],num_natural_variables_,T_->AMP_max_prec);
		this->projections_[num_projections_]->size = num_natural_variables_;
		
		
		if (proj->size != num_natural_variables_) {
			vec_mp tempvec;
			init_vec_mp2(tempvec,num_natural_variables_,T_->AMP_max_prec); tempvec->size = num_natural_variables_;
			for (int kk=0; kk<num_natural_variables_; kk++) {
				set_mp(&tempvec->coord[kk], &proj->coord[kk]);
			}
			vec_cp_mp(projections_[num_projections_], tempvec);
			clear_vec_mp(tempvec);
		}
		else
		{
			vec_cp_mp(projections_[num_projections_], proj);
		}
		

		num_projections_++;
		
		return num_projections_;
	}
	
	
	
	/**
	 \brief single target mpi send
	 
	 Send a VertexSet to a single target process.
	 
	 \param target who to send to
	 \param mpi_config the current mpi configuration
	 */
	void send(int target, ParallelismConfig & mpi_config) const;
	
	
	
	/**
	 \brief single source mpi receive
	 
	 Receive a VertexSet from a single source.
	 
	 \param source the source of the receive
	 \param mpi_config the current mpi configuration
	 */
	void receive(int source, ParallelismConfig & mpi_config);
	
	
	/**
	 reset the set to empty.
	 */
	void reset()
	{
		for (int ii=0; ii<num_projections_; ii++) {
			clear_vec_mp(projections_[ii]);
		}
		num_projections_ = 0;
		
		filenames_.resize(0);
		
		vertices_.resize(0);
		num_vertices_ = 0;
		
		clear();
		init();
	}
protected:
	
	void init()
	{
		T_ = NULL;
		
		same_point_tolerance_ = 1e-5;
		num_projections_ = 0;
		projections_ = NULL;
		curr_projection_ = -1;
		
		curr_input_index_ = -2;
		
		this->num_vertices_ = 0;
		this->num_natural_variables_ = 0;
		
		init_vec_mp(checker_1_,0);
		init_vec_mp(checker_2_,0);
		
		
		
		init_mp(this->diff_);

		mpf_init(abs_);
		mpf_init(zerothresh_);
		mpf_set_d(zerothresh_, 1e-8);
	}
	
	
	void copy(const VertexSet &other)
	{
		set_tracker_config(other.T());
		
		
		this->curr_projection_ = other.curr_projection_;
		if (this->num_projections_==0 && other.num_projections_>0) {
			projections_ = (vec_mp *) br_malloc(other.num_projections_*sizeof(vec_mp));
		}
		else if(this->num_projections_>0 && other.num_projections_>0) {
			projections_ = (vec_mp *) br_realloc(projections_,other.num_projections_*sizeof(vec_mp));
		}
		else if (this->num_projections_>0 && other.num_projections_==0){
			for (int ii=0; ii<this->num_projections_; ii++) {
				clear_vec_mp(projections_[ii]);
			}
			free(projections_);
		}
		
		for (int ii=0; ii<other.num_projections_; ii++) {
			if (ii>=this->num_projections_){
				init_vec_mp2(projections_[ii],1,1024); projections_[ii]->size = 1;
			}
			vec_cp_mp(projections_[ii],other.projections_[ii]);
		}
		
		this->num_projections_ = other.num_projections_;
		
		
		curr_input_index_ = other.curr_input_index_;
		filenames_ = other.filenames_;
		
		
		this->num_vertices_ = other.num_vertices_;
		this->num_natural_variables_ = other.num_natural_variables_;
		
		this->vertices_ = other.vertices_;
		
		vec_cp_mp(this->checker_1_,other.checker_1_);
		vec_cp_mp(this->checker_2_,other.checker_2_);
		
	}
	
	void clear()
	{
		clear_vec_mp(checker_1_);
		clear_vec_mp(checker_2_);
		
		clear_mp(diff_);
		
		mpf_clear(abs_);
		mpf_clear(zerothresh_);
		
		for (int ii=0; ii<num_projections_; ii++) {
			clear_vec_mp(projections_[ii]);
		}
		free(projections_);
	}

};




