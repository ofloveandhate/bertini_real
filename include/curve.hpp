#ifndef CURVE_H
#define CURVE_H

/** \file curve.hpp */
#include <cstddef>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include <mpfr.h>
#include <mpf2mpfr.h>


#include <set>





#include "bertini_headers.hpp"

#include "checkSelfConjugate.hpp"

#include "fileops.hpp"

#include "data_type.hpp"
#include "isosingular.hpp"
#include "programConfiguration.hpp"



#include "nullspace.hpp"
#include "solver_sphere.hpp"




/**
 \brief 1-cell.
 
 the edge data type.  has three indices: left, right, midpt.
 */
class Edge : public Cell
{
	int left_;  ///< index into vertices
	int right_; ///< index into vertices
	
	std::vector< int > removed_points_;
	
	
public:
	
	
	typedef std::vector< int >::iterator removed_iterator;
	typedef std::vector< int >::const_iterator removed_const_iterator;
	
	/**
	 \return get an iterator to the beginning of the set of removed points, which were merged away in the merge operation.
	 */
	inline removed_iterator removed_begin(){return removed_points_.begin();}
	
	/**
	 \return get an iterator to the beginning of the set of removed points, which were merged away in the merge operation.
	 */
	inline removed_const_iterator removed_begin() const {return removed_points_.begin();}
	
	/**
	 \return get an iterator to the end of the set of removed points, which were merged away in the merge operation.
	 */
	inline removed_iterator removed_end(){return removed_points_.end();}
	
	/**
	 \return get an iterator to the end of the set of removed points, which were merged away in the merge operation.
	 */
	inline removed_const_iterator removed_end() const {return removed_points_.end();}
	
	
	
	
	/**
	 \brief adds a point as a removed point.  tacks on to the end of the vector
	 
	 \param new_removed_point the index of the point to add
	 \return the index of the point
	 */
	int add_removed_point(int new_removed_point)
	{
		removed_points_.push_back(new_removed_point);
		return new_removed_point;
	}
	
	
	/**
	 \brief get the right point
	 
	 \return the index of the right point
	 */
	inline int right() const
	{
		return right_;
	}
	
	
	/**
	 \brief set the left point
	 \param new_right the new index of the left point
	 \return the index of the left point
	 */
	int right(int new_right)
	{
		return right_ = new_right;
	}
	
	
	
	
	
	/**
	 \brief get the left point
	 
	 \return the index of the left point
	 */
	inline int left() const
	{
		return left_;
	}
	
	
	/**
	 \brief set the left point
	 \param new_left the new index of the left point
	 \return the index of the left point
	 */
	int left(int new_left)
	{
		return left_ = new_left;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	 default constructor
	 */
	Edge() : Cell()
	{
		left_ = right_ = -1;
	}
	
	
	/**
	 \brief construct edge from left mid and right indices
	 \param new_left the new left index for constructed edge
	 \param new_midpt the new mid index for constructed edge
	 \param new_right the new right index for constructed edge
	 */
	Edge(int new_left, int new_midpt, int new_right)
	{
		left(new_left);
		right(new_right);
		midpt(new_midpt);
	}
	
	
	
	
	/**
	 check whether the edge is degenerate
	 \return true if the edge is degenerate, false if not.
	 */
	inline bool is_degenerate()
	{
		if ((left() == right()) && (left()==midpt()) && (right()==midpt()))
			return true;
		else
			return false;
	}
	
	
	
	/**
	 \brief send to a single target
	 \param target to whom to send this edge
	 \param mpi_config the current state of parallelism
	 */
	void send(int target, ParallelismConfig & mpi_config)
	{
		int * buffer = new int[4];
		
		buffer[0] = left();
		buffer[1] = midpt();
		buffer[2] = right();
		buffer[3] = removed_points_.size();
		
		MPI_Send(buffer, 4, MPI_INT, target, EDGE, mpi_config.comm());
		
		delete [] buffer;
		
		buffer = new int[removed_points_.size()];
		for (unsigned int ii=0; ii!=removed_points_.size(); ii++) {
			buffer[ii] = removed_points_[ii];
		}
		MPI_Send(buffer, removed_points_.size(), MPI_INT, target, EDGE, mpi_config.comm());
		delete [] buffer;
		
		
	}
	
	
	/**
	 \brief receive from a single source
	 \param source from whom to receive this edge
	 \param mpi_config the current state of parallelism
	 */
	void receive(int source, ParallelismConfig & mpi_config)
	{
		MPI_Status statty_mc_gatty;
		int * buffer = new int[4];
		MPI_Recv(buffer, 4, MPI_INT, source, EDGE, mpi_config.comm(), &statty_mc_gatty);
		
		left(buffer[0]);
		midpt(buffer[1]);
		right(buffer[2]);
		int temp_num_removed = buffer[3];
		
		
		delete [] buffer;
		
		buffer = new int[temp_num_removed];
		MPI_Recv(buffer, temp_num_removed, MPI_INT, source, EDGE, mpi_config.comm(), &statty_mc_gatty);
		for (int ii=0; ii<temp_num_removed; ii++) {
			removed_points_.push_back(buffer[ii]);
		}
		
		delete [] buffer;
		
	}
	
	
	/**
	 \brief get edge from input stream.  this function is defunct, and needs implementation apparently.
	 \param is the stream from whom to read
	 */
	virtual void read_from_stream( std::istream &is )
	{
		
		is >> left_ >> midpt_ >> right_;
	}
	
	friend std::ostream & operator<<(std::ostream & out, const Edge & E){
		out << E.left_ << " " << E.midpt_ << " " << E.right_;
		return out;
	}
	
};









/**
 \brief A bertini_real cell curve Decomposition.
 
 includes methods to add vertices, look up vertices, etc
 */
class Curve : public Decomposition
{
	
	std::vector<unsigned int> num_samples_each_edge_; ///< the number of samples on each edge, where there is a strict correspondence between elements of this vector and edges.
	
	std::vector< std::vector<int > > sample_indices_; ///< the indices of the vertices for the samples on the edges.
	
	
	std::vector<Edge> edges_; ///< The edges (1-cells) computed by Bertini_real
	size_t      num_edges_;  ///< How many edges this Decomposition currently has.  This could also be inferred from edges.size()
	
	

	
public:
	
	
	/**
	 query how many samples there are on an edge
	 
	 \throws out of range if there aren't as many edges as needed to do the lookup.
	 
	 \param edge_index the index of the edge, for which to query.
	 \return the number of samples on edge i
	 */
	unsigned int num_samples_on_edge(unsigned int edge_index)
	{
		if (edge_index >= num_samples_each_edge_.size()) {
			throw std::out_of_range("trying to access num_samples_each_edge_ of out of range index");
		}
		
		return num_samples_each_edge_[edge_index];
	}
	
	
	
	
	/**
	 get the index of sample on edge, at index position
	 
	 \throws out of range if there aren't as many edges as needed to do the lookup, or if there aren't that many samples..
	 
	 \param edge_index the index of the edge, for which to query.
	 \param vertex_index the position of the vertex you want the index of.
	 \return the index of the point, in the vertex set
	 */
	int sample_index(unsigned int edge_index, unsigned int vertex_index)
	{
		if (edge_index>=num_edges()) {
			std::stringstream ss;
			ss << "when accessing sample index, trying to access EDGE index out of range" << "\n" << "edge_index: "  << edge_index << " vertex_index: " << vertex_index;
			throw std::out_of_range(ss.str());

		}
		
		if (vertex_index >= sample_indices_[edge_index].size()) {
			std::stringstream ss;
			ss << "when accessing sample index, trying to access VERTEX index out of range" << "\n" << "edge_index: "  << edge_index << " vertex_index: " << vertex_index;
			throw std::out_of_range(ss.str());
		}
		
		return sample_indices_[edge_index][vertex_index];
		
	}
	
	/** 
	 \brief get the number of edges in the curve
	 
	 \return the number of edges in the Decomposition
	 */
	unsigned int num_edges() const
	{
		return num_edges_;
	}
	
	
	/**
	 \brief get an edge of the curve.  they are stored in the order in which they were computed.
	 
	 \param index the index of the edge you want.
	 
	 \return the ith edge
	 */
	Edge get_edge(unsigned int index) const
	{
		
		if (index>=num_edges_) {
			throw std::out_of_range("trying to get Edge out of range");
		}
		
		return edges_[index];
	}
	
	
	
	/**
	 \brief Get the set of all Vertex indices to which this Decomposition refers.
	 
	 This set is computed by reading all the edges.
	 
	 \return a std::set of ALL the Vertex indices occuring in this Decomposition.
	 */
	std::set< int > all_edge_indices()
    {
        std::set< int > index_set;
        
        for (auto iter=edges_.begin(); iter!=edges_.end(); iter++) {
            index_set.insert(iter->left());
            index_set.insert(iter->midpt());
            index_set.insert(iter->right());
        }
        
        return index_set;
    }
	
	
	
	
	/**
	 \brief Commit an edge to the curve.
	 
	 Increments the edge counter, and pushes the edge to the pack of the vector of edges.
	 
	 \return the index of the edge just added.
	 \param new_edge the edge to add.
	 */
	int add_edge(Edge new_edge)
	{
		num_edges_++;
		edges_.push_back(new_edge);
		return num_edges_-1; // -1 to correct for the fencepost problem
	}
	
	/**
	 \brief Read a Decomposition from text file.
	 
	 This function takes in the name of a folder containing a Decomposition to set up, and parses the two folders "containing_folder/decomp" and "containing_folder/E.edge".  Note that the VertexSet is set up separately from decompositions.
	 
	 \return the number 1.  Seems stupid.
	 
	 \param containing_folder The name of the folder containing the Decomposition.
	 */
	int setup(boost::filesystem::path containing_folder);
	
	
	
	/**
	 \brief open a folder as an edges file, and parse it into the vector of edges in this Decomposition.
	 
	 \return the number of edges added
	 \param INfile the name of the file to open and read as an edge file.
	 */
	int setup_edges(boost::filesystem::path INfile);
	
	
	/**
	 
	 \brief Write all the edges in this Decomposition to a file.
	 
	 the format is 
	 
	 [
	 num_edges
	 
	 left mid right
	 ]
	 
	 \param outputfile The name of the file to write to.
	 */
	void print_edges(boost::filesystem::path outputfile);
	
	
	/**
	 \brief print the complete curve Decomposition, including the base Decomposition class.
	 
	 \param base The name of the base folder to print the Decomposition to.
	 */
	void print(boost::filesystem::path base);
	
	
	
	
	
	
	
	/**
	 \brief Search this Decomposition for an nondegenerate edge with input point index as midpoint.
	 
	 \param ind The index to search for.
	 \return The index of the found edge, or -10 if no edge had the point as midpoint.
	 */
	int nondegenerate_edge_w_midpt(int ind)
	{
		
		for (unsigned int ii=0; ii<num_edges_; ii++){
			if (this->edges_[ii].midpt() == ind){
				
				if (edges_[ii].is_degenerate()) {
					continue;
				}
				
				return ii;
			}
		}
		
		return -10;
	}
	
	
	
	/**
	 \brief Search this Decomposition for an nondegenerate edge with input point index as left point.
	 
	 \param ind The index to search for.
	 \return The index of the found edge, or -11 if no edge had the point as left point.
	 */
	int nondegenerate_edge_w_left(int ind)
	{
		
		for (unsigned int ii=0; ii<num_edges_; ii++){
			if (this->edges_[ii].left() == ind){
				
				if (edges_[ii].is_degenerate()) {
					continue;
				}
				
				return ii;
			}
		}
		
		return -11;
	}
	
	
	/**
	 \brief Search this Decomposition for an nondegenerate edge with input point index as right point.
	 
	 \param ind The index to search for.
	 \return The index of the found edge, or -12 if no edge had the point as right point.
	 */
	int nondegenerate_edge_w_right(int ind)
	{
		
		for (unsigned int ii=0; ii<num_edges_; ii++){
			if (this->edges_[ii].right() == ind){
				
				if (edges_[ii].is_degenerate()) {
					continue;
				}
				
				return ii;
			}
		}
		
		return -12;
	}
	
	
	
	/**
	 \brief Search this Decomposition for any edge (degenerate or not) with input point index as mid point.
	 
	 \param ind The index to search for.
	 \return The index of the found edge, or -10 if no edge had the point as mid point.
	 */
	int edge_w_midpt(int ind)
	{
		
		for (unsigned int ii=0; ii<num_edges_; ii++){
			if (this->edges_[ii].midpt() == ind){
				return ii;
			}
		}
		
		return -10;
	}
	
	
	/**
	 \brief Search this Decomposition for any edge (degenerate or not) with input point index as left point.
	 
	 \param ind The index to search for.
	 \return The index of the found edge, or -11 if no edge had the point as left point.
	 */
	int edge_w_left(int ind)
	{
		
		for (unsigned int ii=0; ii<num_edges_; ii++){
			if (this->edges_[ii].left() == ind){
				return ii;
			}
		}
		
		return -11;
	}
	
	
	/**
	 \brief Search this Decomposition for any edge (degenerate or not) with input point index as right point.
	 
	 \param ind The index to search for.
	 \return The index of the found edge, or -12 if no edge had the point as right point.
	 */
	int edge_w_right(int ind)
	{
		
		for (unsigned int ii=0; ii<num_edges_; ii++){
			if (this->edges_[ii].right() == ind){
				return ii;
			}
		}
		
		return -12;
	}
	
	
	
	/**
	 \brief Search this Decomposition for any edge (degenerate or not) with input point index as a removed point.
	 
	 \param ind The index to search for.
	 \return The index of the found edge, or -13 if no edge had the point as a removed point.
	 */
	int edge_w_removed(int ind)
	{
		
		for (unsigned int ii=0; ii<num_edges_; ii++){
			
			for (auto iter=edges_[ii].removed_begin(); iter!=edges_[ii].removed_end(); ++iter) {
				if (*iter == ind){
					return ii;
				}
			}
			
			
		}
		
		return -13;
	}
	
	
	/**
	 \brief  Find candidates for merging, by finding vertices with type NEW.
	 
	 \return a vector of integers containing the indices of edges to merge together into a single new edge.
	 \param V The Vertex set to which the Decomposition refers.
	 */
	std::vector<int> get_merge_candidate(const VertexSet & V);
	
	
	/**
	 \brief Merge method for curves, to remove edges with NEW type points.
	 
	 \todo Make this function callable even outside the interslice method.
	 
	 \param W_midpt Skeleton witness set coming from the precedin(containing) interslice call.
	 \param V the Vertex set to which this Decomposition refers.
	 \param pi_in The linear projection being used to decompose.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current solver configuration.
	 */
	void merge(WitnessSet & W_midpt, VertexSet & V,
			   vec_mp * pi_in,
			   BertiniRealConfig & program_options,
			   SolverConfiguration & solve_options);
	
	
	/**
	 \brief The main call for decomposing an indepentent curve.
	 
	 Includes deflation, etc.
	 
	 \param V The Vertex set which runs through everything.
	 \param W_curve The input witness set for the curve, computed probably by Bertini's tracktype:1.
	 \param pi A set of pointers to vec_mp's containing the projection being used.
	 \param program_options The current state of Bertini_real.
	 \param solve_options The current state of the solver routines.
	 */
	void main(VertexSet & V,
			  WitnessSet & W_curve,
			  vec_mp *pi,
			  BertiniRealConfig & program_options,
			  SolverConfiguration & solve_options);
	
	
	
	
	
	/**
	\brief the main function for computing the self-intersection of a non-self-conjugate dimension 1 component.
	 
	 \param W		the witness set
	 \param V		the Vertex set being passed around.  C indexes into here.
	 \param num_vars		the total number of variables in the problem.
	 \param program_options		main structure holding configuration
	 \param solve_options the current state of the solver.
	 */
	void 	computeCurveNotSelfConj(const WitnessSet		& W,
									VertexSet				&V,
									int						num_vars,
									BertiniRealConfig		& program_options,
									SolverConfiguration & solve_options);
	
	
	
	/**
	 \brief the main function for computing cell decom for a curve.  only for use on a self-conjugate component.
	 
	 \param W	WitnessSet containing linears, patches, points, etc.  much info and calculation performed from this little guy.
	 \param pi the set of projections to use.  in this case, there should be only 1.
	 \param V		Vertex set structure into which to place collected data.
	 \param options program configuration.
	 \param solve_options solver configuration.
	 */
	void computeCurveSelfConj(const WitnessSet & W,
							  vec_mp *pi,
							  VertexSet &V,
							  BertiniRealConfig & options,
							  SolverConfiguration & solve_options);
	
	
	
	
	/**
	 \brief The main method for slicing above critical points, halfway between them, and connecting the dots.
	 
	 This method peforms three functions.  It slices the curve halfway between the critical points, and tracks the midpoints-upstairs to the bounding critical values, added new points where necessary.  It also calls merge if desired.
	 
	 \return SUCCESSFUL
	 \param W_curve The main witness set for the curve.
	 \param W_crit_real Witness set containing the real critical points for the curve, previously computed.
	 \param pi The projection being used to decompose.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 \param V The Vertex set being used to contain the computed points.
	 */
	int interslice(const WitnessSet & W_curve,
				   const WitnessSet & W_crit_real,
				   vec_mp *pi,
				   BertiniRealConfig & program_options,
				   SolverConfiguration & solve_options,
				   VertexSet & V);
	
	
	
	/**
	 \brief Compute actual critical points of the curve, by regeneration and the left-nullspace method.
	 
	 \return SUCCESSFUL
	 \param W_curve Witness set for the curve, from which we will regenerate
	 \param pi the linear projection being used.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 \param W_crit_real The computed value, containing the real critical points of the curve.
	 */
	int compute_critical_points(const WitnessSet & W_curve,
								vec_mp *pi,
								BertiniRealConfig & program_options,
								SolverConfiguration & solve_options,
								WitnessSet & W_crit_real);
	
	
	
	
	
	/**
	 \brief Compute the intersection of the curve with a sphere containing the supplied critical points.
	 
	 Calling the sphere intersection solver, this method finds the points of intersection between the sphere and curve.
	 
	 \return SUCCESSFUL flag.
	 \param W_crit_real The previously computed real critical points.
	 \param W supplied witness set for the curve.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	int get_sphere_intersection_pts(WitnessSet *W_crit_real,
							   const WitnessSet & W,
							   BertiniRealConfig & program_options,
							   SolverConfiguration & solve_options);
	
	
	/**
	 \brief send a curve to a single MPI target
	 
	 \param target Who to send to.
	 \param mpi_config The state of MPI.
	 */
	void send(int target, ParallelismConfig & mpi_config);
	
	
	/**
	 \brief receive a curve from a single MPI source
	 
	 \param source Who to receive from.
	 \param mpi_config The state of MPI.
	 */
	void receive(int source, ParallelismConfig & mpi_config);
	
	
	
	
	
	
	/**
	 \brief Initialize a curve for sampling by the adaptive method.
	 
	 \see Curve::adaptive_sampler
	 
	 \ingroup samplermethods
	 
	 \return The number 0.
	 */
	int adaptive_set_initial_sample_data();
	
	
	
	/**
	 \brief Sample a curve using adaptive method based on convergence of the distance between next estimated point, and next computed point.
	 
	 \todo Add a better summary of the adaptive curve method.
	 
	 \ingroup samplermethods
	 
	 \param V the Vertex set containing the points of the Decomposition.
	 \param sampler_options The current state of the sampler program.
	 \param solve_options The current state of the solver and tracker configuration.
	 */
	void adaptive_sampler_movement(VertexSet &V,
								   sampler_configuration & sampler_options,
								   SolverConfiguration & solve_options);
	
	
	
	/**
	\brief Sample a curve using adaptive method based on distance between computed samples, and a maximum number of refinement passes.
	 
	 \todo Add a better summary of the adaptive curve method.
	
	 \ingroup samplermethods
	 
	 \param V the Vertex set containing the points of the Decomposition.
	 \param sampler_options The current state of the sampler program.
	 \param solve_options The current state of the solver and tracker configuration.
	*/
	void adaptive_sampler_distance(VertexSet &V,
						  sampler_configuration & sampler_options,
						  SolverConfiguration & solve_options);
	
	
	
	/**
	 \brief sets up refinement flags to YES for every interval, for first pass of adaptive refinement.
	 
	 \param num_refinements The number of intervals to refine.
	 \param refine_flags The mutable vector of bools indicating whether to refine a particular interval.
	 \param current_indices Indices of points between which to refine (or not).
	 \param V The Vertex set storing the points.
	 \param current_edge The index of which edge is being refined.
	 \param sampler_options The current state of the program.
	 */
	void adaptive_set_initial_refinement_flags(int & num_refinements,
											   std::vector<bool> & refine_flags,
											   std::vector<int> & current_indices,
											   VertexSet &V,
											   int current_edge, sampler_configuration & sampler_options);
	
	
	/**
	 \brief Initialize for the fixed-number curve sampler method.
	 
	 \param target_num_samples The number of samples per edge, including boundary points.
	 \return The number 0.
	 */
	int fixed_set_initial_sample_data(int target_num_samples);
	
	
	/**
	 \brief Sample a curve so it has an equal number of points per edge, including boundary points.
	 
	 \todo Add a description ofthe method with picture right here.
	 
	 \param V the Vertex set containing the points computed.
	 \param sampler_options The current state of the sampler program.
	 \param solve_options The current state of the solver.
	 \param target_num_samples The number of points to get on each edge, including boundary points.
	 */
	void fixed_sampler(VertexSet &V,
					   sampler_configuration & sampler_options,
					   SolverConfiguration & solve_options,
					   int target_num_samples);
	
	
	/**
	 \brief Dump the curve's sampler data to a file.
	 
	 \param samplingName The name of the file to write.
	 */
	void  output_sampling_data(boost::filesystem::path samplingName);
	
	
	
	
	
	
	
	/**
	 \brief Reset the curve to an empty set.
	 */
	void reset()
	{
		
		
		Decomposition::reset();
		
		this->clear();
		this->init();
		
		
	}
	
	
	/**
	 constructor
	 */
	Curve() : Decomposition()
	{
		this->init();
	}

	/**
	 
	assignment operator
	 */
	Curve & operator=(const Curve& other){
		this->init();
		this->copy(other);
		return *this;
	}
	
	
	/** 
	 copy constructor
	 */
	Curve(const Curve & other){
		this->init();
		this->copy(other);
	}
	
	
	/** 
	 destructor
	 */
	~Curve()
	{
		this->clear();
	}
	
	
	
	
protected:
	
	/**
	 clear the contents.  not a complete reset.
	 */
	void clear()
	{
		edges_.clear();
		num_edges_ = 0;
	}
	
	
	/**
	 get ready!
	 */
	void init(){
		
		num_edges_ = 0;
		set_dimension(1);
	}
	
	/**
	 copy the contents from one to the other.
	 */
	void copy(const Curve & other)
	{
		Decomposition::copy(other);
		this->edges_ = other.edges_;
		this->num_edges_ = other.num_edges_;
	}
	
}; // end Curve










/**
 writes the input file for the diagonal homotopy used in curve case to find the intersection points.
 
 \param outputFile the name of the file to create
 \param funcInputx	the name of the first input file
 \param funcInputy	the name of the second input file
 \param configInput	the name of the config file to use
 \param L			the linears for the problem
 \param num_vars		the number of variables in the problem, including the homogeneous ones.
 */
void 	diag_homotopy_input_file(boost::filesystem::path outputFile,
								 boost::filesystem::path funcInputx,
								 boost::filesystem::path funcInputy,
								 boost::filesystem::path configInput,
								 vec_mp L,
								 int   num_vars);




/**
 writes the start file for the diagonal homotopy used in curve case to find the intersection points.
 
 \param startFile the name of the start file to write
 \param W the input witness set
 */
void 	diag_homotopy_start_file(boost::filesystem::path startFile,
								 const WitnessSet & W);








/**
 \brief Check whether a projection is valid or not, by testing the rank of the jacobian at a random point.
 
 \return a boolean integer indicating whether the projection is valid.
 \param W witness set containing required information
 \param projection The linear projection being checked.
 \param solve_options The current state of the solver.
 */
int verify_projection_ok(const WitnessSet & W,
						 vec_mp * projection,
						 SolverConfiguration & solve_options);


/**
 \brief Check whether a projection is valid or not, by testing the rank of the jacobian at a random point.
 
 \return a boolean integer indicating whether the projection is valid.
 \param W witness set containing required information
 \param randomizer The way the system is randomized. 
 \param projection The linear projection being checked.
 \param solve_options The current state of the solver.
 */
int verify_projection_ok(const WitnessSet & W,
						 std::shared_ptr<SystemRandomizer> randomizer,
						 vec_mp * projection,
						 SolverConfiguration & solve_options);





#endif
