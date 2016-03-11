
#ifndef SURFACE_H
#define SURFACE_H


/** \file surface.hpp */

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


#include <boost/filesystem.hpp>

#include "bertini_headers.hpp"


#include "fileops.hpp"


#include "nullspace.hpp"

#include "solver_midpoint_tracker.hpp"

#include "checkSelfConjugate.hpp"
#include "curve.hpp"
#include "data_type.hpp"
#include "isosingular.hpp"
#include "programConfiguration.hpp"

#include "solver_sphere.hpp"



/**
 \brief the 2-Cell
 
 The Face data type.  Extends the functionality of the basic Cell type, adding left and right indices of connected edges, top and bottom indices and system names, left and right projection values, and an index for which midslice it came from.  The midpoint member is inherited from the base class Cell.
 */
class Face : public Cell
{
	
	std::vector<int>	left_edges_;  ///< index into vertices
	std::vector<int>	right_edges_; ///< index into vertices
	
	int top_edge_index_;			///<  index of the top edge in the appropriate Decomposition, as indicated by system_name_top.
	int bottom_edge_index_;			///< index of the bottom edge in the appropriate Decomposition, as indicated by system_name_bottom.
	
	
	std::string system_name_bottom_; ///< the plain text name of the bottom edge's system.  e.g. crit_curve
	std::string system_name_top_; ///< the plain text name of the top edge's system.  e.g. crit_curve
	
	int crit_slice_index_; ///< which midpoint slice this Face came from.
	
	
	
	comp_mp left_crit_val_; ///< the pi_0 projection value of the left crit slice
	comp_mp right_crit_val_; ///< the pi_0 projection value of the right crit slice
	
public:
	
	/**
	 \brief get a pointer to the comp_mp of the left critical value.
	 \return pointer to a comp_mp of the critical value of the left edge.
	 */
	comp_mp * left_crit_val()
	{
		return & left_crit_val_;
	}
	
	/**
	 \brief set the projection value of the left critical edge
	 \param new_left_crit_val the new value of the projection for left edge
	 */
	void set_left_crit_val(comp_mp new_left_crit_val)
	{
		set_mp(left_crit_val_,new_left_crit_val);
	}
	
	/**
	 \brief get a pointer to the comp_mp of the right critical value.
	 \return pointer to a comp_mp of the critical value of the right edge.
	 */
	comp_mp * right_crit_val()
	{
		return & right_crit_val_;
	}
	
	
	/**
	 \brief set the projection value of the right critical edge
	 \param new_right_crit_val the new value of the projection for right edge
	 */
	void set_right_crit_val(comp_mp new_right_crit_val)
	{
		set_mp(right_crit_val_,new_right_crit_val);
	}
	
	
	/**
	 \brief add an edge to set of found edges.
	 \param index the index of the edge
	 */
	void add_left_edge(int index)
	{
		left_edges_.push_back(index);
	}
	
	
	
	/**
	 \brief add an edge to set of found edges.
	 \param index the index of the edge
	 */
	void add_right_edge(int index)
	{
		right_edges_.push_back(index);
	}
	
	
	
	
	/**
	 \brief get the bottom edge index
	 \return the index of the bottom edge
	 */
	int bottom_edge() const
	{
		return bottom_edge_index_;
	}
	
	
	/**
	 \brief set the bottom edge index
	 \param index the new index to set
	 \return the index of the bottom edge
	 */
	int set_bottom_edge(int index)
	{
		return bottom_edge_index_ = index;
	}
	
	
	/**
	 \brief get the top edge index
	 \return the index of the top edge
	 */
	int top_edge() const
	{
		return top_edge_index_;
	}
	
	
	/**
	 \brief get the top edge index
	 \param index the new index
	 \return the index of the top edge
	 */
	int set_top_edge(int index)
	{
		return top_edge_index_ = index;
	}
	
	
	
	
	
	
	
	/**
	 \brief set the name of the top system
	 \param new_name the new name of the top system
	 */
	void system_name_top(std::string new_name)
	{
		system_name_top_ = new_name;
	}
	
	/**
	 \brief get the name of the top system
	 \return the string name of the top system
	 */
	std::string system_name_top() const
	{
		return system_name_top_;
	}
	
	
	
	/**
	 \brief set the name of the top system
	 \param new_name the new name of the bottom system
	 */
	void system_name_bottom(std::string new_name)
	{
		system_name_bottom_ = new_name;
	}
	
	/**
	 \brief get the name of the bottom system
	 \return the string name of the bottom system
	 */
	std::string system_name_bottom() const
	{
		return system_name_bottom_;
	}
	
	
	
	/**
	 \brief get the index of the crit slice to which this corresponds
	 \return the index of the critical slice to which this Cell corresponds
	 */
	int crit_slice_index() const
	{
		return crit_slice_index_;
	}
	
	
	/**
	 \brief set the index of the crit slice to which this corresponds
	 \param new_index the new index
	 \return the index of the critical slice to which this Cell corresponds
	 */
	int crit_slice_index(int new_index)
	{
		return crit_slice_index_ = new_index;
	}
	
	
	/**
	 \brief get the index of the left edge
	 \return the index
	 \param index the index to look up.
	 */
	int left_edge(unsigned int index) const
	{
		if (index>=left_edges_.size()) {
			throw std::out_of_range("trying to access left edge, index out of range");
			return -1;
		}
		else
		{
			return left_edges_[index];
		}
	}

	
	/**
	 \brief get the index of the right edge
	 \return the index
	 \param index the index to look up.
	 \throws out_of_range, if the index exceeds the size of the right edges
	 */
	int right_edge(unsigned int index) const
	{
		if (index>=right_edges_.size()) {
			throw std::out_of_range("trying to access right edge, index out of range");
			return -1;
		}
		else
		{
			return right_edges_[index];
		}
	}
	
	
	
	
	
	/**
	 \brief get the number of left edges
	 \return the number of left edges
	 */
	unsigned int num_left() const
	{ return left_edges_.size();	 ///< the number of left mapped edges.
	}
	
	/**
	 \brief get the number of right edges
	 \return the number of right edges
	 */
	unsigned int num_right() const
	{ return right_edges_.size();	 ///< the number of right mapped edges.
	}
	
	
	

	
	

	

	
	friend std::ostream & operator<<(std::ostream &os, const Face & f)
	{
		os << f.midpt() << std::endl;
		os << f.crit_slice_index_ << std::endl << f.top_edge_index_ << " " << f.bottom_edge_index_ << std::endl;
		os << f.system_name_top_ << " " << f.system_name_bottom_ << std::endl;
		
		os << f.num_left() << std::endl;
		for (auto jj=f.left_edges_.begin(); jj!=f.left_edges_.end(); jj++) {
			os << *jj << " ";
		}
		os << std::endl;
		
		os << f.num_right() << std::endl;
		for (auto jj=f.right_edges_.begin(); jj!=f.right_edges_.end(); jj++) {
			os << *jj << " ";
		}
		os << std::endl << std::endl;
		
		return os;
	}
	
	
	friend std::istream & operator>>(std::istream &os, Face & f)
	{
		
		f.read_from_stream(os);
		return os;
	}
	
	/**
	 \brief read a Face from an input stream.
	 
	 \param os the istream input stream to pass into this Face.
	 */
	virtual void read_from_stream( std::istream &os )
	{
		
		
		int tmp;
		os >> tmp;  midpt(tmp);
		os >> crit_slice_index_ >> top_edge_index_ >> bottom_edge_index_;
		os >> system_name_top_ >> system_name_bottom_;
		
		unsigned int temp_num_left;
		os >> temp_num_left;
		left_edges_.resize(temp_num_left);
		for (unsigned int jj=0; jj<temp_num_left; jj++) {
			os >> left_edges_[jj];
		}
		

		unsigned temp_num_right;
		os >> temp_num_right;
		right_edges_.resize(temp_num_right);
		for (unsigned int jj=0; jj<temp_num_right; jj++) {
			os >> right_edges_[jj];
		}
	}
	
	
	Face() : Cell()
	{
		init();
	}
	
	~Face(){
		clear_mp(left_crit_val_);
		clear_mp(right_crit_val_);
	}
	
	Face(const Face & other){
		init();
		copy(other);
	}
	
	Face& operator=(const Face & other){

		init();
		
		copy(other);
		return *this;
	}
	
	
	/**
	 \brief single-target MPI send of a Face.
	 
	 
	 \param target The ID of the MPI target to send the Face to.
	 \param mpi_config The current state of MPI
	 */
	void send(int target, ParallelismConfig & mpi_config);
	
	
	
	/**
	 \brief single-source MPI receive of a Face.
	 
	 \param source The ID of the MPI source from whom to receive a Face.
	 \param mpi_config The current state of MPI
	 */
	void receive(int source, ParallelismConfig & mpi_config);
	
	/**
	 \brief Test a Face for degeneracy, based on which crit slice the Face is from.
	 
	 \return Boolean, true if degenerate (crit_slice_index < 0), false if not.
	 */
	bool is_degenerate()
	{
		if (crit_slice_index_<0) {
			return true;
		}
		else{
			return false;
		}
	}
	
private:
	
	void init()
	{
		system_name_top_ = "UNSET_TOP";
		system_name_bottom_ = "UNSET_BOTTOM";
		
		
		init_mp2(left_crit_val_,1024);
		init_mp2(right_crit_val_,1024);
		
		left_edges_.resize(0);
		right_edges_.resize(0);
		top_edge_index_ = bottom_edge_index_ = -1;
		crit_slice_index_ = -1;

		
	}
	
	
	
	void copy(const Face & other)
	{

				
		Cell::copy(other);
		
		this->system_name_bottom_ = other.system_name_bottom_;
		this->system_name_top_ = other.system_name_top_;
		
		
		this->crit_slice_index_ = other.crit_slice_index_;
		
		this->left_edges_.clear();
		for (auto ii = other.left_edges_.begin(); ii!=other.left_edges_.end(); ii++) {
			this->left_edges_.push_back(*ii);
		}
		
		this->right_edges_.clear();
		for (auto ii = other.right_edges_.begin(); ii!=other.right_edges_.end(); ii++) {
			this->right_edges_.push_back(*ii);
		}
		
		this->top_edge_index_ = other.top_edge_index_;
		this->bottom_edge_index_ = other.bottom_edge_index_;
		
		
		set_mp(this->left_crit_val_, other.left_crit_val_);
		set_mp(this->right_crit_val_, other.right_crit_val_);
	}
	
};











/**
 \brief A simple Triangle class, holding three integers referring to points in a vertex set.
 */
class Triangle
{
	long long v1_; ///< point index 1.
	long long v2_; ///< point index 2.
	long long v3_; ///< point index 3.
public:
	
	void set_v(long long new_v1, long long new_v2, long long new_v3)
	{
		v1_ = new_v1;
		v2_ = new_v2;
		v3_ = new_v3;
		
	}
	
	void set_v1(long long new_v){
		v1_ = new_v;
	}
	
	void set_v2(long long new_v){
		v2_ = new_v;
	}
	
	void set_v3(long long new_v){
		v3_ = new_v;
	}
	
	
	inline long long v1() const
	{return v1_;}

	inline long long v2() const
	{return v2_;}
	
	inline long long v3() const
	{return v3_;}
	
	Triangle(){};
	
	Triangle(long long new_v1,long long new_v2,long long new_v3)
	{
		v1_ = new_v1;
		v2_ = new_v2;
		v3_ = new_v3;
		
	}
	
	
	friend std::ostream & operator<<(std::ostream &os, const Triangle & t)
	{
		
		
		os << t.v1_ << " " << t.v2_ << " " << t.v3_;
		
		return os;
	}
	
	
	/**
	 \brief test a Triangle for degeneracy.
	 
	 \return true if degenerate, false if not.
	 */
	bool is_degenerate()
	{
		if (v1_==v2_) {
			return true;
		}
		
		if (v2_==v3_) {
			return true;
		}
		
		if (v1_==v3_) {
			return true;
		}
		
		return false;
	}
};





class SingularObjectMetadata
{
	unsigned int multiplicity_;
	unsigned int index_;
	
public:
	
	
	SingularObjectMetadata(){multiplicity_ = index_ = 0;}
	
	/**
	 constructor for setting both multiplicity and index
	 \param mult the multiplicity of the object
	 \param ind the index of the object.  this should start at 0, in principle.
	 */
	SingularObjectMetadata(unsigned int mult, unsigned int ind)
	{
		multiplicity_ = mult;
		index_ = ind;
	}
	
	/**
	 get the multiplicity
	 \return the multiplicity
	 */
	unsigned int multiplicity() const
	{return multiplicity_;}
	
	/**
	 \brief set the multiplicity
	 \param new_multiplicity the new multiplicity
	 */
	void set_multiplicity(unsigned int new_multiplicity){multiplicity_ = new_multiplicity;};
	
	
	
	/**
	 get the index
	 \return the multiplicity
	 */
	unsigned int index() const
	{return index_;}
	
	
	/**
	 \brief set the index
	 \param new_index the new index
	 */
	void set_index(unsigned int new_index){index_ = new_index;};
	
	
	/**
	 basic equality comparitor.  both multiplicity and index must be the same.
	 \return true if multiplicity and index are == for rhs and lhs
	 \param lhs left hand side of comparison
	 \param rhs right hand side of comparison
	 */
	friend bool operator==(const SingularObjectMetadata & lhs, const SingularObjectMetadata & rhs)
	{
		if (lhs.multiplicity_==rhs.multiplicity_ && lhs.index_==rhs.index_) {
			return true;
		}
		else{
			return false;}
	}
	
	
	/**
	 basic inequality comparitor.  either multiplicity or index must be the different.
	 \return true if multiplicity or index are != for rhs and lhs.  false if both ==.
	 \param lhs left hand side of comparison
	 \param rhs right hand side of comparison
	 */
	friend bool operator!=(const SingularObjectMetadata & lhs, const SingularObjectMetadata & rhs)
	{
		if (lhs.multiplicity_!=rhs.multiplicity_ || lhs.index_!=rhs.index_) {
			return true;
		}
		else{
			return false;}
	}

	
	
	
	/**
	 basic less than comparitor.
	 \return true if lhs.multiplicity<rhs.multiplicity.  false if rhs.multiplicity<lhs.multiplicity.  otherwise, multiplicities are equal, and we sort on index.  true if lhs.index<rhs.index, else false.
	 \param lhs left hand side of comparison
	 \param rhs right hand side of comparison
	 */
	friend bool operator<(const SingularObjectMetadata & lhs, const SingularObjectMetadata & rhs)
	{
		if (lhs.multiplicity_<rhs.multiplicity_){
			return true;
		}
		else if (rhs.multiplicity_<lhs.multiplicity_){
			return false;
		}
		else if (lhs.index_<rhs.index_)
		{
			return true;
		}
		else{
			return false;
		}
	}

	
	
	/**
	 basic less than or equal to comparitor.
	 \return !(rhs<lhs)
	 \param lhs left hand side of comparison
	 \param rhs right hand side of comparison
	 */
	friend bool operator<=(const SingularObjectMetadata & lhs, const SingularObjectMetadata & rhs)
	{
		return !(rhs<lhs);
	}

	
	/**
	 basic greater than comparitor.
	 \return (rhs<lhs)
	 \param lhs left hand side of comparison
	 \param rhs right hand side of comparison
	 */
	friend bool operator>(const SingularObjectMetadata & lhs, const SingularObjectMetadata & rhs)
	{
		return rhs<lhs;
	}

	
	
	/**
	 basic greater than or equal to comparitor.
	 \return !(lhs<rhs)
	 \param lhs left hand side of comparison
	 \param rhs right hand side of comparison
	 */
	friend bool operator>=(const SingularObjectMetadata & lhs, const SingularObjectMetadata & rhs)
	{
		return  !(lhs < rhs);
	}

	
	
	
	
};


/**
 \brief Bertini_real surface Decomposition.
 
 includes methods to add vertices, look up vertices, etc
 */
class Surface : public Decomposition
{
	
	std::vector<Face> faces_;
	//these counters keep track of the number of things
	
	unsigned int      num_faces_; ///< how many faces are there
	
	
	std::vector<int> num_samples_each_face_; ///< the number of samples per Face
	
	std::vector< std::vector< Triangle > > samples_; ///< refined triangulation of the surface.
	
	
	
	std::vector< Curve > mid_slices_; ///< the mid slice curves, occurring halfway between critical points of critical curves.
	std::vector< Curve > crit_slices_; ///< critical slice curves, occurring at critial points of critical curves.
	Curve crit_curve_; ///< the main critical curve, being nonsingular.
	Curve sphere_curve_; ///< the intersection of the surface and a sphere.
	
	std::map<SingularObjectMetadata,Curve> singular_curves_; ///<  the singular curves, which are formally but not really part of the critical curve.
	size_t num_singular_curves_; ///< how many singular curves are there.
	
	
	
	
public:
	
	/**
	 \brief get the number of faces in the surface.
	 
	 \return the number of faces.
	 */
	unsigned int num_faces()
	{
		return num_faces_;
	}
	
	
	/**
	 \brief query the number of singular curves in the Decomposition.
	 \return the number of singular curves in the Decomposition.
	 */
	unsigned int num_singular_curves() const
	{
		return num_singular_curves_;
	}
	
	
	const std::map<SingularObjectMetadata,Curve>::const_iterator singular_curves_iter_begin() const
	{
		return singular_curves_.begin();
	}
	
	
	
	const std::map<SingularObjectMetadata,Curve>::const_iterator singular_curves_iter_end() const
	{
		return singular_curves_.end();
	}
	
	
	std::map<SingularObjectMetadata,Curve>::iterator singular_curves_iter_begin()
	{
		return singular_curves_.begin();
	}
	
	
	
	std::map<SingularObjectMetadata,Curve>::iterator singular_curves_iter_end()
	{
		return singular_curves_.end();
	}
	
	
	
	
	std::vector< Curve >::iterator mid_slices_iter_begin()
	{
		return mid_slices_.begin();
	}
	
	std::vector< Curve >::iterator mid_slices_iter_end()
	{
		return mid_slices_.end();
	}
	
	
	
	std::vector< Curve >::iterator crit_slices_iter_begin()
	{
		return crit_slices_.begin();
	}
	
	std::vector< Curve >::iterator crit_slices_iter_end()
	{
		return crit_slices_.end();
	}
	
	
	
	/**
	 \brief get a const reference to the critical curve.
	 \return a reference to the critical curve Decomposition.
	 */
	const Curve & crit_curve() const
	{return crit_curve_;}
	
	/**
	 \brief get a reference to the critical curve.
	 \return a mutable reference to the critical curve Decomposition.
	 */
	Curve & crit_curve()
	{return crit_curve_;}
	
	
	
	
	
	
	/**
	 \brief get a const reference to the sphere curve.
	 \return a reference to the sphere curve Decomposition.
	 */
	const Curve & sphere_curve() const
	{return sphere_curve_;}
	
	
	
	/**
	 \brief get a mutable reference to the sphere curve.
	 \return a reference to the sphere curve Decomposition.
	 */
	Curve & sphere_curve()
	{return sphere_curve_;}
	
	
	/**
	 \brief set up the faces in a surface Decomposition from a file
	 
	 \param load_from_me the name of the file to read.
	 */
	void read_faces(boost::filesystem::path load_from_me);
	
	
	/**
	 \brief set up a Surface from a set of files
	 
	 \param base the folder in which the Decomposition lives.
	 */
	void setup(boost::filesystem::path base);
	
	
	/**
	 \brief print a surface Decomposition to a set of files located in folder base.
	 
	 \param base the name of the folder to which to print the Decomposition.
	 */
	void print(boost::filesystem::path base);
	
	/**
	 \brief print all the faces in a Decomposition to a file named outputfile.
	 
	 \param outputfile the name of the file into which to dump the faces (and a bit of header).
	 */
	void print_faces(boost::filesystem::path outputfile);
	
	

	

	/**
	 \brief add a Face after it has been constructed.
	 
	 also increments a counter.
	 
	 \param F the constructed Face to add to the Decomposition.
	 */
	void add_face(const Face & F)
	{
		faces_.push_back(F);
		num_faces_++;
	}
	

	
	

	
	/**
	 \brief The main outer method for decomposing a surface in Bertini_real.
	 
	 \param V the VertexSet into which to deposit the found points, along with their metadata.
	 \param W the input witness set for the surface, complete with points, linears, and patches.
	 \param pi the two linear projection vectors to use.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	void main(VertexSet & V,
						const WitnessSet & W,
						vec_mp *pi,
						BertiniRealConfig & program_options,
						SolverConfiguration & solve_options);
	
	
	/**
	\brief Perform a few preliminaries for the surface, including deflation, parsing of the input file, and creation of a randomization method.
	 
	 \param W_surf the witness set for the surface.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	*/
	void beginning_stuff(const WitnessSet & W_surf,
						BertiniRealConfig & program_options,
						 SolverConfiguration & solve_options);
	
	
	/**
	 \brief deflate the input system for the surface with respect to computed critical witness points, and split the points according to deflation system.
	 
	 \param split_sets an output variable, contains the witness sets split according to deflation system.
	 \param higher_multiplicity_witness_sets The witness sets computed by the nullspace method for critical conditions, and split according to 'multiplicity', or the number of times they appeared as solutions to the nullspace system.
	 \param program_options The current state of Bertini_real.
	 \param solve_options The current state of the solver.
	 */
	void deflate_and_split(std::map< SingularObjectMetadata, WitnessSet > & split_sets,
						   std::map<int, WitnessSet > & higher_multiplicity_witness_sets,
						   WitnessSet & points_which_needed_no_deflation,
						   BertiniRealConfig & program_options,
						   SolverConfiguration & solve_options);
	
	/**
	 \brief compute the critical points of all singular curves.
	 
	 \see deflate_and_split
	 
	 \param W_singular_crit A return parameter, contains the computed critical points.
	 \param split_sets The split witness sets for the critical curves, produced by deflate_and_split
	 \param V the total vertex set into which to deposit the computed points on the surfaces.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the Solver.
	 */
	void compute_singular_crit(WitnessSet & W_singular_crit,
							   const std::map<SingularObjectMetadata, WitnessSet> & split_sets,
							   VertexSet & V,
							   BertiniRealConfig & program_options,
							   SolverConfiguration & solve_options);
	
	
	/**
	 \brief decompose all singular curves.
	 
	 using the computed witness points and critical points, decompose the singular curves.
	 
	 \param W_total_crit the computed critical points of all critical curves, including the critical curve itself, the sphere curve, and the singular curves.
	 \param split_sets The witness sets for the singular curves, produced by nullspace left.
	 \param V the vertex set into which the points are stored.
	 \param program_options The current state of Bertini_real.
	 \param solve_options The current state of the solver.
	 */
	void compute_singular_curves(const WitnessSet & W_total_crit,
								 const std::map< SingularObjectMetadata, WitnessSet> & split_sets,
								 VertexSet & V,
								 BertiniRealConfig & program_options,
								 SolverConfiguration & solve_options);
	
	
	
	/**
	 \brief slice the surface at the supplied projection values downstairs.
	 
	 \param W_surf the witness set for the surface itself.
	 \param V the vertex set into which the computed points are committed.
	 \param projection_values_downstairs The pi_0 projection values of the critical points.
	 \param slices the computed slice curve decompositions.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 \param kindofslice A string indicating what kind of slice you are decompositon -- purely for screen output.
	 */
	void compute_slices(const WitnessSet W_surf,
											VertexSet & V,
											vec_mp projection_values_downstairs,
											std::vector< Curve > & slices,
											BertiniRealConfig & program_options,
											SolverConfiguration & solve_options,
											std::string kindofslice);
	
	
    /**
	 \brief compute a witness set for the critical curve(s), including the critcurve itself, as well as singular curves.
	 
	 \param W_critcurve One of two computed parameters for this functions, this one contains regular witness points for the deflation-free critcurve.
	 \param higher_multiplicity_witness_sets The other computed returned parameter to this function, this is a map between multiplicities and witness sets of that multiplicity.  These must be deflated.
	 \param W_surf The input witness set, for the surface.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
    void compute_critcurve_witness_set(WitnessSet & W_critcurve,
									   std::map<int, WitnessSet> & higher_multiplicity_witness_sets,
                                       const WitnessSet & W_surf,
                                       BertiniRealConfig & program_options,
                                       SolverConfiguration & solve_options);
    
	
	/**
	 \brief compute critical points of the regular critical curve.
	 
	 \param W_critcurve_crit The computed value, containing the critical points of the critical curve.
	 \param W_critcurve Witness set for the critical curve.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
    void compute_critcurve_critpts(WitnessSet & W_critcurve_crit, // the computed value
                                   WitnessSet & W_critcurve,
                                   BertiniRealConfig & program_options,
                                   SolverConfiguration & solve_options);
    
    
	/**
	 \brief Decompose the critical curve with respect to pi_0.
	 
	 Assumes have the critical points and witness points already computed.  This essentially performs interslice.
	 
	 \see Curve::interslice
	 
	 
	 \param W_critcurve Witness set for the critical curve.
	 \param W_critpts The computed critical point of all critical curves including the singular curves.
	 \param V The vertex set into which to deposit the computed points.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	void compute_critical_curve(const WitnessSet & W_critcurve,
                                const WitnessSet & W_critpts,
                                VertexSet & V,
                                BertiniRealConfig & program_options,
                                SolverConfiguration & solve_options);
	
	
	

	
	
	/**
	 \brief Obtain a witness set for the sphere intersection curve.
	 
	 
	 \param W_surf Witness set for the surface.
	 \param W_intersection_sphere The computed value, returning the witness set for the sphere intersection curve.
	 \param program_options The current state of Bertini_real.
	 \param solve_options The current state of the solver.
	 */
	void compute_sphere_witness_set(const WitnessSet & W_surf,
									WitnessSet & W_intersection_sphere,//
									BertiniRealConfig & program_options,//
									SolverConfiguration & solve_options);
	
	
	/**
	 \brief Compute critical points of the sphere curve, from the previously computed witness points.
	 
	 \see Surface::compute_sphere_witness_set
	 
	 \param W_intersection_sphere The witness set for the curve, computed by compute_sphere_witness_set
	 \param W_sphere_crit The computed returned value, the critical points of the sphere curve.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	void compute_sphere_crit(const WitnessSet & W_intersection_sphere,
													 WitnessSet & W_sphere_crit,
													 BertiniRealConfig & program_options,
													 SolverConfiguration & solve_options);
	

	/**
	 \brief Performs interslice on the sphere curve, using the critical points and the witness points.
	 
	 \see Surface::compute_sphere_crit
	 \see Surface::compute_sphere_witness_set
	 
	 
	 \param W_intersection_sphere Input witness set for the curve.
	 \param W_crit The input critical points.
	 \param V The vertex set into which to commit the computed points.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	void compute_bounding_sphere(const WitnessSet & W_intersection_sphere,
															 const WitnessSet & W_crit,
															 VertexSet & V,
															 BertiniRealConfig & program_options,
															SolverConfiguration & solve_options);
	
	
	/**
	 \brief Main outer method for obtaining faces of the Decomposition.
	 
	 Switches between serial_connect and master_connect depending on parallel mode.  Produces the faces of the Decomposition by calling the midpoint tracker.
	 
	 \todo Remove pi as an input for this function.
	 
	 \param V The vertex set into which all the points have been collected.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	void connect_the_dots(VertexSet & V,
						  BertiniRealConfig & program_options,
						  SolverConfiguration & solve_options);
	
	
	/**
	 \brief The serial mode for producing the faces. 
	 
	 In a loop, calls make_face.
	 
	 \param V the vertex set into which all the decompositions index.
	 \param md_config The already-set-up midpoint tracker config object.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	void serial_connect(VertexSet & V,
						MidpointConfiguration & md_config,
						SolverConfiguration & solve_options,
						BertiniRealConfig & program_options);
	
	
	/**
	 \brief The master parallel mode for producing the faces.
	 
	 In a loop, sends workers indices of faces to make, and receives the completed faces.  Then adds the Face to this Decomposition.
	 
	 \param V the vertex set into which all the decompositions index.
	 \param md_config The already-set-up midpoint tracker config object.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	void master_connect(VertexSet & V,
						MidpointConfiguration & md_config,
						SolverConfiguration & solve_options,
						BertiniRealConfig & program_options);
	
	/**
	 \brief The worker mode for producing faces.
	 
	 Worker receives the md_config, surface to be connected, and vertex set from the master, then waits for an index to connect up.  It connects, and then sends back the completed Face to the master.
	 
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	void worker_connect(SolverConfiguration & solve_options,
						BertiniRealConfig & program_options);
	
	
	/**
	 \brief Send the next worker the index of the Face to make.
	 
	 \param ii which critslice the Face will correspond to.
	 \param jj which edge on the critslice the Face will correspond to.
	 \param next_worker The MPI ID of the worker to request to make the Face.
	 \param mpi_config the current state of MPI.
	 */
	void master_face_requester(int ii, int jj, int next_worker, ParallelismConfig & mpi_config);
	
	
	/**
	 \brief As a worker, receive the indices of the Face to make.
	 
	 \param ii which critslice the Face will correspond to.
	 \param jj which edge on the critslice the Face will correspond to.
	 \param mpi_config the current state of MPI.
	 */
	void worker_face_requester(int & ii, int & jj, ParallelismConfig & mpi_config);
	
	
	/**
	 \brief  Actually make the Face for the ii-th midslice and the jj-th edge, calling the midtrack solver.
	 
	 \return The constructed Face (which may be degenerate).
	 \param ii which critslice the Face will correspond to.
	 \param jj which edge on the critslice the Face will correspond to.
	 \param V the vertex set containing the vertices for the total Decomposition.
	 \param md_config The already-set-up midtrack config.
	 \param solve_options The current state of the solver.
	 \param program_options The current state of Bertini_real.
	 */
	Face make_face(int ii, int jj, VertexSet & V,
				   MidpointConfiguration & md_config,
				   SolverConfiguration & solve_options, BertiniRealConfig & program_options);
	
	
	
	
	
	
	/**
	 \brief Method for sampling a surface so that each Face contains approximately the same number of triangles.
	 
	 \todo Parallelize the sampler.
	 
	 \param V The vertex set holding the vertices.
	 \param sampler_options The current state of the sampler program.
	 \param solve_options The current state of the solver.
	 */
	void fixed_sampler(VertexSet &V,
					   sampler_configuration & sampler_options,
					   SolverConfiguration & solve_options);
	
	
	
	
	/**
	 \brief Write the results of a sampling run to a folder.
	 
	 \param base_path The name of the folder into which to write the sampling data.
	 */
	void  output_sampling_data(boost::filesystem::path base_path);
	
	
	
	
	
	
	Surface() : Decomposition()
	{
		this->init();
	}
	
	Surface & operator=(const Surface& other){
		this->init();
		this->copy(other);
		return *this;
	}
	
	Surface(const Surface & other){
		this->init();
		this->copy(other);
	}
	
	
	~Surface()
	{
		this->clear();
	}
	
	
	
	/**
	 \brief Find a curve with the input name as its input_filename.
	 
	 \return A pointer to the curve with the name, or NULL if it is not found.
	 \param findme The name you want to appear as the input_filename of the curve.
	 */
	Curve * curve_with_name(const std::string & findme)
	{
		
		if (findme.compare(crit_curve_.input_filename().filename().string())==0) {
			return &crit_curve_;
		}
		
		if (findme.compare(sphere_curve_.input_filename().filename().string())==0) {
			return &sphere_curve_;
		}
		for (auto iter = singular_curves_.begin(); iter!=singular_curves_.end(); ++iter) {
			if (findme.compare(iter->second.input_filename().filename().string())==0) {
				return &(iter->second);
			}
		}
		
		
		for (auto iter = mid_slices_.begin(); iter!=mid_slices_.end(); ++iter) {
			if (findme.compare(iter->input_filename().filename().string())==0) {
				return &(*iter);
			}
		}
		
		
		for (auto iter = crit_slices_.begin(); iter!=crit_slices_.end(); ++iter) {
			if (findme.compare(iter->input_filename().filename().string())==0) {
				return &(*iter);
			}
		}
	
		std::cout << "failed to find curve with name " << findme << std::endl;
				
		return NULL;
	}
	
	
	
	/**
	 \brief Single-target MPI send function for a Surface..
	 
	 \param target The MPI ID of the process to send it to.
	 \param mpi_config The current state of MPI.
	 */
	void send(int target, ParallelismConfig & mpi_config);
	
	/**
	 \brief Single-source MPI receive function.
	 
	 \param source The MPI ID of the person from whom to receive the Surface.
	 \param mpi_config The current state of MPI.
	 */
	void receive(int source, ParallelismConfig & mpi_config);
	
protected:
	
	
	void copy(const Surface & other)
	{
		Decomposition::copy(other);
		
		this->faces_ = other.faces_;
		
		this->num_faces_ = other.num_faces_;
		
		this->mid_slices_ = other.mid_slices_;
		this->crit_slices_ = other.crit_slices_;
		this->crit_curve_ = other.crit_curve_;
		this->singular_curves_ = other.singular_curves_;
		this->num_singular_curves_ = other.num_singular_curves_;
	}
	
	
	void init()
	{
		num_singular_curves_ = 0;
		num_faces_ = 0;
		set_dimension(2);
	}
	
	void clear()
	{
		
	}
};





//
/**
 
 will compute a randomizer matrix since you don't provide one. must have current PPD in solve_options for this to work correctly
 Assumes the input file for W is already parsed.
 
 \return SUCCESSFUL, unless W has no points, in which case returns TOLERABLE_FAILURE.
 \param W_match A computed value, this contains all the points in W which satisfy its input file.
 \param W_reject A computed value, this contains all points in W which do NOT satisfy its input file.
 \param W Input witness set, which you want to split in terms of which points do/not satisfy its input file (which was produced by isosingular deflation).
 \param solve_options The current state of the solver.
 */
int find_matching_singular_witness_points(WitnessSet & W_match,
										  WitnessSet & W_reject,
										  const WitnessSet & W,
										  SolverConfiguration & solve_options);


/**
 \brief Write a Bertini input file which slices at input linears.
 
 \param input_file The name of the input file to which we are appending functions.
 \param output_file The name of the output file, which will have the sliced system in it.
 \param linears The linear functions to slice with.
 \param num_to_add The number of linears we are adding on.  This should match the number of vec_mp in linears.
 \param W The input witness set.
 */
void create_sliced_system(boost::filesystem::path input_file, boost::filesystem::path output_file,
													vec_mp * linears, int num_to_add,
													const WitnessSet & W);



/**
 \brief Write a Bertini input file corresponding to the intersection of the input system with a sphere of a given center and radius.
 
 \param input_file The name of the input file to which we are appending functions.
 \param output_file The name of the output file, which will have the sphere system in it.
 \param sphere_radius The radius of the sphere.
 \param sphere_center The center of the sphere
 \param W the input witness set.
 */
void create_sphere_system(boost::filesystem::path input_file, boost::filesystem::path output_file,
													comp_mp sphere_radius,
													vec_mp sphere_center,
													const WitnessSet & W);


/**
 \brief Read the file named INfile, and get the S.surf information from it.
 
 
 \see Surface::print
 
 This function is used when reading a Decomposition back in from files.
 
 \param singular_multiplicities A returned value, containing the multiplicity information for the Decomposition.
 \param temp_num_mid The number of midslices.
 \param temp_num_crit The number of critslices.  Should be temp_num_mid+1.
 \param INfile The name of the file to read.
 */
void read_summary(std::vector<SingularObjectMetadata > & singular_multiplicities,
				  int & temp_num_mid,
				  int & temp_num_crit,
				  boost::filesystem::path INfile);


#endif


