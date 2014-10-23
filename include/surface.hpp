
#ifndef SURFACE_H
#define SURFACE_H


/** \file surface.hpp */


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


#include "nullspace_left.hpp"

#include "solver_midpoint_tracker.hpp"

#include "checkSelfConjugate.hpp"
#include "curve.hpp"
#include "data_type.hpp"
#include "isosingular.hpp"
#include "programConfiguration.hpp"

#include "solver_sphere.hpp"



/**
 \brief the 2-cell
 
 The face data type.  Extends the functionality of the basic cell type, adding left and right indices of connected edges, top and bottom indices and system names, left and right projection values, and an index for which midslice it came from.  The midpoint member is inherited from the base class cell.
 */
class face : public cell
{
public:
	
  std::vector<int>	left;  ///< index into vertices
  std::vector<int>	right; ///< index into vertices
	
	int top;			///<  index of the top edge in the appropriate decomposition, as indicated by system_name_top.
	int bottom;			///< index of the bottom edge in the appropriate decomposition, as indicated by system_name_bottom.
	
	unsigned int num_left;		///< the number of left mapped edges.
	unsigned int num_right;	 ///< the number of right mapped edges.
	
	
	comp_mp left_crit_val; ///< the pi_0 projection value of the left crit slice
	comp_mp right_crit_val; ///< the pi_0 projection value of the right crit slice
	
	std::string system_name_bottom; ///< the plain text name of the bottom edge's system.  e.g. crit_curve
	std::string system_name_top; ///< the plain text name of the top edge's system.  e.g. crit_curve

	
	int crit_slice_index; ///< which midpoint slice this face came from.
	
	
	friend std::ostream & operator<<(std::ostream &os, const face & f)
	{
		os << f.midpt() << std::endl;
		os << f.crit_slice_index << std::endl << f.top << " " << f.bottom << std::endl;
		os << f.system_name_top << " " << f.system_name_bottom << std::endl;
		
		os << f.left.size() << std::endl;
		for (int jj=0; jj<int(f.left.size()); jj++) {
			os << f.left[jj] << " ";
		}
		os << std::endl;
		
		os << f.right.size() << std::endl;
		for (int jj=0; jj<int(f.right.size()); jj++) {
			os << f.right[jj] << " ";
		}
		os << std::endl << std::endl;
		
		return os;
	}
	
	
	friend std::istream & operator>>(std::istream &os, face & f)
	{
		
		f.read_from_stream(os);
		return os;
	}
	
	/**
	 \brief read a face from an input stream.
	 
	 \param os the istream input stream to pass into this face.
	 */
	virtual void read_from_stream( std::istream &os )
	{
		int tmp;
		os >> tmp;  midpt(tmp);
		os >> crit_slice_index >> top >> bottom;
		os >> system_name_top >> system_name_bottom;
		
		os >> num_left;
		left.resize(num_left);
		for (unsigned int jj=0; jj<num_left; jj++) {
			os >> left[jj];
		}
		

		os >> num_right;
		right.resize(num_right);
		for (unsigned int jj=0; jj<num_right; jj++) {
			os >> right[jj];
		}
	}
	
	
	face() : cell()
	{
		init();
	}
	
	~face(){
		clear_mp(left_crit_val);
		clear_mp(right_crit_val);
	}
	
	face(const face & other){
		init();
		copy(other);
	}
	
	face& operator=(const face & other){

		init();
		
		copy(other);
		return *this;
	}
	
	
	/**
	 \brief single-target MPI send of a face.
	 
	 
	 \param target The ID of the MPI target to send the face to.
	 \param mpi_config The current state of MPI
	 */
	void send(int target, parallelism_config & mpi_config);
	
	/**
	 \brief single-source MPI receive of a face.
	 
	 \param source The ID of the MPI source from whom to receive a face.
	 \param mpi_config The current state of MPI
	 */
	void receive(int source, parallelism_config & mpi_config);
	
	/**
	 \brief Test a face for degeneracy, based on which crit slice the face is from.
	 
	 \return Boolean, true if degenerate (crit_slice_index < 0), false if not.
	 */
	bool is_degenerate()
	{
		if (crit_slice_index<0) {
			return true;
		}
		else{
			return false;
		}
	}
	
private:
	
	
	
	void init()
	{
		system_name_top = "UNSET_TOP";
		system_name_bottom = "UNSET_BOTTOM";
		
		
		init_mp2(left_crit_val,1024);
		init_mp2(right_crit_val,1024);
		
		num_left = num_right = 0;
		top = bottom = -1;
		crit_slice_index = -1;
	}
	
	
	void copy(const face & other)
	{
		cell::copy(other);
		
		this->system_name_bottom = other.system_name_bottom;
		this->system_name_top = other.system_name_top;
		
		
		this->crit_slice_index = other.crit_slice_index;
		
		this->left = other.left;
		this->right = other.right;
		
		this->top = other.top;
		this->bottom = other.bottom;
		
		this->num_left = other.num_left;
		this->num_right = other.num_right;
		
		set_mp(this->left_crit_val, other.left_crit_val);
		set_mp(this->right_crit_val, other.right_crit_val);
	}
	
};


/**
 \brief A simple triangle class, holding three integers referring to points in a vertex set.
 */
class triangle
{
public:
	long long v1; ///< point index 1.
	long long v2; ///< point index 2.
	long long v3; ///< point index 3.
	

	triangle(){};
	
	triangle(long long _v1,long long _v2,long long _v3)
	{
		v1 = _v1;
		v2 = _v2;
		v3 = _v3;
		
	}
	
	
	friend std::ostream & operator<<(std::ostream &os, const triangle & t)
	{
		
		
		os << t.v1 << " " << t.v2 << " " << t.v3;
		
		return os;
	}
	
	
	/**
	 \brief test a triangle for degeneracy.
	 
	 \return true if degenerate, false if not.
	 */
	bool is_degenerate()
	{
		if (v1==v2) {
			return true;
		}
		
		if (v2==v3) {
			return true;
		}
		
		if (v1==v3) {
			return true;
		}
		
		return false;
	}
};



/**
 \brief Bertini_real surface decomposition.
 
 includes methods to add vertices, look up vertices, etc
 */
class surface_decomposition : public decomposition
{
	
	std::vector<face> faces;
	//these counters keep track of the number of things
	
	unsigned int      num_faces; ///< how many faces are there
	
	
	std::vector<int> num_samples_each_face; ///< the number of samples per face
	
	std::vector< std::vector< triangle >> samples; ///< refined triangulation of the surface.
	
	
public:
	
	std::vector< curve_decomposition > mid_slices; ///< the mid slice curves, occurring halfway between critical points of critical curves.
	std::vector< curve_decomposition > crit_slices; ///< critical slice curves, occurring at critial points of critical curves.
	curve_decomposition crit_curve; ///< the main critical curve, being nonsingular.
	curve_decomposition sphere_curve; ///< the intersection of the surface and a sphere.
	
	std::map<std::pair<int,int>,curve_decomposition> singular_curves; ///<  the singular curves, which are formally but not really part of the critical curve.
	int num_singular_curves; ///< how many singular curves are there.
	
	
	/**
	 \brief set up the faces in a surface decomposition from a file
	 
	 \param load_from_me the name of the file to read.
	 */
	void read_faces(boost::filesystem::path load_from_me);
	
	
	/**
	 \brief set up a surface_decomposition from a set of files
	 
	 \param base the folder in which the decomposition lives.
	 */
	void setup(boost::filesystem::path base);
	
	
	/**
	 \brief print a surface decomposition to a set of files located in folder base.
	 
	 \param base the name of the folder to which to print the decomposition.
	 */
	void print(boost::filesystem::path base);
	
	/**
	 \brief print all the faces in a decomposition to a file named outputfile.
	 
	 \param outputfile the name of the file into which to dump the faces (and a bit of header).
	 */
	void print_faces(boost::filesystem::path outputfile);
	
	

	

	/**
	 \brief add a face after it has been constructed.
	 
	 also increments a counter.
	 
	 \param F the constructed face to add to the decomposition.
	 */
	void add_face(const face & F)
	{
		faces.push_back(F);
		num_faces++;
	}
	

	
	

	
	/**
	 \brief The main outer method for decomposing a surface in Bertini_real.
	 
	 \param V the vertex_set into which to deposit the found points, along with their metadata.
	 \param W the input witness set for the surface, complete with points, linears, and patches.
	 \param pi the two linear projection vectors to use.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	void main(vertex_set & V,
						const witness_set & W,
						vec_mp *pi,
						BR_configuration & program_options,
						solver_configuration & solve_options);
	
	
	/**
	\brief Perform a few preliminaries for the surface, including deflation, parsing of the input file, and creation of a randomization method.
	 
	 \param W_surf the witness set for the surface.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	*/
	void beginning_stuff(const witness_set & W_surf,
						BR_configuration & program_options,
						 solver_configuration & solve_options);
	
	
	/**
	 \brief deflate the input system for the surface with respect to computed critical witness points, and split the points according to deflation system.
	 
	 \param split_sets an output variable, contains the witness sets split according to deflation system.
	 \param higher_multiplicity_witness_sets The witness sets computed by the nullspace method for critical conditions, and split according to 'multiplicity', or the number of times they appeared as solutions to the nullspace system.
	 \param program_options The current state of Bertini_real.
	 \param solve_options The current state of the solver.
	 */
	void deflate_and_split(std::map< std::pair<int,int>, witness_set > & split_sets,
						   std::map<int, witness_set > & higher_multiplicity_witness_sets,
						   BR_configuration & program_options,
						   solver_configuration & solve_options);
	
	/**
	 \brief compute the critical points of all singular curves.
	 
	 \see deflate_and_split
	 
	 \param W_singular_crit A return parameter, contains the computed critical points.
	 \param split_sets The split witness sets for the critical curves, produced by deflate_and_split
	 \param V the total vertex set into which to deposit the computed points on the surfaces.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the Solver.
	 */
	void compute_singular_crit(witness_set & W_singular_crit,
							   const std::map<std::pair<int,int>, witness_set> & split_sets,
							   vertex_set & V,
							   BR_configuration & program_options,
							   solver_configuration & solve_options);
	
	
	/**
	 \brief decompose all singular curves.
	 
	 using the computed witness points and critical points, decompose the singular curves.
	 
	 \param W_total_crit the computed critical points of all critical curves, including the critical curve itself, the sphere curve, and the singular curves.
	 \param split_sets The witness sets for the singular curves, produced by nullspace left.
	 \param V the vertex set into which the points are stored.
	 \param program_options The current state of Bertini_real.
	 \param solve_options The current state of the solver.
	 */
	void compute_singular_curves(const witness_set & W_total_crit,
								 const std::map< std::pair<int,int>, witness_set> & split_sets,
								 vertex_set & V,
								 BR_configuration & program_options,
								 solver_configuration & solve_options);
	
	
	
	/**
	 \brief slice the surface at the supplied projection values downstairs.
	 
	 \param W_surf the witness set for the surface itself.
	 \param V the vertex set into which the computed points are committed.
	 \param projection_values_downstairs The pi_0 projection values of the critical points.
	 \param slices the computed slice curve decompositions.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 \param rerun_empty A flag of whether to adjust the tracker config and rerun any empty decompositions.
	 \param kindofslice A string indicating what kind of slice you are decompositon -- purely for screen output.
	 */
	void compute_slices(const witness_set W_surf,
											vertex_set & V,
											vec_mp projection_values_downstairs,
											std::vector< curve_decomposition > & slices,
											BR_configuration & program_options,
											solver_configuration & solve_options,
											bool rerun_empty,
											std::string kindofslice);
	
	
    /**
	 \brief compute a witness set for the critical curve(s), including the critcurve itself, as well as singular curves.
	 
	 \param W_critcurve One of two computed parameters for this functions, this one contains regular witness points for the deflation-free critcurve.
	 \param higher_multiplicity_witness_sets The other computed returned parameter to this function, this is a map between multiplicities and witness sets of that multiplicity.  These must be deflated.
	 \param W_surf The input witness set, for the surface.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
    void compute_critcurve_witness_set(witness_set & W_critcurve,
									   std::map<int, witness_set> & higher_multiplicity_witness_sets,
                                       const witness_set & W_surf,
                                       BR_configuration & program_options,
                                       solver_configuration & solve_options);
    
	
	/**
	 \brief compute critical points of the regular critical curve.
	 
	 \param W_critcurve_crit The computed value, containing the critical points of the critical curve.
	 \param W_surf Witness point for the surface itself.  I believe this should be unnecessary.
	 \param W_critcurve Witness set for the critical curve.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
    void compute_critcurve_critpts(witness_set & W_critcurve_crit, // the computed value
                                   const witness_set & W_surf, // input witness set
                                   witness_set & W_critcurve,
                                   BR_configuration & program_options,
                                   solver_configuration & solve_options);
    
    
	/**
	 \brief Decompose the critical curve with respect to pi_0.
	 
	 Assumes have the critical points and witness points already computed.  This essentially performs interslice.
	 
	 \see curve_decomposition::interslice
	 
	 
	 \param W_critcurve Witness set for the critical curve.
	 \param W_critpts The computed critical point of all critical curves including the singular curves.
	 \param V The vertex set into which to deposit the computed points.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	void compute_critical_curve(const witness_set & W_critcurve,
                                const witness_set & W_critpts,
                                vertex_set & V,
                                BR_configuration & program_options,
                                solver_configuration & solve_options);
	
	
	

	
	
	/**
	 \brief Obtain a witness set for the sphere intersection curve.
	 
	 
	 \param W_surf Witness set for the surface.
	 \param W_intersection_sphere The computed value, returning the witness set for the sphere intersection curve.
	 \param program_options The current state of Bertini_real.
	 \param solve_options The current state of the solver.
	 */
	void compute_sphere_witness_set(const witness_set & W_surf,
									witness_set & W_intersection_sphere,//
									BR_configuration & program_options,//
									solver_configuration & solve_options);
	
	
	/**
	 \brief Compute critical points of the sphere curve, from the previously computed witness points.
	 
	 \see surface_decomposition::compute_sphere_witness_set
	 
	 \param W_intersection_sphere The witness set for the curve, computed by compute_sphere_witness_set
	 \param W_sphere_crit The computed returned value, the critical points of the sphere curve.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	void compute_sphere_crit(const witness_set & W_intersection_sphere,
													 witness_set & W_sphere_crit,
													 BR_configuration & program_options,
													 solver_configuration & solve_options);
	

	/**
	 \brief Performs interslice on the sphere curve, using the critical points and the witness points.
	 
	 \see surface_decomposition::compute_sphere_crit
	 \see surface_decomposition::compute_sphere_witness_set
	 
	 
	 \param W_intersection_sphere Input witness set for the curve.
	 \param W_crit The input critical points.
	 \param V The vertex set into which to commit the computed points.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	void compute_bounding_sphere(const witness_set & W_intersection_sphere,
															 const witness_set & W_crit,
															 vertex_set & V,
															 BR_configuration & program_options,
															solver_configuration & solve_options);
	
	
	/**
	 \brief Main outer method for obtaining faces of the decomposition.
	 
	 Switches between serial_connect and master_connect depending on parallel mode.  Produces the faces of the decomposition by calling the midpoint tracker.
	 
	 \todo Remove pi as an input for this function.
	 
	 \param V The vertex set into which all the points have been collected.
	 \param pi the projections.  Seems unnecessary considering they have already been stored in this object.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	void connect_the_dots(vertex_set & V,
												 vec_mp *pi,
												 BR_configuration & program_options,
												 solver_configuration & solve_options);
	
	
	/**
	 \brief The serial mode for producing the faces. 
	 
	 In a loop, calls make_face.
	 
	 \param V the vertex set into which all the decompositions index.
	 \param md_config The already-set-up midpoint tracker config object.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	void serial_connect(vertex_set & V,
						midpoint_config & md_config,
						solver_configuration & solve_options,
						BR_configuration & program_options);
	
	
	/**
	 \brief The master parallel mode for producing the faces.
	 
	 In a loop, sends workers indices of faces to make, and receives the completed faces.  Then adds the face to this decomposition.
	 
	 \param V the vertex set into which all the decompositions index.
	 \param md_config The already-set-up midpoint tracker config object.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	void master_connect(vertex_set & V,
						midpoint_config & md_config,
						solver_configuration & solve_options,
						BR_configuration & program_options);
	
	/**
	 \brief The worker mode for producing faces.
	 
	 Worker receives the md_config, surface to be connected, and vertex set from the master, then waits for an index to connect up.  It connects, and then sends back the completed face to the master.
	 
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	void worker_connect(solver_configuration & solve_options,
						BR_configuration & program_options);
	
	
	/**
	 \brief Send the next worker the index of the face to make.
	 
	 \param ii which critslice the face will correspond to.
	 \param jj which edge on the critslice the face will correspond to.
	 \param next_worker The MPI ID of the worker to request to make the face.
	 \param mpi_config the current state of MPI.
	 */
	void master_face_requester(int ii, int jj, int next_worker, parallelism_config & mpi_config);
	
	
	/**
	 \brief As a worker, receive the indices of the face to make.
	 
	 \param ii which critslice the face will correspond to.
	 \param jj which edge on the critslice the face will correspond to.
	 \param mpi_config the current state of MPI.
	 */
	void worker_face_requester(int & ii, int & jj, parallelism_config & mpi_config);
	
	
	/**
	 \brief  Actually make the face for the ii-th midslice and the jj-th edge, calling the midtrack solver.
	 
	 \return The constructed face (which may be degenerate).
	 \param ii which critslice the face will correspond to.
	 \param jj which edge on the critslice the face will correspond to.
	 \param V the vertex set containing the vertices for the total decomposition.
	 \param md_config The already-set-up midtrack config.
	 \param solve_options The current state of the solver.
	 \param program_options The current state of Bertini_real.
	 */
	face make_face(int ii, int jj, vertex_set & V,
				   midpoint_config & md_config,
				   solver_configuration & solve_options, BR_configuration & program_options);
	
	
	
	
	
	
	/**
	 \brief Method for sampling a surface so that each face contains approximately the same number of triangles.
	 
	 \todo Parallelize the sampler.
	 
	 \param V The vertex set holding the vertices.
	 \param sampler_options The current state of the sampler program.
	 \param solve_options The current state of the solver.
	 */
	void fixed_sampler(vertex_set &V,
					   sampler_configuration & sampler_options,
					   solver_configuration & solve_options);
	
	
	/**
	 \brief Outer method for connecting the computed samples into triangles.
	 
	 \return A triangulation of the face described by the input rib_indices.
	 \param rib_indices The indices into V for the points.
	 \param V The vertex set holding the points into which the ribs index.
	 \param bin_it_by_projection Indicator of whether to bin by projection or distance.  True is projection, false is distance binning.
	 */
	std::vector< triangle > stitch_triangulation(const std::vector< std::vector< int > > & rib_indices,
												 vertex_set & V,
												 bool bin_it_by_projection);
	
	
	/**
	 \brief Write the results of a sampling run to a folder.
	 
	 \param base_path The name of the folder into which to write the sampling data.
	 */
	void  output_sampling_data(boost::filesystem::path base_path);
	
	
	
	
	
	
	surface_decomposition() : decomposition()
	{
		this->init();
	}
	
	surface_decomposition & operator=(const surface_decomposition& other){
		this->init();
		this->copy(other);
		return *this;
	}
	
	surface_decomposition(const surface_decomposition & other){
		this->init();
		this->copy(other);
	}
	
	
	~surface_decomposition()
	{
		this->clear();
	}
	
	
	
	/**
	 \brief Find a curve with the input name as its input_filename.
	 
	 \return A pointer to the curve with the name, or NULL if it is not found.
	 \param findme The name you want to appear as the input_filename of the curve.
	 */
	curve_decomposition * curve_with_name(const std::string & findme)
	{
		
		if (findme.compare(crit_curve.input_filename().filename().string())==0) {
			return &crit_curve;
		}
		
		if (findme.compare(sphere_curve.input_filename().filename().string())==0) {
			return &sphere_curve;
		}
		for (auto iter = singular_curves.begin(); iter!=singular_curves.end(); ++iter) {
			if (findme.compare(iter->second.input_filename().filename().string())==0) {
				return &(iter->second);
			}
		}
		
		
		for (auto iter = mid_slices.begin(); iter!=mid_slices.end(); ++iter) {
			if (findme.compare(iter->input_filename().filename().string())==0) {
				return &(*iter);
			}
		}
		
		
		for (auto iter = crit_slices.begin(); iter!=crit_slices.end(); ++iter) {
			if (findme.compare(iter->input_filename().filename().string())==0) {
				return &(*iter);
			}
		}
	
		std::cout << "failed to find curve with name " << findme << std::endl;
				
		return NULL;
	}
	
	
	
	/**
	 \brief Single-target MPI send function for a surface_decomposition..
	 
	 \param target The MPI ID of the process to send it to.
	 \param mpi_config The current state of MPI.
	 */
	void send(int target, parallelism_config & mpi_config);
	
	/**
	 \brief Single-source MPI receive function.
	 
	 \param source The MPI ID of the person from whom to receive the surface_decomposition.
	 \param mpi_config The current state of MPI.
	 */
	void receive(int source, parallelism_config & mpi_config);
	
protected:
	
	
	void copy(const surface_decomposition & other)
	{
		decomposition::copy(other);
		
		this->faces = other.faces;
		
		this->num_faces = other.num_faces;
		
		this->mid_slices = other.mid_slices;
		this->crit_slices = other.crit_slices;
		this->crit_curve = other.crit_curve;
		this->singular_curves = other.singular_curves;
		this->num_singular_curves = other.num_singular_curves;
	}
	
	
	void init()
	{
		num_singular_curves = 0;
		num_faces = 0;
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
int find_matching_singular_witness_points(witness_set & W_match,
										  witness_set & W_reject,
										  const witness_set & W,
										  solver_configuration & solve_options);


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
													const witness_set & W);



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
													const witness_set & W);


/**
 \brief Read the file named INfile, and get the S.surf information from it.
 
 
 \see surface_decomposition::print
 
 This function is used when reading a decomposition back in from files.
 
 \param singular_multiplicities A returned value, containing the multiplicity information for the decomposition.
 \param temp_num_mid The number of midslices.
 \param temp_num_crit The number of critslices.  Should be temp_num_mid+1.
 \param INfile The name of the file to read.
 */
void read_summary(std::vector<std::pair<int,int>> & singular_multiplicities,
				  int & temp_num_mid,
				  int & temp_num_crit,
				  boost::filesystem::path INfile);


#endif


