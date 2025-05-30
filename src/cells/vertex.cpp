#include "cells/vertex.hpp"

void Vertex::set_point(const vec_mp new_point)
{
	change_prec_vec_mp(this->pt_mp_, new_point->curr_prec);
	vec_cp_mp(this->pt_mp_, new_point);
}




void Vertex::send(int target, ParallelismConfig & mpi_config) const
{
	// std::cout << "send vertex to " << target << std::endl;

	send_vec_mp(pt_mp_, target);

	send_vec_mp(projection_values_, target);

	int * buffer = new int[4];
	buffer[0] = type_;
	buffer[1] = input_filename_index_;
	buffer[2] = input_filename_indices_.size();
	buffer[3] = path_numbers_ending_here_.size();

	MPI_Send(buffer, 4, MPI_INT, target, VERTEX, mpi_config.comm());

	MPI_Send(&input_filename_indices_.front(),     input_filename_indices_.size(), MPI_INT, target, VERTEX, mpi_config.comm());
	
	MPI_Send(&path_numbers_ending_here_.front(), path_numbers_ending_here_.size(), MPI_INT, target, VERTEX, mpi_config.comm());



	delete[] buffer;

}


void Vertex::receive(int source, ParallelismConfig & mpi_config)
{
	// std::cout << "recv vertex from " << source << std::endl;

	MPI_Status statty_mc_gatty;

	receive_vec_mp(pt_mp_, source);
	receive_vec_mp(projection_values_, source);


	int * buffer = new int[4];

	MPI_Recv(buffer, 4, MPI_INT, source, VERTEX, mpi_config.comm(), &statty_mc_gatty);

	type_ = static_cast<VertexType>(buffer[0]);
	input_filename_index_ = buffer[1];

	int num_filename_indices = buffer[2];
	int num_paths = buffer[3];

	input_filename_indices_.resize(num_filename_indices);
	MPI_Recv(&input_filename_indices_.front(), num_filename_indices, MPI_INT, source, VERTEX, mpi_config.comm(), &statty_mc_gatty);
		
	path_numbers_ending_here_.resize(num_paths);
	MPI_Recv(&path_numbers_ending_here_.front(), num_paths,          MPI_INT, source, VERTEX, mpi_config.comm(), &statty_mc_gatty);


	delete[] buffer;
}
