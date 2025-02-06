%  the main function that Chebfun should call, for interfacing with
%  Bertini_real curves.

% this function assumes that a bertini input file has already been written,
% and that it uses parameter homotopy.  it further assumes that the start
% point has already been written to start_parameters

% all this function does is write the target projection value to the
% `final_parameters` file, call bertini, and read the solution file.

% this function is intended to be bound by use of an anonymous function
% call, leaving target_proj_vals as the only unbound value.

% cf coord assumed to be between 0 and 1.
function solutions = call_bertini_chebfun(cf_coord, pi_out, pi_mid, cycle_num, coordinate_indices, num_vars, point_map, loglevel, logfile)

	solutions = zeros(length(cf_coord),length(coordinate_indices));
	if loglevel>=1
		fprintf(logfile, '\t\t\tcomputing %i additional points\n',length(cf_coord));
	end
	
	for ii = 1:length(cf_coord) % this loop could be replaced by a call to paramotopy's step2
		
		p = cf_coord(ii);
		
		if isKey(point_map,p)
			temp = point_map(p);
			
			if isempty(temp)
				error('empty retrieved point')
			end
			
			solutions(ii,:) = temp(coordinate_indices);
			
			if loglevel>=2
				fprintf(logfile, '\t\tcache lookup successful p=%f\n',p);
			end
		else
			pi = pi_out + (pi_mid - pi_out) * (1-p)^cycle_num;
			
			if loglevel>=2
				fprintf(logfile, '\t\tcomputing point at p=%f\n',p);
				
				if loglevel>=3
					fprintf(logfile, '\t\t\tscaled p=%f\n',pi);
					fprintf(logfile, '\t\t\tout p=%f\n',pi_out);
				end
			end
			
			fid = fopen('final_parameters','w');
			fprintf(fid,'1\n\n');

			fprintf(fid,'%1.17e 0.0\n',pi); 
			fclose(fid);

			bertini('filename','input_BrCf_slice','stifle');  % call from the bertini_tropical code
			s = get_generic_solns('real_finite_solutions',num_vars); 

			if isempty(s)
				error('empty point from real_finite_solutions')
			end

			try
				point_map(p) = real(s);
				solutions(ii,:) = real(s(coordinate_indices))';
			catch
				error('bertini returned no solutions classified as real');
			end
		end
		
		
	end
end