function preprocess_data(br_plotter)
	
	if isfield(br_plotter.BRinfo,'compressed_data')
		br_plotter.data.raw.vertices = br_plotter.BRinfo.compressed_data.fv.vertices;
	else
		disp('matrix-izing data.  consider re-gathering the data to skip this time consuming step.')
		br_plotter.data.raw.vertices = zeros(br_plotter.BRinfo.num_vertices,br_plotter.BRinfo.num_variables-1);
		for ii=1:br_plotter.BRinfo.num_vertices
			br_plotter.data.raw.vertices(ii,:) = transpose(real(br_plotter.BRinfo.vertices(ii).point(1:br_plotter.BRinfo.num_variables-1)));
		end
	
	end
	
	
	if br_plotter.options.use_custom_projection
		try
			try
				proj = br_plotter.options.custom_projection;
				br_plotter.data.space.vertices = proj(br_plotter.data.raw.vertices);
			catch
				disp('unsuccesful attempt to call function naturally on entire point set, rows as points. attempting vectorize')
				proj = br_plotter.options.custom_projection;
				proj = vectorize(proj);
				br_plotter.data.space.vertices = proj(br_plotter.data.raw.vertices);
			end
		catch ME
			disp('unsuccesful attempt to call function naturally or vectorized.  calling one by one.  this can be slow for large data sets')
			fprintf('caught error %s\n',ME.message)

			% preallocate
			test_point = br_plotter.options.custom_projection(transpose(br_plotter.BRinfo.vertices(1).point(1:br_plotter.BRinfo.num_variables-1)));
			
			br_plotter.data.space.vertices = zeros(br_plotter.BRinfo.num_vertices,size(test_point,2));
			for ii = 1:br_plotter.BRinfo.num_vertices
				br_plotter.data.space.vertices(ii,:) = br_plotter.options.custom_projection(transpose(real(br_plotter.BRinfo.vertices(ii).point(1:br_plotter.BRinfo.num_variables-1))));
			end
			br_plotter.BRinfo.num_variables = length(br_plotter.options.custom_projection(br_plotter.BRinfo.vertices(1).point(1:br_plotter.BRinfo.num_variables-1)))+1;
		end
		
		br_plotter.BRinfo.var_names = {};
		for ii = 1:length(br_plotter.BRinfo.vertices(1).point)
			br_plotter.BRinfo.var_names{ii} = ['proj_' num2str(ii)];
		end
	else
		br_plotter.data.space = br_plotter.data.raw;
	end
    
	

	
	get_indices(br_plotter);

	br_plotter.fv.vertices = br_plotter.data.space.vertices(:,br_plotter.indices);
	% faces will be made later, if necessary
	
	if br_plotter.options.use_colorfn
		try
			if br_plotter.options.colorfn_uses_raw
				br_plotter.data.color.vertices = br_plotter.options.colorfn(br_plotter.data.raw.vertices);
			else
				br_plotter.data.color.vertices = br_plotter.options.colorfn(br_plotter.data.space.vertices);
			end
		catch
			if br_plotter.options.colorfn_uses_raw
				out_size = length(br_plotter.options.colorfn(br_plotter.data.raw.vertices(1,:)));
				br_plotter.data.color.vertices = zeros(br_plotter.BRinfo.num_vertices,out_size);
				for ii = 1:br_plotter.BRinfo.num_vertices
					br_plotter.data.color.vertices(ii,:) = br_plotter.options.colorfn(br_plotter.data.raw.vertices(ii,:));
				end
			else
				out_size = length(br_plotter.options.colorfn(br_plotter.data.space.vertices(1,:)));
				br_plotter.data.color.vertices = zeros(br_plotter.BRinfo.num_vertices,out_size);
				for ii = 1:br_plotter.BRinfo.num_vertices
					br_plotter.data.color.vertices(ii,:) = br_plotter.options.colorfn(br_plotter.data.space.vertices(ii,:));
				end
			end
		end
	end

end