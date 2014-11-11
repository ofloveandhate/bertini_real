function get_indices(br_plotter)





% set up the plotting indices -- which data to plot

% 		first scan for zero columns.
if br_plotter.dimension == 1
		if ~isempty(br_plotter.BRinfo.sampler_data)
			tmpdata = zeros(sum(br_plotter.BRinfo.sampler_data.sample_sizes),br_plotter.BRinfo.num_variables-1);
			counter = 1;
			for ii = 1:br_plotter.BRinfo.num_edges
				for jj = 1:br_plotter.BRinfo.sampler_data.sample_sizes(ii)
					tmpdata(counter,:) = br_plotter.BRinfo.vertices(br_plotter.BRinfo.sampler_data.edge(ii).samples(jj)+1).point;
					counter = counter+1;
				end
			end
		else
			tmpdata = zeros(br_plotter.BRinfo.num_vertices,br_plotter.BRinfo.num_variables-1);
			for ii = 1:br_plotter.BRinfo.num_vertices
				tmpdata(ii,:) = br_plotter.BRinfo.vertices(ii).point(1:br_plotter.BRinfo.num_variables-1);
			end
		end
elseif br_plotter.dimension ==2 
	num_pts = min(1000,br_plotter.BRinfo.num_vertices);
	
	tmpdata = zeros(num_pts,br_plotter.BRinfo.num_variables-1);
	
	for ii = 1:num_pts
		tmpdata(ii,:) = br_plotter.BRinfo.vertices(ii).point(1:br_plotter.BRinfo.num_variables-1);
	end
end



indices_of_nonconst_cols = find_constant_vars(tmpdata);


if length(indices_of_nonconst_cols)>=4
	br_plotter.indices = get_user_indices(indices_of_nonconst_cols,br_plotter.BRinfo);
else
	if br_plotter.dimension == 1
		br_plotter.indices = indices_of_nonconst_cols;
	else
		if length(indices_of_nonconst_cols)==3
			br_plotter.indices = indices_of_nonconst_cols;
		else
			br_plotter.indices = get_user_indices(1:br_plotter.BRinfo.num_variables-1,br_plotter.BRinfo);
		end
	end
end

end





