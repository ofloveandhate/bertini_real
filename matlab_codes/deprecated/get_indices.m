function [plot_indices] = get_indices(BRinfo)


% set up the plotting indices -- which data to plot

% 		first scan for zero columns.
if ~isempty(BRinfo.sampler_data)
	tmpdata = zeros(sum(BRinfo.sampler_data.sample_sizes),BRinfo.num_variables-1);
	counter = 1;
	for ii = 1:BRinfo.num_edges
		for jj = 1:BRinfo.sampler_data.sample_sizes(ii)
			tmpdata(counter,:) = BRinfo.vertices(BRinfo.sampler_data.edge(ii).samples(jj)+1).point;
			counter = counter+1;
		end
	end
else
	tmpdata = zeros(BRinfo.num_vertices,BRinfo.num_variables-1);
	for ii = 1:BRinfo.num_vertices
		tmpdata(ii,:) = BRinfo.vertices(ii).point(1:BRinfo.num_variables-1);
	end
end

indices_of_nonconst_cols = find_constant_vars(tmpdata);

if length(indices_of_nonconst_cols)<4
	plot_indices = indices_of_nonconst_cols;
else
	plot_indices = get_user_indices(indices_of_nonconst_cols,BRinfo);
end

end




