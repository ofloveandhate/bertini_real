function get_indices(br_plotter)
	if ~isempty(br_plotter.options.user_indices)
		user_set_ind(br_plotter)
	elseif br_plotter.BRinfo.num_variables-1<=3 % assumes the number of variables includes the homvar...
		vanilla_set_ind(br_plotter)
	else
		complicated_set_ind(br_plotter)
	end

end


function user_set_ind(br_plotter)

	% first do error checking
	if length(br_plotter.options.user_indices) > 3
		error('cannot plot an object using %i coordinates, sorry',length(br_plotter.options.user_indices));
	elseif and(length(br_plotter.options.user_indices) < 3,br_plotter.BRinfo.dimension==2)
		error('cannot plot a surface using only %i coordinates, sorry.  use more variables (3, to be exact)',length(br_plotter.options.user_indices));
	elseif and(length(br_plotter.options.user_indices) < 2,br_plotter.BRinfo.dimension==1)
		error('cannot plot a curve using only %i coordinates, sorry.  use more variables (2 or 3, to be exact)',length(br_plotter.options.user_indices));
	elseif any(br_plotter.options.user_indices>br_plotter.BRinfo.num_variables-1)
		error('attempting to plot with indices exceeding the number available');
	end
	
	br_plotter.indices = br_plotter.options.user_indices;
end


function vanilla_set_ind(br_plotter)
	br_plotter.indices = 1:br_plotter.BRinfo.num_variables-1;
end


function complicated_set_ind(br_plotter)

	indices_of_nonconst_cols = get_nonconstind(br_plotter);

	if length(indices_of_nonconst_cols)>=4
		br_plotter.indices = get_user_indices(indices_of_nonconst_cols,br_plotter.BRinfo);
	else
		br_plotter.indices = indices_of_nonconst_cols;
	end

end




function indices_of_nonconst_cols = get_nonconstind(br_plotter)

	num_pts = min(1000,br_plotter.BRinfo.num_vertices);
	
	tmpdata = zeros(num_pts,br_plotter.BRinfo.num_variables-1);
	
	for ii = 1:num_pts
		tmpdata(ii,:) = br_plotter.BRinfo.vertices(ii).point(1:br_plotter.BRinfo.num_variables-1);
	end
	
	indices_of_nonconst_cols = find_constant_vars(tmpdata);
end