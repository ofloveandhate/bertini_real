function br_plotter = plot_edge_points(br_plotter)

if isfield(br_plotter.options,'use_fixed_linestyle')
    style = br_plotter.options.linestyle;
else
    if isempty(br_plotter.BRinfo.sampler_data)
        style = '-';
    else
        style = '--';
    end
end


%
curr_axis = br_plotter.axes.main;

ind = br_plotter.indices;


degen = zeros(1,br_plotter.BRinfo.num_edges);
for ii = 1:br_plotter.BRinfo.num_edges
	
	left_plot  = br_plotter.BRinfo.edges(ii,1);
	mid_plot   = br_plotter.BRinfo.edges(ii,2);
	right_plot = br_plotter.BRinfo.edges(ii,3);
	
	if and(left_plot==mid_plot, mid_plot==right_plot)
		degen(ii) = 1;
	else
		num_non_degen = num_non_degen+1;
	end
end
num_non_degen = sum(degen);



colors = br_plotter.options.colormap(num_non_degen);

nondegen_edge_ind = 1;
for ii = 1:br_plotter.BRinfo.num_edges
	
	
	
	left_plot = br_plotter.data.space.vertices(br_plotter.BRinfo.edges(ii,1),:);
	mid_plot = br_plotter.data.space.vertices(br_plotter.BRinfo.edges(ii,2),:);
	right_plot = br_plotter.data.space.vertices(br_plotter.BRinfo.edges(ii,3),:);
	
	if degen(ii)
		continue;
	end
	switch length(ind)
		case {2}
			plotme = [left_plot(ind(1)) left_plot(ind(2));...
				mid_plot(ind(1)) mid_plot(ind(2));...
				right_plot(ind(1)) right_plot(ind(2))];
			plotme = real(plotme);
			
			h = main_plot_function(plotme,[1 2], curr_axis);
			
		case {3}
			
			plotme = [left_plot(ind(1)) left_plot(ind(2)) left_plot(ind(3));...
				mid_plot(ind(1)) mid_plot(ind(2)) mid_plot(ind(3));...
				right_plot(ind(1)) right_plot(ind(2)) right_plot(ind(3))];
			plotme = real(plotme);
			h = main_plot_function(plotme,[1 2 3], curr_axis);
	end
	br_plotter.handles.edges(nondegen_edge_ind) = h;
	set(h,'Color',colors(nondegen_edge_ind,:),'LineWidth',3,'LineStyle',style);
	nondegen_edge_ind = nondegen_edge_ind+1;
end


end
