%curve specific code

function curve_plot(BRinfo,plot_indices)
global plot_params
%

num_non_degen = 0;
for ii = 1:BRinfo.num_edges
	
	left_plot = BRinfo.vertices(BRinfo.edges(ii,1)).point;
	mid_plot = BRinfo.vertices(BRinfo.edges(ii,2)).point;
	right_plot = BRinfo.vertices(BRinfo.edges(ii,3)).point;
	
	if and(left_plot==mid_plot, mid_plot==right_plot)
		continue;
	else
		num_non_degen = num_non_degen+1;
	end
	
end

colors = jet(num_non_degen);








create_axes_br();


label_axes(plot_indices,BRinfo,plot_params.axes.main);

% actually plot
plot_projection(BRinfo,plot_indices);



vertices = plot_vertices(plot_indices, BRinfo);

plot_params.init_cam_pos = adjust_axes_br(vertices,plot_params.axes.main);


plot_params.handles.sample_edges = [];
plot_params.handles.edges = [];

plot_edge_points(plot_indices,BRinfo,0.8*colors);


if ~isempty(BRinfo.sampler_data)
	plot_params.handles.sample_edges = plot_sampler_data(plot_indices, BRinfo.vertices,BRinfo.sampler_data,colors);
end

sphere_plot(BRinfo);


% colorbar
% make_discrete_colorbar(1,num_non_degen,plot_params.fontsize,'edge number');
end













function plot_edge_points(ind, BRinfo,colors)
global plot_params

%
curr_axis = plot_params.axes.main;



nondegen_edge_ind = 1;
for ii = 1:BRinfo.num_edges
	
	left_plot = BRinfo.vertices(BRinfo.edges(ii,1)).point;
	mid_plot = BRinfo.vertices(BRinfo.edges(ii,2)).point;
	right_plot = BRinfo.vertices(BRinfo.edges(ii,3)).point;
	
	if and(left_plot==mid_plot, mid_plot==right_plot)
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
	plot_params.handles.edges(ii) = h;
	set(h,'Color',colors(nondegen_edge_ind,:),'LineWidth',3);
	if isempty(BRinfo.sampler_data)
		set(h,'LineStyle','-');
	else
		set(h,'LineStyle','--');
	end
	nondegen_edge_ind = nondegen_edge_ind+1;
end


end












