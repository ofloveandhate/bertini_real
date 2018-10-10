function br_plotter = plot_edge_points(br_plotter)

curr_axes = br_plotter.axes.main;
text_color = 'k';

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
	end
end
num_non_degen = br_plotter.BRinfo.num_edges - sum(degen);


colors = br_plotter.options.colormap(num_non_degen);


num_nondegen = br_plotter.BRinfo.num_edges - sum(degen);

text_locations = zeros(num_nondegen,length(ind));
text_labels = cell(num_nondegen,1);

%ugh, we now have two counters
nondegen_index = 0;
nondegen_edge_ind = 1;
	
for ii = 1:br_plotter.BRinfo.num_edges
	
	
		
			
	left_plot = br_plotter.data.space.vertices(br_plotter.BRinfo.edges(ii,1),:);
	mid_plot = br_plotter.data.space.vertices(br_plotter.BRinfo.edges(ii,2),:);
	right_plot = br_plotter.data.space.vertices(br_plotter.BRinfo.edges(ii,3),:);
	
	if degen(ii)
		continue;
	else
		nondegen_index = nondegen_index+1;
	end
	
	switch length(ind)
		case {2}
			plotme = [left_plot(ind(1)) left_plot(ind(2));...
				mid_plot(ind(1)) mid_plot(ind(2));...
				right_plot(ind(1)) right_plot(ind(2))];
			plotme = real(plotme);
			
			h = main_plot_function(plotme,[1 2], curr_axis);
			
			text_locations(nondegen_index,:) = [mid_plot(ind(1)) mid_plot(ind(2))];
		case {3}
			
			plotme = [left_plot(ind(1)) left_plot(ind(2)) left_plot(ind(3));...
				mid_plot(ind(1)) mid_plot(ind(2)) mid_plot(ind(3));...
				right_plot(ind(1)) right_plot(ind(2)) right_plot(ind(3))];
			plotme = real(plotme);
			h = main_plot_function(plotme,[1 2 3], curr_axis);
			
			text_locations(nondegen_index,:) = [mid_plot(ind(1)) mid_plot(ind(2)) mid_plot(ind(3))];
			
	end
	
	
	if br_plotter.options.labels
		text_labels{ii} = sprintf('edge %d  ',ii);
	end
	
	br_plotter.handles.edges(nondegen_edge_ind) = h;
	
	
	set(h,'Color',colors(nondegen_edge_ind,:),'LineWidth',3,'LineStyle',style);
	nondegen_edge_ind = nondegen_edge_ind+1;
end


if br_plotter.options.labels
	switch length(ind)
		case {2}
			text_handle = text(text_locations(:,1),text_locations(:,2), text_labels,...
						'HorizontalAlignment','right',...
						'FontSize',br_plotter.options.fontsizes.labels,...
						'Parent',curr_axes,'Color',text_color);
		case{3}
			text_handle = text(text_locations(:,1),text_locations(:,2),text_locations(:,3), text_labels,...
						'HorizontalAlignment','right',...
						'FontSize',br_plotter.options.fontsizes.labels,...
						'Parent',curr_axes,'Color',text_color);
	end
	
	br_plotter.handles.label.edges = text_handle;
end


end
