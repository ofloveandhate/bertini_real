function br_plotter = plot_vertices(br_plotter)

markersize = br_plotter.options.markersize;
ind = br_plotter.indices;
curr_axis = br_plotter.axes.main;
fontsize = br_plotter.options.fontsizes.labels;


names = {'UNSET', 'CRITICAL', 'SEMICRITICAL', ...
			'MIDPOINT', 'ISOLATED', 'NEW', ...
			'CURVE_SAMPLE_POINT', 'SURFACE_SAMPLE_POINT', 'REMOVED', ...
			'PROBLEMATIC'};
		
markers = {'x','o','s','d','^','v','>','<','p','h'};
colors = [1 0 0;lines(length(names)-1)];

for ii = 1:length(names)
	indexed_markers.(names{ii}).marker = markers{ii};
	indexed_markers.(names{ii}).color = colors(ii,:);
end


catnames = cell(length(names),1);

curr_index = 100; %sadly set manually...  this corresponds to header files in bertini_real.  
  % see data_type.hpp
  
  
  
for ii = 1:length(names)
	catnames{ii} = ['ind_' num2str(curr_index)];
	available_types.(catnames{ii}) = names{ii};
	curr_index = curr_index+1;
end



curr_types = zeros(br_plotter.BRinfo.num_vertices, 1);

if br_plotter.options.labels
	labels = cell(br_plotter.BRinfo.num_vertices,1);
end

br_plotter.fv.vertices = zeros(br_plotter.BRinfo.num_vertices,length(ind));


types = 100:109; %initialize to blank  eww i hate these magic constants.  sorry.

for ii=1:br_plotter.BRinfo.num_vertices
	
    br_plotter.fv.vertices(ii,:) = transpose(real(br_plotter.BRinfo.vertices(ii).point(ind)));
	
	curr_types(ii) = br_plotter.BRinfo.vertices(ii).type;
	
	
% 	if isempty(find(unique(types)==curr_types(ii),1))
% 		types = [types curr_types(ii)];
% 	end
	if br_plotter.options.labels
		labels{ii} = [ '    ' num2str(ii-1)];
	end
end



if br_plotter.options.render_vertices
	rendered_counter = 0;
	for ii = 1:length(types)

		curr_name = names{types(ii)-99};

		curr_indices_logical = curr_types==types(ii);

		if sum(curr_indices_logical)==0
			continue;
		else
			rendered_counter = rendered_counter+1;
		end

		curr_points = br_plotter.fv.vertices(curr_indices_logical,:);


		h = main_plot_function(curr_points,   1:length(ind), curr_axis);


		set(h,'LineStyle','none');
		set(h,'Marker', indexed_markers.(curr_name).marker,'MarkerSize',markersize,'MarkerFaceColor',0.6*indexed_markers.(curr_name).color,'MarkerEdgeColor',indexed_markers.(curr_name).color);
		set(h,'Color',colors(ii,:));

		local_catname = sprintf('ind_%d',types(ii));



		br_plotter.handles.vertices.(available_types.(local_catname)) = h;



		br_plotter.legend.vertices.handles(rendered_counter) = h;
		br_plotter.legend.vertices.types{rendered_counter} = available_types.(local_catname);
		br_plotter.legend.vertices.text{rendered_counter} = ['point type ' available_types.(local_catname)];


		if br_plotter.options.labels
			switch length(ind)

				case 2
					br_plotter.handles.vertex_text.(available_types.(local_catname)) = text(curr_points(:,1), curr_points(:,2), labels(curr_types==types(ii)),...
						'HorizontalAlignment','left','VerticalAlignment','bottom',...
						'FontSize',fontsize,'Parent',curr_axis,'Color',colors(ii,:));
				case 3
					br_plotter.handles.vertex_text.(available_types.(local_catname)) = text(curr_points(:,1), curr_points(:,2), curr_points(:,3), labels(curr_types==types(ii)),...
						'HorizontalAlignment','left','VerticalAlignment','bottom',...
						'FontSize',fontsize,'Parent',curr_axis,'Color',colors(ii,:));
				otherwise
			end
		end

	end
end

end% re function
