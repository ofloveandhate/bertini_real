function vertices = plot_vertices(ind, BRinfo)
global plot_params

markersize = 13;

curr_axis = plot_params.axes.main;

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

curr_index = 100; %sadly set manually...  this corresponds to header files in bertini_real
for ii = 1:length(names)

	catnames{ii} = ['ind_' num2str(curr_index)];
	available_types.(catnames{ii}) = names{ii};
	curr_index = curr_index+1;
end

plotme = zeros(BRinfo.num_vertices, length(ind));



indices = zeros(BRinfo.num_vertices, 1);

labels = cell(1,BRinfo.num_vertices);


vertices = zeros(BRinfo.num_vertices,length(ind));


types = []; %initialize to blank
for ii=1:BRinfo.num_vertices
    vertices(ii,:) = transpose(BRinfo.vertices(ii).point(ind));
	plotme(ii,:) = transpose(BRinfo.vertices(ii).point(ind)); % is this unnecessary?
	indices(ii) = BRinfo.vertices(ii).type;
	if isempty(find(unique(types)==indices(ii),1))
		types = [types indices(ii)];
	end
	
	labels{ii} = [ '    ' num2str(ii-1)];
end


for ii = 1:length(types)
	
	curr_name = names{types(ii)-99};
	
	curr_points = real(plotme(indices==types(ii),:));
	h = main_plot_function(curr_points,1:length(ind), curr_axis);
	
	
	set(h,'LineStyle','none');
	set(h,'Marker',indexed_markers.(curr_name).marker,'MarkerSize',markersize,'MarkerFaceColor',0.6*indexed_markers.(curr_name).color,'MarkerEdgeColor',indexed_markers.(curr_name).color);
	set(h,'Color',colors(ii,:));
	local_catname = ['ind_' num2str(types(ii))];
	
	
	plot_params.legend.vertices.handles(ii) = h;
	plot_params.handles.vertices.(available_types.(local_catname)) = h;
	
	plot_params.legend.vertices.types{ii} = available_types.(local_catname);
	
	plot_params.legend.vertices.text{ii} = ['point type ' available_types.(local_catname)];
	

	
	switch length(ind)

		case 2
			plot_params.handles.vertex_text.(available_types.(local_catname)) = text(curr_points(:,1), curr_points(:,2), labels(indices==types(ii)),...
				'HorizontalAlignment','left','VerticalAlignment','bottom',...
				'FontSize',plot_params.fontsize-4,'Parent',curr_axis,'Color',colors(ii,:));
		case 3
			plot_params.handles.vertex_text.(available_types.(local_catname)) = text(curr_points(:,1), curr_points(:,2), curr_points(:,3), labels(indices==types(ii)),...
				'HorizontalAlignment','left','VerticalAlignment','bottom',...
				'FontSize',plot_params.fontsize-4,'Parent',curr_axis,'Color',colors(ii,:));
		otherwise
			
			


end


vertices = real(vertices);
end
