function handles = plot_curve_samples(br_plotter,curve,style,desiredcolor)

if nargin == 4
	repeat_colors = true;
else
	repeat_colors = false;
end

ind = br_plotter.indices;

sampler_data = curve.sampler_data;

touches = ones(curve.num_edges,1);
if br_plotter.options.touching_edges_only
	num_touching = 0;
	touching_edge_indices = zeros(curve.num_edges,1); %preallocate
	for ii = 1:curve.num_edges
		if edge_touches_faces(br_plotter, ii, curve)
			num_touching = num_touching+1;
			touching_edge_indices(num_touching) = ii;
		else
			touches(ii) = 0;
		end
	end
	touching_edge_indices = touching_edge_indices(1:num_touching); %trim the fat
else
	num_touching = curve.num_edges;
	touching_edge_indices = 1:num_touching;
end




num_non_degen = 0;

nondegen = zeros(curve.num_edges,1);
for ii = 1:size(sampler_data.sample_sizes,1)
	
	left_point = br_plotter.BRinfo.vertices(sampler_data.edge(ii).samples(1)).point(1:br_plotter.BRinfo.num_variables-1);
	right_point = br_plotter.BRinfo.vertices(sampler_data.edge(ii).samples(end)).point(1:br_plotter.BRinfo.num_variables-1);
	
	if norm(left_point-right_point)<1e-8
		continue;
	else
		nondegen(ii) = 1;
		num_non_degen = num_non_degen+1;
	end
	
end

indicator = and(touches, nondegen);
num_to_plot = sum(indicator);

if repeat_colors
	colors = repmat(desiredcolor,[num_to_plot 1]);
else
	colors = br_plotter.options.colormap(num_to_plot);
end


curr_axis = br_plotter.axes.main;
counter = 1;

handles = zeros(num_to_plot,1);

for ii = 1:length(sampler_data.edge)
	
	if indicator(ii)==0
       continue; 
	end

    curr_samples = sampler_data.edge(ii).samples;

	plotme = br_plotter.fv.vertices(curr_samples,:);
	

	h = main_plot_function(plotme,1:length(ind),curr_axis);

	
	set(h,'Color',colors(counter,:));
	set(h,'LineWidth',br_plotter.options.line_thickness);
    set(h,'LineStyle',style);
	
	
    handles(counter) = h;
	
	
	counter = counter+1;
	
end




end













