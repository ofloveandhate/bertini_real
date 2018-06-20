

function plot_surface_samples(br_plotter)

if and(~isempty(br_plotter.BRinfo.sampler_data),br_plotter.options.render_samples)
	br_plotter.fv.faces = plot_surf_samples_main(br_plotter);
end

end



function sampler_faces = plot_surf_samples_main(br_plotter)

total_num_faces = 0;
if ~isempty(br_plotter.BRinfo.sampler_data)
	for cc = 1:length(br_plotter.options.which_faces)
		ii = br_plotter.options.which_faces(cc);
		total_num_faces = total_num_faces+size(br_plotter.BRinfo.sampler_data{ii},1);
	end
end

sampler_faces = -ones(total_num_faces,3);
default_colors = br_plotter.options.colormap(length(br_plotter.options.which_faces));
ind_so_far = 1;

face_colors = -ones(total_num_faces,3);

if ~isempty(br_plotter.BRinfo.sampler_data)
	for cc = 1:length(br_plotter.options.which_faces)
		ii = br_plotter.options.which_faces(cc);
		
		num_this_time = size(br_plotter.BRinfo.sampler_data{ii},1);
		
		rang = ind_so_far:ind_so_far+num_this_time-1;
		% add the faces to the total face obj
		sampler_faces(rang,:) = br_plotter.BRinfo.sampler_data{ii};
		
		if br_plotter.options.monocolor
			face_colors(rang,:) = repmat(br_plotter.options.monocolor_color,[num_this_time,1]);
		elseif ~br_plotter.options.use_colorfn
			face_colors(rang,:) = repmat(default_colors(ii,:),[num_this_time,1]);
		end
		
		ind_so_far = ind_so_far + num_this_time;
	end
end

% actually do the render
h = patch('faces',sampler_faces,'vertices',br_plotter.fv.vertices);
if br_plotter.options.monocolor
	set(h,...
		'FaceColor',br_plotter.options.monocolor_color,...
		'FaceAlpha',br_plotter.options.sample_alpha,...
		'EdgeColor',br_plotter.options.monocolor_color,...
		'EdgeAlpha',br_plotter.options.sample_triangulation_alpha,...
		'LineWidth',br_plotter.options.edge_width);
elseif br_plotter.options.use_colorfn
	set(h,...
		'FaceVertexCData', br_plotter.data.color.vertices,...
		'FaceColor', 'interp',...
		'FaceAlpha',br_plotter.options.sample_alpha,...
		'EdgeColor','interp',...
		'EdgeAlpha',br_plotter.options.sample_triangulation_alpha,...
		'LineWidth',br_plotter.options.edge_width);
else
	set(h,'FaceVertexCData',face_colors,...
		  'FaceAlpha',br_plotter.options.sample_alpha,...
		  'EdgeAlpha',br_plotter.options.sample_triangulation_alpha,...
		  'FaceColor', 'flat',...
		  'EdgeColor', 'none',...
		  'LineWidth',br_plotter.options.edge_width);
end
br_plotter.handles.faces.samples = h;

end

