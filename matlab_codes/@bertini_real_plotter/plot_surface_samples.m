

function plot_surface_samples(br_plotter)

if ~isempty(br_plotter.BRinfo.sampler_data)
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

sampler_faces = zeros(total_num_faces,3);
by_face_colors = br_plotter.options.colormap(length(br_plotter.BRinfo.sampler_data));
ind_so_far = 1;

tmp_fv = br_plotter.fv;


if ~isempty(br_plotter.BRinfo.sampler_data)
	for cc = 1:length(br_plotter.options.which_faces)
		ii = br_plotter.options.which_faces(cc);
		
		
		
		tmp_fv.faces = br_plotter.BRinfo.sampler_data{ii};
		tmp_fv.faces(any(br_plotter.fv.faces<=0,2),:) = []; % omit problematic faces.
		
		num_this_time = size(tmp_fv.faces,1);
		
		%plot the BR face
		h = patch(tmp_fv);
		
		if br_plotter.options.monocolor
			set(h,...
				'FaceColor',br_plotter.options.monocolor_color,...
				'FaceAlpha',br_plotter.options.face_alpha,...
				'EdgeColor',br_plotter.options.monocolor_color,...
				'EdgeAlpha',br_plotter.options.edge_alpha);
		elseif br_plotter.options.use_colorfn
			set(h,...
				'FaceVertexCData', br_plotter.data.color.vertices,...
				'FaceColor', 'interp',...
				'FaceAlpha',br_plotter.options.face_alpha,'EdgeColor','none','EdgeAlpha',br_plotter.options.edge_alpha);
		else
			set(h,'FaceColor',by_face_colors(ii,:),'FaceAlpha',br_plotter.options.face_alpha,'EdgeColor',0.8*by_face_colors(ii,:),'EdgeAlpha',br_plotter.options.edge_alpha);
		end
	
		br_plotter.handles.surface_samples(end+1) = h;
		
		% add the faces to the total face blabal
		sampler_faces(ind_so_far:ind_so_far+num_this_time-1,:) = tmp_fv.faces;
		ind_so_far = ind_so_far + num_this_time;
	end
end


sampler_faces = sampler_faces(1:ind_so_far-1,:);


end

