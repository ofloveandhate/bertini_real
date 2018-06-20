function br_plotter = plot_faces(br_plotter)

plot_faces_main(br_plotter);

end





function plot_faces_main(br_plotter)

num_faces = length(br_plotter.options.which_faces);
ind = br_plotter.indices;
curr_axis = br_plotter.axes.main;

%
num_total_faces = 0;
for ii = 1:num_faces
	curr_face = br_plotter.BRinfo.faces(br_plotter.options.which_faces(ii));
	num_total_faces = num_total_faces + 2*(curr_face.num_left + curr_face.num_right + 2); %the 2 replaces c(urr_face.top>=0 + curr_face.bottom>=0)
end
num_total_faces = num_total_faces*2;
br_plotter.fv.faces = nan(num_total_faces, 3);



txt = cell(br_plotter.BRinfo.num_faces,1);
pos = zeros(br_plotter.BRinfo.num_faces,length(ind));


face_colors = nan(num_total_faces,3); %used for the default colors
default_face_colors = br_plotter.options.colormap(num_faces);


total_face_index = 1;
for cc = 1:length(br_plotter.options.which_faces)
	ii = br_plotter.options.which_faces(cc);
	if br_plotter.BRinfo.faces(ii).midslice_index == -1
		continue
	end
	
	if br_plotter.options.labels
		txt{ii} = ['\newline' num2str(ii)];
		pos(ii,:) = br_plotter.data.space.vertices(br_plotter.BRinfo.faces(ii).midpoint,ind);
	end
	
	
	pass = 1; left_edge_counter = 1;  right_edge_counter = 1;
	
	while 1
		
		switch pass
			case 1  %the top edge
				pass = pass+1;
				
				
				if br_plotter.BRinfo.faces(ii).top<0
					continue;
				end
				
				curr_edge = -10;
				if strcmp(br_plotter.BRinfo.faces(ii).system_top,'input_critical_curve')
					curr_edge = br_plotter.BRinfo.crit_curve.edges(br_plotter.BRinfo.faces(ii).top,:);
				elseif strcmp(br_plotter.BRinfo.faces(ii).system_top,'input_surf_sphere')
					curr_edge = br_plotter.BRinfo.sphere_curve.edges(br_plotter.BRinfo.faces(ii).top,:);
				else
					%do a lookup (slower)
					for zz = 1:length(br_plotter.BRinfo.singular_curves)
						if strcmp(br_plotter.BRinfo.singular_names{zz},br_plotter.BRinfo.faces(ii).system_top)
							curr_edge = br_plotter.BRinfo.singular_curves{zz}.edges(br_plotter.BRinfo.faces(ii).top,:);
						end
					end
					
				end
				
				if curr_edge<0
					continue;
				end
				
				curr_edge = curr_edge([3 2 1]);
				
			case 2  %the bottom edge
				pass = pass+1;
				
				if br_plotter.BRinfo.faces(ii).bottom<0
					continue;
				end
				
				curr_edge = -10;
				if strcmp(br_plotter.BRinfo.faces(ii).system_bottom,'input_critical_curve')
					curr_edge = br_plotter.BRinfo.crit_curve.edges(br_plotter.BRinfo.faces(ii).bottom,:);
				elseif strcmp(br_plotter.BRinfo.faces(ii).system_bottom,'input_surf_sphere')
					curr_edge = br_plotter.BRinfo.sphere_curve.edges(br_plotter.BRinfo.faces(ii).bottom,:);
				else
					%do a lookup
					for zz = 1:length(br_plotter.BRinfo.singular_curves)
						if strcmp(br_plotter.BRinfo.singular_names{zz},br_plotter.BRinfo.faces(ii).system_bottom)
							curr_edge = br_plotter.BRinfo.singular_curves{zz}.edges(br_plotter.BRinfo.faces(ii).bottom,:);
						end
					end
					
				end
				
				
				if curr_edge<0
					continue;
				end
				
			case 3  %the left edges
				if left_edge_counter <= br_plotter.BRinfo.faces(ii).num_left
					if br_plotter.BRinfo.faces(ii).left(left_edge_counter)<0 %an error check
						continue;
					end
					
					slice_ind = br_plotter.BRinfo.faces(ii).midslice_index; 
					edge_ind = br_plotter.BRinfo.faces(ii).left(left_edge_counter); 
					
					curr_edge = br_plotter.BRinfo.critpoint_slices{slice_ind}.edges(edge_ind,:);
					left_edge_counter = left_edge_counter +1; %increment
					
					
				else
					pass = pass+1;
					continue;
				end
			case 4 %the right edges
				if right_edge_counter <= br_plotter.BRinfo.faces(ii).num_right
					
					if br_plotter.BRinfo.faces(ii).right(right_edge_counter)<0
						continue;
					end
					
					slice_ind = br_plotter.BRinfo.faces(ii).midslice_index+1;
					edge_ind = br_plotter.BRinfo.faces(ii).right(right_edge_counter);
					curr_edge = br_plotter.BRinfo.critpoint_slices{slice_ind}.edges(edge_ind,:);
					right_edge_counter = right_edge_counter +1;
					
					curr_edge = curr_edge([3 2 1]);
					
				else
					pass = pass+1;
					continue;
				end
			case 5
				break;
				
		end
		
		t1 = [curr_edge(1) curr_edge(2) br_plotter.BRinfo.faces(ii).midpoint];
		t2 = [curr_edge(2) curr_edge(3) br_plotter.BRinfo.faces(ii).midpoint];
		
		br_plotter.fv.faces(total_face_index,:) = t1;
		br_plotter.fv.faces(total_face_index+1,:) = t2;
		
		if br_plotter.options.monocolor
			face_colors(total_face_index:total_face_index+1,:) = repmat(br_plotter.options.monocolor_color,[2,1]);
		elseif ~br_plotter.options.use_colorfn
			face_colors(total_face_index:total_face_index+1,:) = repmat(default_face_colors(ii,:),[2,1]);
		end
		
		total_face_index = total_face_index+2;
	end
	
end

if br_plotter.options.monocolor
	h = patch(br_plotter.fv,...
		'FaceColor',br_plotter.options.monocolor_color,...
		'FaceAlpha',br_plotter.options.face_alpha,...
		'EdgeColor',br_plotter.options.monocolor_color,...
		'EdgeAlpha',br_plotter.options.raw_triangulation_alpha,...
		'LineWidth',br_plotter.options.edge_width,...
		'Parent',curr_axis);
elseif br_plotter.options.use_colorfn
	h = patch(br_plotter.fv,...
		'FaceVertexCData', br_plotter.data.color.vertices, ...
		'FaceColor', 'interp',...
		'FaceAlpha',br_plotter.options.face_alpha,...
		'EdgeColor','interp',...
		'EdgeAlpha',br_plotter.options.raw_triangulation_alpha,...
		'LineWidth',br_plotter.options.edge_width,...
		'Parent',curr_axis);
else
	h = patch(br_plotter.fv,...
		'FaceVertexCData',face_colors,...
		'FaceColor','flat',...
		'FaceAlpha',br_plotter.options.face_alpha,...% 		
		'EdgeAlpha',br_plotter.options.raw_triangulation_alpha,...
		'LineWidth',br_plotter.options.edge_width,...
		'Parent',curr_axis);
end
br_plotter.handles.faces.raw = h;



if br_plotter.options.labels
	switch length(ind)
		case 2
			br_plotter.handles.face_labels = text(pos(:,1),pos(:,2),txt,'Parent',curr_axis,'HorizontalAlignment','center','VerticalAlignment','top');

		case 3
			br_plotter.handles.face_labels = text(pos(:,1),pos(:,2),pos(:,3),txt,'Parent',curr_axis,'HorizontalAlignment','center','VerticalAlignment','top');

		otherwise
				error('length of ind is not 2 or 3...')
	end

	set(br_plotter.handles.face_labels,'visible','off');
end
end








