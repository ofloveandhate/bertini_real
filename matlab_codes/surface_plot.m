%  surface specific code



function [fv,sampler_faces] = surface_plot(BRinfo,ind)
global plot_params

% material shiny


create_axes_br();

label_axes(ind,BRinfo,plot_params.axes.main);

sphere_plot(BRinfo);

fv.vertices = plot_vertices(ind, BRinfo);

plot_params.init_cam_pos = adjust_axes_br(fv.vertices,plot_params.axes.main);


plot_surface_edges(BRinfo,ind);


fv.faces = plot_faces(BRinfo, ind, fv);


plot_projection(BRinfo,ind);

sampler_faces = plot_surface_samples(BRinfo,fv);

if ~isempty(sampler_faces)
	fv.faces = sampler_faces;
end

sync_axes();


end






function plot_surface_edges(BRinfo,plot_indices)
global plot_params
%
line_thickness = 3;

crit_curve = zeros(3, BRinfo.num_variables-1, BRinfo.crit_curve.num_edges);

handle_counter = 0;
plot_params.handles.critcurve_labels = [];
plot_params.handles.spherecurve_labels = [];
plot_params.handles.refinements.critcurve = [];
plot_params.handles.critcurve = [];
critcurve_counter = 0;

if BRinfo.crit_curve.num_edges>0
	curr_axes = plot_params.axes.main;
	colors = 0.8*jet(BRinfo.crit_curve.num_edges);
	
	
	
	
	for ii =1:BRinfo.crit_curve.num_edges
		
		for jj = 1:3
			crit_curve(jj,:,ii) = real(transpose(BRinfo.vertices(BRinfo.crit_curve.edges(ii,jj)).point(1:BRinfo.num_variables-1)));
		end
		
		h = plot3(crit_curve(:,1,ii),crit_curve(:,2,ii),crit_curve(:,3,ii),'Parent',curr_axes);
		set(h,'Color',colors(ii,:));
		set(h,'LineStyle','-','LineWidth',line_thickness);
		
		critcurve_counter = critcurve_counter+1;
		plot_params.handles.critcurve(critcurve_counter) = h;
		
		
		if ii==1
			handle_counter = handle_counter+1;
			plot_params.legend.surface_edges.handles(handle_counter) = h;
			plot_params.legend.surface_edges.text{handle_counter} = 'critical curve';
		end
		
		new_handles = text(crit_curve(2,1,ii),crit_curve(2,2,ii),crit_curve(2,3,ii),['critcurve ' num2str(ii-1) '  '],'HorizontalAlignment','right','FontSize',plot_params.fontsize-4,...
			'Parent',curr_axes,'Color','r');
		plot_params.handles.critcurve_labels = [plot_params.handles.critcurve_labels; new_handles];
		
		
	end
	
	
	if ~isempty(BRinfo.crit_curve.sampler_data)
		plot_params.handles.refinements.critcurve = plot_sampler_data(plot_indices, BRinfo.vertices,BRinfo.crit_curve.sampler_data,colors);
	end
	
	
end

plot_params.handles.refinements.spherecurve = [];
plot_params.handles.spherecurve = [];
spherecurve_counter = 0;
if isfield(BRinfo,'sphere_curve')
	if BRinfo.sphere_curve.num_edges>0
		
		curr_axes = plot_params.axes.main;
		sphere_curve = zeros(3, BRinfo.num_variables-1, BRinfo.sphere_curve.num_edges);
		
		colors = 0.8*jet(BRinfo.sphere_curve.num_edges);
		
		for ii =1:BRinfo.sphere_curve.num_edges
			
			for jj = 1:3
				sphere_curve(jj,:,ii) = real(transpose(BRinfo.vertices(BRinfo.sphere_curve.edges(ii,jj)).point(1:BRinfo.num_variables-1)));
			end
			
			h = plot3(sphere_curve(:,1,ii),sphere_curve(:,2,ii),sphere_curve(:,3,ii),'Parent',curr_axes);
			set(h,'Color',colors(ii,:));
			set(h,'LineStyle','-.','LineWidth',line_thickness);
			
			
			spherecurve_counter = spherecurve_counter+1;
			plot_params.handles.spherecurve(spherecurve_counter) = h;
			
			if ii==1
				handle_counter = handle_counter+1;
				plot_params.legend.surface_edges.handles(handle_counter) = h;
				plot_params.legend.surface_edges.text{handle_counter} = 'sphere curve';
			end
			
			new_handles = text(sphere_curve(2,1,ii),sphere_curve(2,2,ii),sphere_curve(2,3,ii),['spherecurve ' num2str(ii-1) '  '],'HorizontalAlignment','right','FontSize',plot_params.fontsize-4,...
				'Parent',curr_axes,'Color','c');
			plot_params.handles.spherecurve_labels = [plot_params.handles.spherecurve_labels; new_handles];
			
			
		end
		
		if ~isempty(BRinfo.sphere_curve.sampler_data)
			plot_params.handles.refinements.spherecurve = plot_sampler_data(plot_indices, BRinfo.vertices,BRinfo.sphere_curve.sampler_data,colors);
		end
		
		
	end
end




colors = 0.7*jet(length(BRinfo.midpoint_slices));

plot_params.handles.midtext = [];

plot_params.handles.midslices = [];
plot_params.handles.refinements.midslice = [];
midslice_counter = 0;
curr_axes = plot_params.axes.main;
for kk = 1:length(BRinfo.midpoint_slices)
% 	BRinfo
% 	BRinfo.midpoint_slices{kk}.num_edges
	midslice = zeros(3, BRinfo.num_variables-1, BRinfo.midpoint_slices{kk}.num_edges);
	
	text_positions = zeros(3,BRinfo.midpoint_slices{kk}.num_edges);
	textme = cell(1,BRinfo.midpoint_slices{kk}.num_edges);
	
	
	for ii =1:BRinfo.midpoint_slices{kk}.num_edges
		for jj = 1:3
			midslice(jj,:,ii) = real(transpose(BRinfo.vertices(BRinfo.midpoint_slices{kk}.edges(ii,jj)).point(1:BRinfo.num_variables-1)));
		end
		
		h = plot3(midslice(:,1,ii),midslice(:,2,ii),midslice(:,3,ii),...
			'Parent',curr_axes);
		set(h,'Color',colors(kk,:));
		set(h,'LineStyle',':','LineWidth',line_thickness);
		
		text_positions(:,ii) = midslice(2,1:3,ii);
		textme{ii} = ['mid ' num2str(kk-1) '.' num2str(ii-1) '  '];
		
		midslice_counter = midslice_counter + 1;
		plot_params.handles.midslices(midslice_counter) = h;
		
		if ~isempty(BRinfo.midpoint_slices{kk}.sampler_data)
			refinement_colors = jet(length(BRinfo.midpoint_slices{kk}.edges));
			plot_params.handles.refinements.midslice = horzcat(plot_params.handles.refinements.midslice,plot_sampler_data(plot_indices, BRinfo.vertices,BRinfo.midpoint_slices{kk}.sampler_data,refinement_colors));
		end
		
	end
	
	if kk==1
		handle_counter = handle_counter+1;
		plot_params.legend.surface_edges.handles(handle_counter) = h;
		plot_params.legend.surface_edges.text{handle_counter} = 'midslice';
	end
	
	
	plot_params.handles.midtext = [plot_params.handles.midtext;...
		text(text_positions(1,:),text_positions(2,:),text_positions(3,:), textme,'HorizontalAlignment','right','FontSize',plot_params.fontsize-4,...
		'Parent',curr_axes,'Color','g')];
	
	
end






plot_params.handles.crittext = [];
colors = 0.7*jet(length(BRinfo.critpoint_slices));
plot_params.handles.critslices = [];
critslice_counter = 0;
curr_axes = plot_params.axes.main;
plot_params.handles.refinements.critslice = [];
for kk = 1:length(BRinfo.critpoint_slices)
	critslice = zeros(3, BRinfo.num_variables-1, BRinfo.critpoint_slices{kk}.num_edges);
	
	text_positions = zeros(3,BRinfo.critpoint_slices{kk}.num_edges);
	textme = cell(1,BRinfo.critpoint_slices{kk}.num_edges);
	for ii =1:BRinfo.critpoint_slices{kk}.num_edges
		for jj = 1:3
			critslice(jj,:,ii) = real(transpose(BRinfo.vertices(BRinfo.critpoint_slices{kk}.edges(ii,jj)).point(1:BRinfo.num_variables-1)));
		end
		
		h = plot3(critslice(:,1,ii),critslice(:,2,ii),critslice(:,3,ii),'Parent',curr_axes);
		set(h,'Color',colors(kk,:));
		set(h,'LineStyle','--','LineWidth',line_thickness);
		
		text_positions(:,ii) = critslice(2,1:3,ii);
		textme{ii} = ['crit ' num2str(kk-1) '.' num2str(ii-1) '  '];
		
		critslice_counter = critslice_counter + 1;
		plot_params.handles.critslices(critslice_counter) = h;
		
		
		if ~isempty(BRinfo.critpoint_slices{kk}.sampler_data)
			a = plot_sampler_data(plot_indices, BRinfo.vertices,BRinfo.critpoint_slices{kk}.sampler_data,colors);
			plot_params.handles.refinements.critslice = horzcat(plot_params.handles.refinements.critslice,a);
		end
		
		
	end
	if kk==1
		handle_counter = handle_counter+1;
		plot_params.legend.surface_edges.handles(handle_counter) = h;
		plot_params.legend.surface_edges.text{handle_counter} = 'critslice';
	end
	
	plot_params.handles.crittext = [plot_params.handles.crittext;...
		text(text_positions(1,:),text_positions(2,:),text_positions(3,:), textme,'HorizontalAlignment','right',...
		'FontSize',plot_params.fontsize-4,'Parent',curr_axes,'Color','b')];
	
	
	
end



plot_params.handles.singtext = [];

plot_params.handles.singular_curves = [];

plot_params.handles.refinements.singularcurve = [];

if isfield(BRinfo,'singular_curves')
	colors = 0.7*jet(length(BRinfo.singular_curves));
	singular_counter = 0;
	for kk = 1:length(BRinfo.singular_curves)
		midslice = zeros(3, BRinfo.num_variables-1, BRinfo.singular_curves(kk).num_edges);

		text_positions = zeros(3,BRinfo.singular_curves(kk).num_edges);
		textme = cell(1,BRinfo.singular_curves(kk).num_edges);


		for ii =1:BRinfo.singular_curves(kk).num_edges
			for jj = 1:3
				midslice(jj,:,ii) = real(transpose(BRinfo.vertices(BRinfo.singular_curves(kk).edges(ii,jj)).point(1:BRinfo.num_variables-1)));
			end

			h = plot3(midslice(:,1,ii),midslice(:,2,ii),midslice(:,3,ii),...
				'Parent',curr_axes);
			set(h,'Color','r');
			set(h,'LineStyle','-','LineWidth',line_thickness+2);

			text_positions(:,ii) = midslice(2,1:3,ii);
			textme{ii} = ['sing ' num2str(kk-1) '.' num2str(ii-1) '  '];

			singular_counter = singular_counter + 1;
			plot_params.handles.singular_curves(singular_counter) = h;

			if ~isempty(BRinfo.singular_curves(kk).sampler_data)
				refinement_colors = jet(length(BRinfo.singular_curves(kk).edges));
				plot_params.handles.refinements.singularcurve = horzcat(plot_params.handles.refinements.singularcurve,plot_sampler_data(plot_indices, BRinfo.vertices,BRinfo.singular_curves(kk).sampler_data,refinement_colors));
			end

		end

		if kk==1
			handle_counter = handle_counter+1;
			plot_params.legend.surface_edges.handles(handle_counter) = h;
			plot_params.legend.surface_edges.text{handle_counter} = 'singular';
		end


		plot_params.handles.singtext = [plot_params.handles.singtext;...
			text(text_positions(1,:),text_positions(2,:),text_positions(3,:), textme,'HorizontalAlignment','right','FontSize',plot_params.fontsize-4,...
			'Parent',curr_axes,'Color','g')];


	end
end




end




function [sampler_faces] = plot_surface_samples(BRinfo,fv)
global plot_params



colors = jet(length(BRinfo.sampler_data));
sampler_faces = [];
plot_params.handles.surface_samples = [];
if ~isempty(BRinfo.sampler_data)
	for ii = 1:length(BRinfo.sampler_data)
		
		
		fv.faces = BRinfo.sampler_data{ii}+1;
		fv.faces(any(fv.faces<=0,2),:) = []; % omit problematic faces.
		sampler_faces = [sampler_faces;fv.faces];
		h = patch(fv);
		
		set(h,'FaceColor',colors(ii,:),'FaceAlpha',0.5,'EdgeColor',0.985*colors(ii,:),'EdgeAlpha',0.5);%,'EdgeColor',0.985*colors(ii,:),'EdgeAlpha',0.5
		plot_params.handles.surface_samples(ii) = h;
	end
end


end


















function stl_faces = plot_faces(BRinfo, ind, fv)
global plot_params
%
num_total_faces = 0;
for ii = 1:BRinfo.num_faces
	curr_face = BRinfo.faces(ii);
	num_total_faces = num_total_faces + curr_face.num_left + curr_face.num_right + curr_face.top>=0 + curr_face.bottom>=0;
end
num_total_faces = num_total_faces*2;
stl_faces = zeros(num_total_faces, 3);
curr_face_index = 1;

curr_axis = plot_params.axes.main;

txt = cell(BRinfo.num_faces,1);
pos = zeros(BRinfo.num_faces,length(ind));
plot_params.handles.faces = [];


colors = jet(BRinfo.num_faces);


local.vertices = fv.vertices;

for ii = 1:BRinfo.num_faces
	if BRinfo.faces(ii).midslice_index == -1
		continue
	end
	
	num_triangles= 2*(2 + BRinfo.faces(ii).num_left + BRinfo.faces(ii).num_right);
	
	local.faces = zeros(num_triangles,3);
	
	
	
% 	% set the midpoint of the face for all triangles to be the first row
	pt = transpose(BRinfo.vertices(BRinfo.faces(ii).midpoint+1).point(ind));

	
	txt{ii} = ['\newline' num2str(ii-1)];
	pos(ii,:) = pt(1:length(ind));
	
	

	pass = 1; left_edge_counter = 1;  right_edge_counter = 1;
	
	local_face_index = 1;
	
	while 1
		
		switch pass
			case 1  %the top edge
				pass = pass+1;
				
				
				if BRinfo.faces(ii).top<0
					continue;
				end
				
				curr_edge = -10;
				if strcmp(BRinfo.faces(ii).system_top,'input_critical_curve')
					curr_edge = BRinfo.crit_curve.edges(BRinfo.faces(ii).top+1,:);
				elseif strcmp(BRinfo.faces(ii).system_top,'input_surf_sphere')
					curr_edge = BRinfo.sphere_curve.edges(BRinfo.faces(ii).top+1,:);
				else
					%do a lookup
					for zz = 1:length(BRinfo.singular_curves)
						if strcmp(BRinfo.singular_names{zz},BRinfo.faces(ii).system_top)
							curr_edge = BRinfo.singular_curves(zz).edges(BRinfo.faces(ii).top+1,:);
						end
					end
					
				end
				
				if curr_edge<0
					continue;
				end
				
				curr_edge = curr_edge([3 2 1]);
				
			case 2  %the bottom edge
				pass = pass+1;
				
				if BRinfo.faces(ii).bottom<0
					continue;
				end
				
				curr_edge = -10;
				if strcmp(BRinfo.faces(ii).system_bottom,'input_critical_curve')
					curr_edge = BRinfo.crit_curve.edges(BRinfo.faces(ii).bottom+1,:);
				elseif strcmp(BRinfo.faces(ii).system_bottom,'input_surf_sphere')
					curr_edge = BRinfo.sphere_curve.edges(BRinfo.faces(ii).bottom+1,:);
				else
					%do a lookup
					for zz = 1:length(BRinfo.singular_curves)
						if strcmp(BRinfo.singular_names{zz},BRinfo.faces(ii).system_bottom)
							curr_edge = BRinfo.singular_curves(zz).edges(BRinfo.faces(ii).bottom+1,:);
						end
					end
					
				end
				
				
				if curr_edge<0
					continue;
				end
				
			case 3  %the left edges
				if left_edge_counter <= BRinfo.faces(ii).num_left
					if BRinfo.faces(ii).left(left_edge_counter)<0 %an error check
						continue;
					end
					
					slice_ind = BRinfo.faces(ii).midslice_index+1; %offset by 1.
					edge_ind = BRinfo.faces(ii).left(left_edge_counter)+1; %offset by 1.
					
					curr_edge = BRinfo.critpoint_slices{slice_ind}.edges(edge_ind,:);
					left_edge_counter = left_edge_counter +1; %increment
					
					
				else
					pass = pass+1;
					continue;
				end
			case 4 %the right edges
				if right_edge_counter <= BRinfo.faces(ii).num_right
					
					if BRinfo.faces(ii).right(right_edge_counter)<0
						continue;
					end
					
					slice_ind = BRinfo.faces(ii).midslice_index+2;
					edge_ind = BRinfo.faces(ii).right(right_edge_counter)+1;
					curr_edge = BRinfo.critpoint_slices{slice_ind}.edges(edge_ind,:);
					right_edge_counter = right_edge_counter +1;
					
					curr_edge = curr_edge([3 2 1]);
					
				else
					pass = pass+1;
					continue;
				end
			case 5
				break;
				
		end
		
		
		local.faces(local_face_index,:) = [curr_edge(1) curr_edge(2) BRinfo.faces(ii).midpoint+1];
		local.faces(local_face_index+1,:) = [curr_edge(2) curr_edge(3) BRinfo.faces(ii).midpoint+1];
		local_face_index = local_face_index+2;
		
		stl_faces(curr_face_index,:) = [curr_edge(1) curr_edge(2) BRinfo.faces(ii).midpoint+1];
		stl_faces(curr_face_index+1,:) = [curr_edge(2) curr_edge(3) BRinfo.faces(ii).midpoint+1];
		curr_face_index = curr_face_index+2;
		
		
% 		%midpoint of the  edge
% 		pt = BRinfo.vertices(curr_edge(2)).point(ind);
% 		triangle.x(2,curr_triangle:curr_triangle+1) = [pt(1) pt(1)];
% 		triangle.y(2,curr_triangle:curr_triangle+1) = [pt(2) pt(2)];
% 		triangle.z(2,curr_triangle:curr_triangle+1) = [pt(3) pt(3)];
% 		
% 		
% 		%left point of the  edge
% 		pt = BRinfo.vertices(curr_edge(1)).point(ind);
% 		triangle.x(3,curr_triangle) = pt(1);
% 		triangle.y(3,curr_triangle) = pt(2);
% 		triangle.z(3,curr_triangle) = pt(3);
% 		
% 		
% 		%right point of the  edge
% 		pt = BRinfo.vertices(curr_edge(3)).point(ind);
% 		triangle.x(3,curr_triangle+1) = pt(1);
% 		triangle.y(3,curr_triangle+1) = pt(2);
% 		triangle.z(3,curr_triangle+1) = pt(3);
% 		
% 		curr_triangle = curr_triangle+2;
	end
	
	
% 	triangle.x = real(triangle.x);
% 	triangle.y = real(triangle.y);
% 	triangle.z = real(triangle.z);
	local.faces = local.faces(1:local_face_index-1,:);

	try
		plot_params.handles.faces(ii) = patch(local,'FaceColor',colors(ii,:),'FaceAlpha',0.3,'EdgeColor',0.985*colors(ii,:),'EdgeAlpha',0.2,'Parent',curr_axis);%,'EdgeAlpha',0.05 [0 0.9 0.9]
	catch
		ii
		local_face_index
		local.faces
		pause
	end
	
	
	
end

switch length(ind)
	case 2
		plot_params.handles.face_labels = text(pos(:,1),pos(:,2),txt,'Parent',curr_axis,'HorizontalAlignment','center','VerticalAlignment','top');

	case 3
		plot_params.handles.face_labels = text(pos(:,1),pos(:,2),pos(:,3),txt,'Parent',curr_axis,'HorizontalAlignment','center','VerticalAlignment','top');

	otherwise
			error('length of ind is not 2 or 3...')
end
	
set(plot_params.handles.face_labels,'visible','off');
end






