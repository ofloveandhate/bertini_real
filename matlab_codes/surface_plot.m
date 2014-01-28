%  surface specific code



function [fv] = surface_plot(sampler_data, BRinfo,ind)
global plot_params




create_axes_surface();

label_axes(ind,BRinfo,plot_params.axes.vertices);


fv.vertices = plot_vertices(ind, BRinfo);


plot_surface_edges(BRinfo);


fv.faces = plot_faces(BRinfo, ind);


plot_projection(BRinfo,ind);


sync_axes();


end












function plot_surface_edges(BRinfo)
global plot_params
%
line_thickness = 0.5;

crit_curve = zeros(3, BRinfo.num_variables-1, BRinfo.crit_curve.num_edges);

handle_counter = 0;
plot_params.handles.critcurve_labels = [];
plot_params.handles.spherecurve_labels = [];

plot_params.handles.critcurve = [];
critcurve_counter = 0;

if BRinfo.crit_curve.num_edges>0
	
	colors = 0.8*jet(BRinfo.crit_curve.num_edges);
	
	
	

	for ii =1:BRinfo.crit_curve.num_edges
		
		for jj = 1:3
			crit_curve(jj,:,ii) = real(transpose(BRinfo.vertices(BRinfo.crit_curve.edges(ii,jj)).point(1:BRinfo.num_variables-1)));
		end

		h = plot3(crit_curve(:,1,ii),crit_curve(:,2,ii),crit_curve(:,3,ii),'Parent',plot_params.axes.crit_curve);
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
			'Parent',plot_params.axes.crit_curve,'Color','r');
		plot_params.handles.critcurve_labels = [plot_params.handles.critcurve_labels; new_handles];
		

	end
end

plot_params.handles.spherecurve = [];
spherecurve_counter = 0;
if isfield(BRinfo,'sphere_curve')
if BRinfo.sphere_curve.num_edges>0
	sphere_curve = zeros(3, BRinfo.num_variables-1, BRinfo.sphere_curve.num_edges);

	colors = 0.8*jet(BRinfo.sphere_curve.num_edges);
	
	for ii =1:BRinfo.sphere_curve.num_edges
		
		for jj = 1:3
			sphere_curve(jj,:,ii) = real(transpose(BRinfo.vertices(BRinfo.sphere_curve.edges(ii,jj)).point(1:BRinfo.num_variables-1)));
		end

		h = plot3(sphere_curve(:,1,ii),sphere_curve(:,2,ii),sphere_curve(:,3,ii),'Parent',plot_params.axes.sphere_curve);
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
			'Parent',plot_params.axes.sphere_curve,'Color','c');
		plot_params.handles.spherecurve_labels = [plot_params.handles.spherecurve_labels; new_handles];
		

	end
end
end




colors = 0.7*jet(length(BRinfo.midpoint_slices));

plot_params.handles.midtext = [];

plot_params.handles.midslices = [];
midslice_counter = 0;
for kk = 1:length(BRinfo.midpoint_slices)
		midslice = zeros(3, BRinfo.num_variables-1, BRinfo.midpoint_slices(kk).num_edges);
		
		text_positions = zeros(3,BRinfo.midpoint_slices(kk).num_edges);
		textme = cell(1,BRinfo.midpoint_slices(kk).num_edges);
		
		
		for ii =1:BRinfo.midpoint_slices(kk).num_edges
			for jj = 1:3
				midslice(jj,:,ii) = real(transpose(BRinfo.vertices(BRinfo.midpoint_slices(kk).edges(ii,jj)).point(1:BRinfo.num_variables-1)));
			end
			
			h = plot3(midslice(:,1,ii),midslice(:,2,ii),midslice(:,3,ii),...
				'Parent',plot_params.axes.midedge);
			set(h,'Color',colors(kk,:));
			set(h,'LineStyle',':','LineWidth',line_thickness);

			text_positions(:,ii) = midslice(2,1:3,ii);
			textme{ii} = ['mid ' num2str(kk-1) '.' num2str(ii-1) '  '];
			
            midslice_counter = midslice_counter + 1;
            plot_params.handles.midslices(midslice_counter) = h;
		
		end
		
		if kk==1
			handle_counter = handle_counter+1;
			plot_params.legend.surface_edges.handles(handle_counter) = h;
			plot_params.legend.surface_edges.text{handle_counter} = 'midslice';
		end
		
				
		plot_params.handles.midtext = [plot_params.handles.midtext;...
			text(text_positions(1,:),text_positions(2,:),text_positions(3,:), textme,'HorizontalAlignment','right','FontSize',plot_params.fontsize-4,...
			'Parent',plot_params.axes.midedge,'Color','g')];


end






plot_params.handles.crittext = [];
colors = 0.7*jet(length(BRinfo.critpoint_slices));
plot_params.handles.critslices = [];
critslice_counter = 0;
for kk = 1:length(BRinfo.critpoint_slices)
		critslice = zeros(3, BRinfo.num_variables-1, BRinfo.critpoint_slices(kk).num_edges);
		
		text_positions = zeros(3,BRinfo.critpoint_slices(kk).num_edges);
		textme = cell(1,BRinfo.critpoint_slices(kk).num_edges);
		for ii =1:BRinfo.critpoint_slices(kk).num_edges
			for jj = 1:3
				critslice(jj,:,ii) = real(transpose(BRinfo.vertices(BRinfo.critpoint_slices(kk).edges(ii,jj)).point(1:BRinfo.num_variables-1)));
			end
			
			h = plot3(critslice(:,1,ii),critslice(:,2,ii),critslice(:,3,ii),'Parent',plot_params.axes.critedge);
			set(h,'Color',colors(kk,:));
			set(h,'LineStyle','--','LineWidth',line_thickness);

			text_positions(:,ii) = critslice(2,1:3,ii);
			textme{ii} = ['crit ' num2str(kk-1) '.' num2str(ii-1) '  '];
			
            critslice_counter = critslice_counter + 1;
			plot_params.handles.critslices(critslice_counter) = h;
            
		end
		if kk==1
			handle_counter = handle_counter+1;
			plot_params.legend.surface_edges.handles(handle_counter) = h;
			plot_params.legend.surface_edges.text{handle_counter} = 'critslice';
		end
		
		plot_params.handles.crittext = [plot_params.handles.crittext;...
			text(text_positions(1,:),text_positions(2,:),text_positions(3,:), textme,'HorizontalAlignment','right',...
			'FontSize',plot_params.fontsize-4,'Parent',plot_params.axes.critedge,'Color','b')];

				
				
end



end

















function stl_faces = plot_faces(BRinfo, ind)
global plot_params
%
num_total_faces = 0;
for ii = 1:BRinfo.num_faces
    curr_face = BRinfo.faces(ii);
    num_total_faces = num_total_faces + curr_face.num_left + curr_face.num_right + curr_face.top>=0 + curr_face.bottom>=0;
end
num_total_faces = num_total_faces*2;
stl_faces = zeros(num_total_faces, length(ind));
curr_face_index = 1;

curr_axis = plot_params.axes.faces;

txt = cell(BRinfo.num_faces,1);
pos = zeros(BRinfo.num_faces,3);
plot_params.handles.faces = [];

for ii = 1:BRinfo.num_faces
	
	
	num_triangles= 2*(2 + BRinfo.faces(ii).num_left + BRinfo.faces(ii).num_right);
	
	triangle.x = zeros(3,num_triangles);
	triangle.y = zeros(3,num_triangles);
	triangle.z = zeros(3,num_triangles);
	cdata = ii * ones(3,num_triangles);
	


	% set the midpoint of the face for all triangles to be the first row
	pt = transpose(BRinfo.vertices(BRinfo.faces(ii).midpoint+1).point(ind));
	triangle.x(1,:) = pt(1);
	triangle.y(1,:) = pt(2); % columns are new triangles
	triangle.z(1,:) = pt(3);
	
	
	txt{ii} = ['\newline' num2str(ii-1)];
	pos(ii,:) = pt(1:3);
	
	
	curr_triangle = 1; % counter
	
	pass = 1; left_edge_counter = 1;  right_edge_counter = 1;
	while 1
        
        switch pass
			case 1  %the top edge
				pass = pass+1;
				
				
                if BRinfo.faces(ii).top<0
					continue;
                end
				
                if BRinfo.faces(ii).system_top == -1600
					curr_edge = BRinfo.crit_curve.edges(BRinfo.faces(ii).top+1,:);
				elseif BRinfo.faces(ii).system_top == -1599
					curr_edge = BRinfo.sphere_curve.edges(BRinfo.faces(ii).top+1,:);
				else
					curr_edge = -10;
                end
                curr_edge = curr_edge([3 2 1]);
				if curr_edge<0
					continue;
				end
				
			case 2  %the bottom edge
				pass = pass+1;
				
				if BRinfo.faces(ii).bottom<0
					continue;
				end
				
				if BRinfo.faces(ii).system_bottom == -1600
					curr_edge = BRinfo.crit_curve.edges(BRinfo.faces(ii).bottom+1,:);
				elseif BRinfo.faces(ii).system_bottom == -1599
					curr_edge = BRinfo.sphere_curve.edges(BRinfo.faces(ii).bottom+1,:);
				else
					curr_edge = -10;
				end
				
				
				if curr_edge<0
					continue;
				end
				
			case 3
				if left_edge_counter <= BRinfo.faces(ii).num_left
                    if BRinfo.faces(ii).left(left_edge_counter)<0 %an error check
						continue;
                    end
                    BRinfo.faces(ii)
					slice_ind = BRinfo.faces(ii).midslice_index+1;
					edge_ind = BRinfo.faces(ii).left(left_edge_counter)+1;

                    curr_edge = BRinfo.critpoint_slices(slice_ind).edges(edge_ind,:);
					left_edge_counter = left_edge_counter +1;
                    
                    
				else
					pass = pass+1;
					continue;
				end
			case 4
				if right_edge_counter <= BRinfo.faces(ii).num_right
					
					if BRinfo.faces(ii).right(right_edge_counter)<0
						continue;
					end
					
					slice_ind = BRinfo.faces(ii).midslice_index+1 +1;
					edge_ind = BRinfo.faces(ii).right(right_edge_counter)+1;
					curr_edge = BRinfo.critpoint_slices(slice_ind).edges(edge_ind,:);
					right_edge_counter = right_edge_counter +1;
                    
                    curr_edge = curr_edge([3 2 1]);
                    
				else
					pass = pass+1;
					continue;
				end
			case 5
				break;
				
        end

        
        stl_faces(curr_face_index,:) = [curr_edge(1) curr_edge(2) BRinfo.faces(ii).midpoint+1];
        stl_faces(curr_face_index+1,:) = [curr_edge(2) curr_edge(3) BRinfo.faces(ii).midpoint+1];
		curr_face_index = curr_face_index+2;
        
        
        %midpoint of the  edge
		pt = BRinfo.vertices(curr_edge(2)).point(ind);
		triangle.x(2,curr_triangle:curr_triangle+1) = [pt(1) pt(1)];
		triangle.y(2,curr_triangle:curr_triangle+1) = [pt(2) pt(2)];
		triangle.z(2,curr_triangle:curr_triangle+1) = [pt(3) pt(3)];


		%left point of the  edge
		pt = BRinfo.vertices(curr_edge(1)).point(ind);
		triangle.x(3,curr_triangle) = pt(1);
		triangle.y(3,curr_triangle) = pt(2);
		triangle.z(3,curr_triangle) = pt(3);


		%right point of the  edge
		pt = BRinfo.vertices(curr_edge(3)).point(ind);
		triangle.x(3,curr_triangle+1) = pt(1);
		triangle.y(3,curr_triangle+1) = pt(2);
		triangle.z(3,curr_triangle+1) = pt(3);

		curr_triangle = curr_triangle+2;
	end


	triangle.x = real(triangle.x);
	triangle.y = real(triangle.y);
	triangle.z = real(triangle.z);

	plot_params.handles.faces(ii) = patch(triangle.x,triangle.y,triangle.z,cdata,'FaceAlpha',0.5,'EdgeColor','none','Parent',curr_axis);%,'EdgeAlpha',0.05 [0 0.9 0.9]
	
		
		
	
end

plot_params.handles.face_labels = text(pos(:,1),pos(:,2),pos(:,3),txt,'Parent',curr_axis,'HorizontalAlignment','center','VerticalAlignment','top');
set(plot_params.handles.face_labels,'visible','off');
end


function create_axes_surface()
global plot_params



new_axes = axes('Parent',plot_params.figures.main_plot);
delete( get(new_axes,'Children') );
set(new_axes,'visible','on');
hold(new_axes,'on');
plot_params.axes.vertices = new_axes;


new_axes = copyobj(plot_params.axes.vertices,plot_params.figures.main_plot);
delete( get(new_axes,'Children') );
set(new_axes,'visible','off');
hold(new_axes,'on')
plot_params.axes.projection = new_axes;
	
	
new_axes = copyobj(plot_params.axes.vertices,plot_params.figures.main_plot);
delete( get(new_axes,'Children') );
set(new_axes,'visible','off');
hold(new_axes,'on');
plot_params.axes.critedge = new_axes;
	

new_axes = copyobj(plot_params.axes.vertices,plot_params.figures.main_plot);
delete( get(new_axes,'Children') );
set(new_axes,'visible','off');
hold(new_axes,'on');
plot_params.axes.crit_curve = new_axes;

new_axes = copyobj(plot_params.axes.vertices,plot_params.figures.main_plot);
delete( get(new_axes,'Children') );
set(new_axes,'visible','off');
hold(new_axes,'on');
plot_params.axes.sphere_curve = new_axes;


new_axes = copyobj(plot_params.axes.vertices,plot_params.figures.main_plot);
delete( get(new_axes,'Children') );
set(new_axes,'visible','off');
hold(new_axes,'on');
plot_params.axes.midedge = new_axes;


new_axes = copyobj(plot_params.axes.vertices,plot_params.figures.main_plot);
delete( get(new_axes,'Children') );
set(new_axes,'visible','off');
hold(new_axes,'on');
plot_params.axes.faces = new_axes;

end


