

function br_plotter = plot_surface_edges(br_plotter)




handle_counter = 0;


if br_plotter.BRinfo.crit_curve.num_edges>0
	
	if isfield(br_plotter.options,'use_fixed_linestyle')
			style = br_plotter.options.linestyle;
		else
			style = '-';
	end
		
	[br_plotter.handles.curves.raw.crit, ...
		br_plotter.handles.curves.refinements.crit, ...
		br_plotter.handles.critcurve_labels] ...
		  = plot_subcurve(br_plotter,br_plotter.BRinfo.crit_curve,'critcurve',style,'r');
	
	

	if ~isempty(br_plotter.handles.curves.raw.crit)
		handle_counter = handle_counter+1;
		br_plotter.legend.surface_edges.handles(handle_counter) = br_plotter.handles.curves.raw.crit(1);
		br_plotter.legend.surface_edges.text{handle_counter} = 'critical curve';
	end
end










if isfield(br_plotter.BRinfo,'sphere_curve')
	if br_plotter.BRinfo.sphere_curve.num_edges>0
		
		if isfield(br_plotter.options,'use_fixed_linestyle')
			style = br_plotter.options.linestyle;
		else
			style = '-.';
		end
		
		[br_plotter.handles.curves.raw.sphere, ...
			br_plotter.handles.curves.refinements.sphere, ...
			br_plotter.handles.spherecurve_labels] ...
			  = plot_subcurve(br_plotter,br_plotter.BRinfo.sphere_curve,'spherecurve',style,'c');
	
		if ~isempty(br_plotter.handles.curves.raw.sphere)
			handle_counter = handle_counter+1;
			br_plotter.legend.surface_edges.handles(handle_counter) = br_plotter.handles.curves.raw.sphere(1);
			br_plotter.legend.surface_edges.text{handle_counter} = 'sphere curve';
		end
	end
end



num_midslices = length(br_plotter.BRinfo.midpoint_slices);
num_crit_slices = length(br_plotter.BRinfo.critpoint_slices);

colors = br_plotter.options.colormap(2*num_midslices+1);






firstone = 1;

if isfield(br_plotter.options,'use_fixed_linestyle')
	style = br_plotter.options.linestyle;
else
	style = ':';
end
		
for kk = 1:num_crit_slices
	color_index = 1+2*(kk-1);
	
	[handie_mc_handhand, refinement_handles, label_handles] ...
		= plot_subcurve(br_plotter,br_plotter.BRinfo.critpoint_slices{kk},sprintf('crit.%d.',kk),style,'m',colors(color_index,:));
	
	br_plotter.handles.crittext = [br_plotter.handles.crittext;label_handles];
	br_plotter.handles.curves.raw.critslices = [br_plotter.handles.curves.raw.critslices;handie_mc_handhand];
	br_plotter.handles.curves.refinements.critslices = [br_plotter.handles.curves.refinements.critslices;refinement_handles];
	
	if firstone
		if ~isempty(handie_mc_handhand)
			handle_counter = handle_counter+1;
			br_plotter.legend.surface_edges.handles(handle_counter) = handie_mc_handhand(1);
			br_plotter.legend.surface_edges.text{handle_counter} = 'crit slices';
			firstone = 0;
		end
	end
	
	
end






if isfield(br_plotter.options,'use_fixed_linestyle')
	style = br_plotter.options.linestyle;
else
	style = '--';
end
		

added = false;
for kk = 1:length(br_plotter.BRinfo.midpoint_slices)
	if strcmp( br_plotter.BRinfo.midpoint_slices{kk}.inputfilename,'unset')
		continue;
	end
	color_index = 2*(kk);
	[handie_mc_handhand, refinement_handles, label_handles] ...
		= plot_subcurve(br_plotter,br_plotter.BRinfo.midpoint_slices{kk},sprintf('mid.%d.',kk),style,'g',colors(color_index,:));
	
	br_plotter.handles.midtext = [br_plotter.handles.midtext;label_handles];
	br_plotter.handles.curves.raw.midslices = [br_plotter.handles.curves.raw.midslices;handie_mc_handhand];
	br_plotter.handles.curves.refinements.midslices = [br_plotter.handles.curves.refinements.midslices;refinement_handles];
	
	if and(added==false,~isempty(handie_mc_handhand))
		handle_counter = handle_counter+1;
		br_plotter.legend.surface_edges.handles(handle_counter) = handie_mc_handhand(1);
		br_plotter.legend.surface_edges.text{handle_counter} = 'mid slices';
		added = true;
	end
	
	
end









if isfield(br_plotter.options,'use_fixed_linestyle')
	style = br_plotter.options.linestyle;
else
	style = ':';
end
		
if isfield(br_plotter.BRinfo,'singular_curves')
	
	
	firstone = 1;

	for kk = 1:length(br_plotter.BRinfo.singular_curves)
		
		
		[handie_mc_handhand, refinement_handles, label_handles] = plot_subcurve(br_plotter,br_plotter.BRinfo.singular_curves{kk},sprintf('sing.%d.',kk),style,'b');
	

		br_plotter.handles.singtext = [br_plotter.handles.singtext;label_handles];
		br_plotter.handles.curves.raw.singular = [br_plotter.handles.curves.raw.singular;handie_mc_handhand];
		br_plotter.handles.curves.refinements.singular = [br_plotter.handles.curves.refinements.singular;refinement_handles];

		if firstone
			if ~isempty(handie_mc_handhand)
				handle_counter = handle_counter+1;
				br_plotter.legend.surface_edges.handles(handle_counter) = handie_mc_handhand(1);
				br_plotter.legend.surface_edges.text{handle_counter} = 'singular curves';
				firstone = 0;
			end
		end
	end
end




end








%plots a single surface subcurve
function [edge_handles, refinement_handles, text_handle] = plot_subcurve(br_plotter,curve,name,style,text_color,desiredcolor)

ind = br_plotter.indices;
curr_axes = br_plotter.axes.main;


edge_handles = [];
refinement_handles = [];
text_handle = [];

if br_plotter.options.touching_edges_only
	num_touching = 0;
	touching_edge_indices = zeros(curve.num_edges,1); %preallocate
	for ii = 1:curve.num_edges
		if edge_touches_faces(br_plotter, ii, curve)
			num_touching = num_touching+1;
			touching_edge_indices(num_touching) = ii;
		end
	end
	touching_edge_indices = touching_edge_indices(1:num_touching); %trim the fat
else
	num_touching = curve.num_edges;
	touching_edge_indices = 1:num_touching;
end


num_nondegen = 0;
nondegen_edge_indices = zeros(num_touching,1); %preallocate
for ii = 1:length(touching_edge_indices)
	curr_edge_index = touching_edge_indices(ii);
	if curve.edges(curr_edge_index,1)~=curve.edges(curr_edge_index,3)
		num_nondegen = num_nondegen+1;
		nondegen_edge_indices(num_nondegen) = curr_edge_index;
	end
end
nondegen_edge_indices = nondegen_edge_indices(1:num_nondegen); %trim the fat



if num_nondegen==0
	return;
end



if nargin==6
	colors = repmat(desiredcolor,[num_nondegen 1]);
else
	colors = 0.8*br_plotter.options.colormap(num_nondegen);
end


% edge_handles = zeros(num_nondegen,1);
if br_plotter.options.labels
	text_locations = zeros(num_nondegen,3);
	text_labels = cell(num_nondegen,1);
end
curve_edge_points = zeros(0,3);
for ii =1:num_nondegen
	curr_edge_index = nondegen_edge_indices(ii);
	curr_edge = curve.edges(curr_edge_index,:);
	
	temp = br_plotter.data.space.vertices(curr_edge,ind);
	curve_edge_points = [curve_edge_points; nan nan nan; temp];
	
	if br_plotter.options.labels
		text_locations(ii,:) = temp(2,:);
		text_labels{ii} = sprintf('%s %d  ',name,curr_edge_index);
	end

end



h = plot3(curve_edge_points(:,1),curve_edge_points(:,2),curve_edge_points(:,3),'Parent',curr_axes);

set(h,'Color',colors(ii,:));
set(h,'LineStyle',style,'LineWidth',br_plotter.options.linewidth);
edge_handles = h;


	
	
if br_plotter.options.labels
	text_handle = text(text_locations(:,1),text_locations(:,2),text_locations(:,3), text_labels,...
				'HorizontalAlignment','right',...
				'FontSize',br_plotter.options.fontsizes.labels,...
				'Parent',curr_axes,'Color',text_color);
end


if and(~isempty(br_plotter.BRinfo.sampler_data),br_plotter.options.render_samples)
	if nargin==6
		refinement_handles = plot_curve_samples(br_plotter,curve,style, desiredcolor);
	else
		refinement_handles = plot_curve_samples(br_plotter,curve,style);
	end
end


		

end





