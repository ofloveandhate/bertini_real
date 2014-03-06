function varargout = plot_br_samples()

clear all;
close all;

global plot_params;


if ~isempty(which('getmondim'))
	mondim = getmondim;
	plot_params.window = [20 20 mondim(3)-130 mondim(4)-130];
else
	plot_params.window = [20 20 1280 768];
end
	

plot_params.figures.supermain = figure(gcf);

set(plot_params.figures.supermain,'Position',plot_params.window)
set(plot_params.figures.supermain,'toolbar','figure');
plot_params.figures.main_plot = plot_params.figures.supermain;

plot_params.fontsize = 13;

	


cameratoolbar

button_setup();

load_and_render();


if nargout==1
	varargout{1} = plot_params;
end



end%re: function




function load_and_render(varargin)
global plot_params

BRinfo = [];
sampler_data = [];


if isempty(varargin)

	prev_filenames = dir('BRinfo*.mat');
	max_found = -1;

	for ii = 1:length(prev_filenames)
		curr_name = prev_filenames(ii).name;
		curr_num = str2num(curr_name(7:end-4));
		if max_found < curr_num
			max_found = curr_num;
		end

	end
	filename = ['BRinfo' num2str(max_found) '.mat'];
	load(filename);
	
	[containing, here, ~] = fileparts(pwd);
	plot_params.basename = here;


	
else
	
	
	
	[filename, pathname, filterindex] = uigetfile( ...
		{  '*.mat','MAT-files (*.mat)'; ...
		'*.*',  'All Files (*.*)'}, ...
		'Choose some data', ...
		'MultiSelect', 'off');

	if isnumeric(filename) %filename is 0 if user cancelled
		return;
	end
		

	load([pathname filename]);
	
	[containing, ~, ~] = fileparts(pathname);
	slashes = find(containing=='/');
	plot_params.basename = containing(slashes(end)+1:end);
	
	
	axes_names = fieldnames(plot_params.axes);
	
	for ii = 1:length(axes_names)
		delete(allchild(plot_params.axes.(axes_names{ii})))
		delete(plot_params.axes.(axes_names{ii}));
	end
	
	f = fieldnames(plot_params.checkboxes);
	for ii = 1:length(f)
		delete(plot_params.checkboxes.(f{ii}));
	end
	delete(plot_params.panels.visibility);
	
	
	plot_params = rmfield(plot_params,'axes');
	plot_params = rmfield(plot_params,'legend');
	plot_params = rmfield(plot_params,'handles'); 
	plot_params = rmfield(plot_params,'switches'); 
	plot_params = rmfield(plot_params,'checkboxes'); 
	
end

plot_params.dimension = BRinfo.dimension;


ind = get_indices_interactive(sampler_data,BRinfo);

plot_params.ind = ind;
switch BRinfo.dimension
	case 1
		curve_plot(BRinfo,ind);
		
	case 2	
		plot_params.faces_and_vertices = surface_plot(BRinfo,ind);
		
	otherwise
	
end



render_legends();

controls(BRinfo);

save_routine
end








function [plot_indices] = get_indices_interactive(sampler_data,BRinfo)


% set up the plotting indices -- which data to plot

% 		first scan for zero columns.
if ~isempty(sampler_data)
	tmpdata = zeros(sum(sampler_data.sample_sizes),BRinfo.num_variables-1);
	counter = 1;
	for ii = 1:BRinfo.num_edges
		for jj = 1:sampler_data.sample_sizes(ii)
			tmpdata(counter,:) = BRinfo.vertices(sampler_data.edge(ii).samples(jj)+1).point;
			counter = counter+1;
		end
	end
else
	tmpdata = zeros(BRinfo.num_vertices,BRinfo.num_variables-1);
	for ii = 1:BRinfo.num_vertices
		tmpdata(ii,:) = BRinfo.vertices(ii).point(1:BRinfo.num_variables-1);
	end
end

indices_of_nonconst_cols = find_constant_vars(tmpdata);

if length(indices_of_nonconst_cols)<4
	plot_indices = indices_of_nonconst_cols;
else
	plot_indices = get_user_indices(indices_of_nonconst_cols,BRinfo);
end

end





function ind = get_user_indices(suggested_indices, BRinfo)
%
display('the following variables are in the problem:');

for ii = 1:BRinfo.num_variables-1
	if find(suggested_indices==ii)
		emptystring = '';
	else
		emptystring = '(uniformly zero)';
	end
	display(sprintf('%i: %s %s',ii,BRinfo.var_names{ii},emptystring));
	
end

num_plot_vars = input('how many variables to plot? ');

ind = zeros(1,num_plot_vars);
for ii = 1:num_plot_vars
	ind(ii) = input(sprintf('var %i: ',ii));
end


end

% find any variables which have zero value uniformly.  should be
% improved to find constant variables.
function ind = find_constant_vars(data,dim)
%
zerothresh = 1e-7;

if nargin < 2 
	dim = 1;
end
biggies = abs(max(data)-min(data));
% biggies = max(abs(data),[],dim);
ind = find(biggies>zerothresh);

end




function controls(BRinfo)
%
global plot_params 



set_initial_visibility();
visibility_setup();

twiddle_visibility(plot_params);


if length(plot_params.ind)==3
    camera_setup(BRinfo);
end

end



function button_setup()
global plot_params


h = uipanel('units','pixels','position',[10 10 120 160],'visible','on');
plot_params.panels.buttons = h;


plot_params.buttons.save = uicontrol('Style', 'pushbutton', 'String', 'Save',...
        'Position', [10 10 100 20],...
        'Callback', {@save_routine},'Parent',h); 
	
plot_params.buttons.load = uicontrol('Style', 'pushbutton', 'String', 'Load & Render',...
        'Position', [10 40 100 20],...
        'Callback', {@load_and_render},'Parent',h); 

	
plot_params.buttons.center = uicontrol('Style', 'pushbutton', 'String', 'Center',...
        'Position', [10 70 100 20],...
        'Callback', {@center_camera_on_selected_point},'Parent',h); 

    plot_params.buttons.leap = uicontrol('Style', 'pushbutton', 'String', 'Leap',...
        'Position', [10 100 100 20],...
        'Callback', {@leap_figure_axes_control},'Parent',h); 
    
    plot_params.buttons.stl = uicontrol('Style', 'pushbutton', 'String', 'STL',...
        'Position', [10 130 100 20],...
        'Callback', {@fv_to_stl},'Parent',h); 
end



function camera_setup(BRinfo)
global plot_params 


% view(3);  %yep, switch to view 3 by default for initial view.
% set(plot_params.axes.vertices,'CameraViewAngle', 20);

set(gca,'CameraTargetMode','manual');
plot_params.scene.target = real(BRinfo.center(plot_params.ind));
set(gca,'CameraTarget',plot_params.scene.target);


init_campos = real(BRinfo.center(plot_params.ind));
init_campos = init_campos + 9*BRinfo.radius;

plot_params.scene.campos = init_campos;

set(gca,'CameraPosition',plot_params.scene.campos);
set(gca,'CameraViewAngle',10);
rotate3d off

end







%labels and set visibility

function set_initial_visibility()
global plot_params

plot_params.switches.legend = 0;

plot_params.switches.label_faces = 0;
plot_params.switches.label_spherecurve = 0;
plot_params.switches.label_critcurve = 0;
plot_params.switches.label_critedges = 0;
plot_params.switches.label_midedges = 0;
plot_params.switches.label_singular = 0;
plot_params.switches.show_edges = 0;
plot_params.switches.show_curve_samples = 0;

plot_params.switches.display_vertices = 0;
plot_params.switches.label_vertices = 0;

f = plot_params.legend.vertices.types;

for ii = 1:length(f)
    if strcmp('CRITICAL',f{ii})
        plot_params.switches.vertex_set.(f{ii}) = 1;
    else
        plot_params.switches.vertex_set.(f{ii}) = 0;
    end
end


if plot_params.dimension==2
%         plot_params.switches.display_faces = 1;
%         plot_params.switches.display_face_samples = 0;
    if isempty(plot_params.handles.surface_samples)
        plot_params.switches.display_faces = 1;
        plot_params.switches.display_face_samples = 0; 
    else
        plot_params.switches.display_faces = 0;
        plot_params.switches.display_face_samples = 1;
    end
        
	
else
    
    if isempty(plot_params.handles.sample_edges)
        plot_params.switches.show_edges = 1;
        plot_params.switches.show_curve_samples = 0; 
    else
        plot_params.switches.show_edges = 0;
        plot_params.switches.show_curve_samples = 1;
    end
    
    plot_params.switches.display_faces = 0;
    plot_params.switches.display_face_samples = 0;

end
plot_params.switches.curve_refinements = 0;
plot_params.switches.display_projection = 0;


plot_params.switches.main_axes = 1;

plot_params.switches.show_critcurve = 0;
plot_params.switches.show_spherecurve = 0;
plot_params.switches.show_critslices = 0;
plot_params.switches.show_midslices = 0;
plot_params.switches.show_singular = 0;
end


function visibility_setup()

global plot_params

button_pos = get(plot_params.panels.buttons,'Position');
init_y = button_pos(2)+button_pos(4)+5;


plot_params.panels.visibility = uipanel('units','pixels','Position',[20 init_y 120 400],'visible','on');

checkbox.y_start = 0;
checkbox.y_factor = 21;

checkbox.x = 0;
checkbox.w = 100;
checkbox.h = 20;

num_checkboxes = 0;








pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
plot_params.checkboxes.legend = uicontrol(...
	'style','checkbox','units','pixels','position',pos,...
	'String','legend','Value',plot_params.switches.legend,'callback',{@flip_switch,'legend'},...
	'Parent',plot_params.panels.visibility); 
num_checkboxes = num_checkboxes+1;





pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
plot_params.checkboxes.main_axes_visibility = uicontrol(...
	'style','checkbox','units','pixels','position',pos,...
	'String','main axes','Value',plot_params.switches.main_axes,'callback',{@flip_switch,'main_axes'},...
	'Parent',plot_params.panels.visibility); 
num_checkboxes = num_checkboxes+1;





% plot_params.panels.labels = uipanel('units','pixels','Position',[5 5+num_checkboxes*checkbox.y_factor 110 390],'visible','on','Parent',plot_params.panels.visibility);





pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
plot_params.checkboxes.vertex_label_visibility = uicontrol(...
	'style','checkbox','units','pixels','position',pos,...
	'String','vertex labels','Value',plot_params.switches.label_vertices,'callback',{@flip_switch,'label_vertices'},...
	'Parent',plot_params.panels.visibility); 
num_checkboxes = num_checkboxes+1;

pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
plot_params.checkboxes.vertex_marker_visibility = uicontrol(...
	'style','checkbox','units','pixels','position',pos,...
	'String','vertex markers','Value',plot_params.switches.display_vertices,'callback',{@flip_switch,'display_vertices'},...
	'Parent',plot_params.panels.visibility); 
num_checkboxes = num_checkboxes+1;


pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
plot_params.checkboxes.projection_visibility = uicontrol(...
	'style','checkbox','units','pixels','position',pos,...
	'String','projection','Value',plot_params.switches.display_projection,'callback',{@flip_switch,'display_projection'},...
	'Parent',plot_params.panels.visibility); 
num_checkboxes = num_checkboxes+1;


f = plot_params.legend.vertices.types;
for ii = 1:length(f)
	pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
	plot_params.checkboxes.vertex_set_visibility(ii) = uicontrol(...
		'style','checkbox','units','pixels','position',pos,...
		'String',f{ii},'Value',plot_params.switches.vertex_set.(f{ii}),...
		'callback',{@flip_switch,f{ii},'vertex_set'},...
		'Parent',plot_params.panels.visibility);

	new_color = get(plot_params.handles.vertex_text.(f{ii}));
	new_color = new_color.Color;
	set(plot_params.checkboxes.vertex_set_visibility(ii),...
		'ForeGroundColor',new_color(1,:));
	
	num_checkboxes = num_checkboxes+1;
end



switch plot_params.dimension
	case 1
		pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
		plot_params.checkboxes.critcurve_visibility = uicontrol(...
			'style','checkbox','units','pixels','position',pos,...
			'String','edges','Value',plot_params.switches.show_edges,'callback',{@flip_switch,'show_edges'},...
			'Parent',plot_params.panels.visibility,...
			'ForeGroundColor','k'); 
		num_checkboxes = num_checkboxes+1;
        
        if ~isempty(plot_params.handles.sample_edges)
            pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
            plot_params.checkboxes.critcurve_visibility = uicontrol(...
                'style','checkbox','units','pixels','position',pos,...
                'String','refinement','Value',plot_params.switches.show_curve_samples,'callback',{@flip_switch,'show_curve_samples'},...
                'Parent',plot_params.panels.visibility,...
                'ForeGroundColor','k'); 
            num_checkboxes = num_checkboxes+1;
        end
        
	case 2
		
		pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
		plot_params.checkboxes.critcurve_visibility = uicontrol(...
			'style','checkbox','units','pixels','position',pos,...
			'String','face labels','Value',plot_params.switches.label_faces,'callback',{@flip_switch,'label_faces'},...
			'Parent',plot_params.panels.visibility,...
			'ForeGroundColor','k'); 
		num_checkboxes = num_checkboxes+1;
		
		
        pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
		plot_params.checkboxes.face_visibility = uicontrol(...
			'style','checkbox','units','pixels','position',pos,...
			'String','faces','Value',plot_params.switches.display_faces,'callback',{@flip_switch,'display_faces'},...
			'Parent',plot_params.panels.visibility); 
		num_checkboxes = num_checkboxes+1;
        
        pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
		plot_params.checkboxes.midslice_visibility = uicontrol(...
			'style','checkbox','units','pixels','position',pos,...
			'String','face samples','Value',plot_params.switches.display_face_samples,'callback',{@flip_switch,'display_face_samples'},...
			'Parent',plot_params.panels.visibility); 
		num_checkboxes = num_checkboxes+1;
        
        
        
		
        if ~isempty(plot_params.handles.critcurve_labels)
            new_color = get(plot_params.handles.critcurve_labels);
            new_color = new_color.Color;
        else
            new_color = [0 0 0];
        end
	
		pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
		plot_params.checkboxes.critcurve_visibility = uicontrol(...
			'style','checkbox','units','pixels','position',pos,...
			'String','critcurve labels','Value',plot_params.switches.label_critcurve,'callback',{@flip_switch,'label_critcurve'},...
			'Parent',plot_params.panels.visibility,...
			'ForeGroundColor',new_color(1,:)); 
		num_checkboxes = num_checkboxes+1;

		
		
		if ~isempty(plot_params.handles.spherecurve_labels)
			new_color = get(plot_params.handles.spherecurve_labels);
			new_color = new_color.Color;
		else
			new_color = [0 0 0];
		end
		
		pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
		plot_params.checkboxes.critcurve_visibility = uicontrol(...
			'style','checkbox','units','pixels','position',pos,...
			'String','spherecurve labels','Value',plot_params.switches.label_spherecurve,'callback',{@flip_switch,'label_spherecurve'},...
			'Parent',plot_params.panels.visibility,...
			'ForeGroundColor',new_color(1,:)); 
		num_checkboxes = num_checkboxes+1;
		
		
		if ~isempty(plot_params.handles.midtext)
			new_color = get(plot_params.handles.midtext);
			new_color = new_color.Color;
		else
			new_color = [0 0 0];
		end
		pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
		plot_params.checkboxes.midedge_visibility = uicontrol(...
			'style','checkbox','units','pixels','position',pos,...
			'String','midedge labels','Value',plot_params.switches.label_midedges,'callback',{@flip_switch,'label_midedges'},...
			'Parent',plot_params.panels.visibility,...
			'ForeGroundColor',new_color(1,:)); 
		num_checkboxes = num_checkboxes+1;

		
		if ~isempty(plot_params.handles.crittext)
			new_color = get(plot_params.handles.crittext);
			new_color = new_color.Color;
		else
			new_color = [0 0 0];
		end
		

		
		pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
		plot_params.checkboxes.critedge_visibility = uicontrol(...
			'style','checkbox','units','pixels','position',pos,...
			'String','critedge labels','Value',plot_params.switches.label_critedges,'callback',{@flip_switch,'label_critedges'},...
			'Parent',plot_params.panels.visibility,...
			'ForeGroundColor',new_color(1,:)); 
		num_checkboxes = num_checkboxes+1;
		
		


        
        if ~isempty(plot_params.handles.singtext)
			new_color = get(plot_params.handles.singtext);
			new_color = new_color.Color;
		else
			new_color = [0 0 0];
		end
		

		
		pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
		plot_params.checkboxes.singular_visibility = uicontrol(...
			'style','checkbox','units','pixels','position',pos,...
			'String','singedge labels','Value',plot_params.switches.label_singular,'callback',{@flip_switch,'label_singular'},...
			'Parent',plot_params.panels.visibility,...
			'ForeGroundColor',new_color(1,:)); 
		num_checkboxes = num_checkboxes+1;
        
        
        
        
        pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
		plot_params.checkboxes.critcurve_visibility = uicontrol(...
			'style','checkbox','units','pixels','position',pos,...
			'String','crit curve','Value',plot_params.switches.show_critcurve,'callback',{@flip_switch,'show_critcurve'},...
			'Parent',plot_params.panels.visibility); 
		num_checkboxes = num_checkboxes+1;
        
        pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
		plot_params.checkboxes.spherecurve_visibility = uicontrol(...
			'style','checkbox','units','pixels','position',pos,...
			'String','sphere curve','Value',plot_params.switches.show_spherecurve,'callback',{@flip_switch,'show_spherecurve'},...
			'Parent',plot_params.panels.visibility); 
		num_checkboxes = num_checkboxes+1;
        
        pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
		plot_params.checkboxes.critslice_visibility = uicontrol(...
			'style','checkbox','units','pixels','position',pos,...
			'String','critslices','Value',plot_params.switches.show_critslices,'callback',{@flip_switch,'show_critslices'},...
			'Parent',plot_params.panels.visibility); 
		num_checkboxes = num_checkboxes+1;
        
        pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
		plot_params.checkboxes.midslice_visibility = uicontrol(...
			'style','checkbox','units','pixels','position',pos,...
			'String','midslices','Value',plot_params.switches.show_midslices,'callback',{@flip_switch,'show_midslices'},...
			'Parent',plot_params.panels.visibility); 
		num_checkboxes = num_checkboxes+1;


        pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
		plot_params.checkboxes.singular_visibility = uicontrol(...
			'style','checkbox','units','pixels','position',pos,...
			'String','singular curves','Value',plot_params.switches.show_singular,'callback',{@flip_switch,'show_singular'},...
			'Parent',plot_params.panels.visibility); 
		num_checkboxes = num_checkboxes+1;

		pos =[checkbox.x,checkbox.y_start+checkbox.y_factor*num_checkboxes,checkbox.w,checkbox.h];
		plot_params.checkboxes.refinement_visibility = uicontrol(...
			'style','checkbox','units','pixels','position',pos,...
			'String','curve refinements','Value',plot_params.switches.curve_refinements,'callback',{@flip_switch,'curve_refinements'},...
			'Parent',plot_params.panels.visibility); 
		num_checkboxes = num_checkboxes+1;

		%end 3d case
	otherwise
		
end


set(plot_params.panels.visibility,'units','pixels','position',[10 init_y 120 checkbox.y_factor*num_checkboxes])
	
	
end


function twiddle_visibility(plot_params)



if plot_params.switches.legend == 0
	set(plot_params.handles.legend,'visible','off');
else
	set(plot_params.handles.legend,'visible','on');
end


if plot_params.switches.main_axes == 0
	set(plot_params.axes.main,'visible','off');
else
	set(plot_params.axes.main,'visible','on');
end


%vertices
if plot_params.switches.display_vertices == 0
	f = fieldnames(plot_params.switches.vertex_set);
	for ii = 1:length(f)
		set(plot_params.handles.vertices.(f{ii}),'visible','off');
		set(plot_params.handles.vertex_text.(f{ii}),'visible','off');
	end
	
else

	f = fieldnames(plot_params.switches.vertex_set);
	for ii = 1:length(f)
		if plot_params.switches.vertex_set.(f{ii}) == 0
			set(plot_params.handles.vertices.(f{ii}),'visible','off')
		else
			set(plot_params.handles.vertices.(f{ii}),'visible','on')
		end
	end

	for ii = 1:length(f)
		if plot_params.switches.label_vertices == 0;

				set(plot_params.handles.vertex_text.(f{ii}),'visible','off')

		else
			if plot_params.switches.vertex_set.(f{ii}) == 0
				set(plot_params.handles.vertex_text.(f{ii}),'visible','off');
			else
				set(plot_params.handles.vertex_text.(f{ii}),'visible','on');
			end
		end
	end
end



if plot_params.switches.display_projection ==0
	set(plot_params.handles.projection(:),'visible','off');
else
	set(plot_params.handles.projection(:),'visible','on');
end
	
	
	
	
	
	
	
	
	
	
	
	switch plot_params.dimension
		case 1
			
			
            if plot_params.switches.show_edges == 0
                set(plot_params.handles.edges(:),'visible','off')
            else
                set(plot_params.handles.edges(:),'visible','on')
            end
			
			if plot_params.switches.show_curve_samples == 0
                set(plot_params.handles.sample_edges(:),'visible','off')
            else
                set(plot_params.handles.sample_edges(:),'visible','on')
            end
			
			
		case 2
            
            
            
            if plot_params.switches.show_critcurve == 0
                set(plot_params.handles.critcurve(:),'visible','off')
                set(plot_params.handles.critcurve_labels,'visible','off');
            else
                set(plot_params.handles.critcurve(:),'visible','on')
                if (plot_params.switches.label_critcurve == 0)
                    set(plot_params.handles.critcurve_labels,'visible','off');
                else
                    set(plot_params.handles.critcurve_labels,'visible','on');
                end
            end
            
            
            
            if plot_params.switches.show_spherecurve == 0
                set(plot_params.handles.spherecurve(:),'visible','off')
                set(plot_params.handles.spherecurve_labels,'visible','off');
            else
                set(plot_params.handles.spherecurve(:),'visible','on')
                if (plot_params.switches.label_spherecurve == 0)
                    set(plot_params.handles.spherecurve_labels,'visible','off');
                else
                    set(plot_params.handles.spherecurve_labels,'visible','on');
                end
            end
            
            if plot_params.switches.show_critslices == 0
                set(plot_params.handles.critslices(:),'visible','off')
                set(plot_params.handles.crittext,'visible','off');
            else
                
                set(plot_params.handles.critslices(:),'visible','on')

                if (plot_params.switches.label_critedges == 0)
                    set(plot_params.handles.crittext,'visible','off');
                else
                    set(plot_params.handles.crittext,'visible','on');
                end
            end
            
            
            
            if plot_params.switches.show_singular == 0
                set(plot_params.handles.singular_curves(:),'visible','off')
                set(plot_params.handles.singtext,'visible','off');
            else
                set(plot_params.handles.singular_curves(:),'visible','on')
                
                if (plot_params.switches.label_singular == 0)
                    set(plot_params.handles.singtext,'visible','off');
                else
                    set(plot_params.handles.singtext,'visible','on');
                end
            
                
                
            end
            
            
            if plot_params.switches.show_midslices == 0
                set(plot_params.handles.midslices(:),'visible','off')
                set(plot_params.handles.midtext,'visible','off');
            else
                set(plot_params.handles.midslices(:),'visible','on')
                
                if (plot_params.switches.label_midedges == 0)
                    set(plot_params.handles.midtext,'visible','off');
                else
                    set(plot_params.handles.midtext,'visible','on');
                end
            
                
                
            end
            
            
			if plot_params.switches.label_faces == 0
				set(plot_params.handles.face_labels,'visible','off');
			else
				set(plot_params.handles.face_labels,'visible','on');
			end
			
			if plot_params.switches.curve_refinements == 0
				set(plot_params.handles.refinements.critslice(:),'visible','off');
				set(plot_params.handles.refinements.midslice(:),'visible','off');
				set(plot_params.handles.refinements.spherecurve,'visible','off');
				set(plot_params.handles.refinements.critcurve,'visible','off');
				set(plot_params.handles.refinements.singularcurve(:),'visible','off');
			else
				set(plot_params.handles.refinements.critslice(:),'visible','on');
				set(plot_params.handles.refinements.midslice(:),'visible','on');
				set(plot_params.handles.refinements.spherecurve,'visible','on');
				set(plot_params.handles.refinements.critcurve,'visible','on');
				set(plot_params.handles.refinements.singularcurve(:),'visible','on');
			end




			

			
			
			
			%faces
			if plot_params.switches.display_faces == 0
				for ii = 1:length(plot_params.handles.faces)
					set(plot_params.handles.faces(ii),'visible','off');
				end
			else
				for ii = 1:length(plot_params.handles.faces)
					set(plot_params.handles.faces(ii),'visible','on');
				end
            end
            
            %faces
			if plot_params.switches.display_face_samples == 0
				for ii = 1:length(plot_params.handles.surface_samples)
					set(plot_params.handles.surface_samples(ii),'visible','off');
				end
			else
				for ii = 1:length(plot_params.handles.surface_samples)
					set(plot_params.handles.surface_samples(ii),'visible','on');
				end
            end
	
            
		otherwise
			
	end
	

	
	
	

	
	
	
	
	
end

function flip_switch(varargin)
global plot_params
% length(varargin)

switch_name = varargin{3};

if length(varargin)==3
	plot_params.switches.(switch_name) = mod(plot_params.switches.(switch_name)+1,2);
elseif length(varargin)>3
	plot_params.switches.(varargin{4}).(switch_name) = mod(plot_params.switches.(varargin{4}).(switch_name)+1,2);
end




twiddle_visibility(plot_params)
end
















function save_routine(varargin) % additional arguments go afer the first two, which must be put here ,even if they are not used.
%perhaps get info from other calls here?

global plot_params
f = fieldnames(plot_params.panels);

for ii = 1:length(f)
	set( findall(plot_params.panels.(f{ii}), '-property', 'visible'), 'visible', 'off')
	set(plot_params.panels.(f{ii}),'visible','off');
	
end


if or(plot_params.switches.display_faces == 1,plot_params.switches.display_face_samples == 1)
    plot_params.format = 'png';
    plot_params.format_flag = 'png';
else
    plot_params.format = 'eps';
    plot_params.format_flag = 'epsc2';
end


render_into_file(plot_params,'-r150');

for ii = 1:length(f)
	set(plot_params.panels.(f{ii}),'visible','on');
	set( findall(plot_params.panels.(f{ii}), '-property', 'visible'), 'visible', 'on')
end

end









function render_legends()
global plot_params
f = fieldnames(plot_params.legend);

a = [];
b = [];


for ii = 1:length(f)
	a = [a plot_params.legend.(f{ii}).handles]; 
	b = [b plot_params.legend.(f{ii}).text];
end



plot_params.handles.legend = legend(a,b); 
set(plot_params.handles.legend,'Location','NorthEastOutside','Interpreter','none');
%,'Location','SouthEast','Parent',plot_params.axes.vertices
end






















function fv_to_stl(varargin)
global plot_params

fv = plot_params.faces_and_vertices;
save('fv.mat','fv');



%delete the bad degenerate faces.  
degen = any(diff(fv.faces(:,[1:3 1]),[],2)==0,2);
fv.faces(degen,:)
  fv.faces(degen,:) = [];

  fv.faces(1,:) = fv.faces(1,[3 2 1]);
  
  % Fix non-uniform face orientations
  fv2 = unifyMeshNormals(fv,'alignTo',1);
  % Solidify
  thickened_fv = surf2solid(fv2,'thickness',0.09);
  
  
  
stlwrite('br_surf.stl',fv2);




stlwrite('br_thickened_surf.stl',thickened_fv);


figure, patch(thickened_fv,'FaceColor','g','FaceAlpha',0.5,'EdgeAlpha',0.5), axis image

end



function center_camera_on_selected_point(source, event)
global plot_params

curr_axes = plot_params.axes.main;

dcm_obj = datacursormode(gcf);
set(dcm_obj,'DisplayStyle','datatip',...
'SnapToDataVertex','on','Enable','on');

set(curr_axes,'CameraPositionMode','manual');

curr_pos = get(curr_axes,'CameraPosition');


cursor_info = getCursorInfo(dcm_obj);


if isfield(cursor_info,'Position')
	

	plot_params.scene.target = cursor_info.Position;
	set(curr_axes,'CameraTarget',cursor_info.Position);
    set(curr_axes,'CameraPosition',curr_pos);
else
	display('select a point');
end






	
end

%view functions

function update_sliders(src, evt,az_handle, el_handle)

global plot_params
[plot_params.azimuth, plot_params.elevation] = view;
set(az_handle,'Value',plot_params.azimuth)
set(el_handle,'Value',plot_params.elevation)

update_view_text;
end


function elevation_changer(src, evt)
global plot_params
plot_params.elevation = get(src,'value');
change_current_view;
end


function azimuth_changer(src, evt)
global plot_params
plot_params.azimuth = get(src,'value');
change_current_view;
end


function change_current_view
global plot_params
view(plot_params.azimuth, plot_params.elevation);
update_view_text;
end


function update_view_text()
global plot_params
set(plot_params.controls.viewdisplaytext,'string',['az ' num2str(plot_params.azimuth) ' el ' num2str(plot_params.elevation)]);
end

