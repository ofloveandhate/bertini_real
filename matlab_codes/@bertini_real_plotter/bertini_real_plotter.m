%
% a class for plotting the data from a bertini_real run of any computable
% dimension
%
%
% options: 
%	'autosave'          - bool [true]
%	'vertices', 'vert'  - bool [true]
%	'filename', 'file'   - string [BRinfo*.mat]
%	'proj'              - handle to function.  no default
%	'mono', 'monocolor' - color or RGB triple.  not on by default
%	'labels'            - bool [true]
%  	'colormap'          - handle to colormap generating function [@jet]
%  	'linestyle'         - string, sets all curves to have this style. 
%                            no default, curves have style by type.
%	'curves', 'curve'   - bool [true]
%	'faces'             - bool [true]
%
%
%
%


% daniel brake
% danielthebrake@gmail.com
% university of notre dame
% applied and computational mathematics and statistics
% 2014, 2015



classdef bertini_real_plotter < handle
	
	
	properties
		BRinfo = [];
% 		scene = scene_manipulator();
		
		window = [20 20 1024 768];
		figures = []
		axes = [];
		handles = [];
		legend = [];
		filename = [];
		
		panels = [];
		buttons = [];
		switches = [];
		checkboxes = [];
		
		scene = [];
		
		dimension = -1;
		
		indices = [];
		
		options = [];
		
		is_bounded = [];
		fv = [];
	end
	
	methods
		
		
		function br_plotter = bertini_real_plotter(varargin)
			
			set_default_options(br_plotter);
			
			set_options(br_plotter,varargin);
			
			
			load_data(br_plotter);
			
			if br_plotter.options.use_custom_projection
				preprocess_data(br_plotter);
			end
			
			get_indices(br_plotter);
			
			plot(br_plotter);
		end %re: bertini_real_plotter() constructor
		
		
		
		function set_default_options(br_plotter)
			br_plotter.options.use_custom_projection = false;
			br_plotter.options.markersize = 10;
			br_plotter.options.sample_alpha = 1;
			br_plotter.options.face_alpha = 1;
			br_plotter.options.edge_alpha = 0.4;
			br_plotter.options.fontsizes.legend = 12;
			br_plotter.options.fontsizes.labels = 16;
			br_plotter.options.fontsizes.axis = 20;
			br_plotter.options.line_thickness = 6;
			br_plotter.options.autosave = true;
			
			br_plotter.options.labels = true;
			br_plotter.options.monocolor = false;
			
			br_plotter.options.render_vertices = true;
			br_plotter.options.render_curves = true;
			br_plotter.options.render_faces = true;
			
            if isempty(which('parula'))
                br_plotter.options.colormap = @jet;
            else
                br_plotter.options.colormap = @parula;
            end
		end
		
		
		function set_options(br_plotter,command_line_options)
			
			
			if mod(length(command_line_options),2)~=0
				error('must have option-value pairs');
			end
			
			
			for ii = 1:2:length(command_line_options)-1
				
				switch command_line_options{ii}
					case 'autosave'
						
						tentative_arg = command_line_options{ii+1};
						
						if ischar(tentative_arg)
							switch tentative_arg
								case {'y','yes','true'}
									br_plotter.options.autosave = true;
								case {'n','no','false'}
									br_plotter.options.autosave = false;
								otherwise
									error('bad option %s for autosave',tentative_arg);
							end
							
						else
							if tentative_arg==1
								br_plotter.options.autosave = true;
							elseif tentative_arg==0
								br_plotter.options.autosave = false;
							else
								error('bad option %f for autosave',tentative_arg);
							end
						end
					
					case {'curves','curve'}
						tentative_arg = command_line_options{ii+1};
						
						if ischar(tentative_arg)
							switch tentative_arg
								case {'y','yes','true'}
									br_plotter.options.render_curves = true;
								case {'n','no','none','false'}
									br_plotter.options.render_curves = false;
								otherwise
									error('bad option %s for curves',tentative_arg);
							end
							
						else
							if tentative_arg==1
								br_plotter.options.render_curves = true;
							elseif tentative_arg==0
								br_plotter.options.render_curves = false;
							else
								error('bad option %f for curves',tentative_arg);
							end
						end
						
						
					case {'vertices','vert'}
						
						tentative_arg = command_line_options{ii+1};
						
						if ischar(tentative_arg)
							switch tentative_arg
								case {'y','yes','true'}
									br_plotter.options.render_vertices = true;
								case {'n','no','none','false'}
									br_plotter.options.render_vertices = false;
								otherwise
									error('bad option %s for vertices',tentative_arg);
							end
							
						else
							if tentative_arg==1
								br_plotter.options.render_vertices = true;
							elseif tentative_arg==0
								br_plotter.options.render_vertices = false;
							else
								error('bad option %f for vertices',tentative_arg);
							end
						end
						
						
					case {'filename','file'}
						br_plotter.filename = command_line_options{ii+1};
						if ~ischar(br_plotter.filename)
							error('filename argument must be a filename')
						end
					case 'proj'


						tmp = command_line_options{ii+1};

						if isa(tmp,'function_handle')
							br_plotter.options.custom_projection = tmp;
							br_plotter.options.use_custom_projection = true;
						elseif strcmpi(tmp,'natural')
	
						else
							error('value for ''proj'' must be a function handle or ''natural''');
						end

					case {'colormap'}
						tmp = command_line_options{ii+1};
						if isa(tmp,'function_handle')
							br_plotter.options.colormap = tmp;
						else
							error('value for ''colormap'' must be a handle to a function generating a colormap for an integer number of colors; e.g. @jet');
						end
						
						
					case {'mono','monocolor'}
						br_plotter.options.monocolor = true;
						
						tentative_color = command_line_options{ii+1};
						
						
						
						if ischar(tentative_color)
							switch tentative_color
								case 'r'
									br_plotter.options.monocolor_color = [1 0 0];
								case 'g'
									br_plotter.options.monocolor_color = [0 1 0];
								case 'b'
									br_plotter.options.monocolor_color = [0 0 1];
								case 'm'
									br_plotter.options.monocolor_color = [1 0 1];
								case 'c'
									br_plotter.options.monocolor_color = [0 1 1];
								case 'y'
									br_plotter.options.monocolor_color = [1 1 0];
								case 'k'
									br_plotter.options.monocolor_color = [0 0 0];
								
								otherwise
									error('input color string must be one of r g b m c y k.  you can also specify a 1x3 RGB color vector');
							end
						else
							[m,n] = size(tentative_color);
							if and(m==1,n==3)
								br_plotter.options.monocolor_color = tentative_color;
							else
								error('explicit color for monocolor surfaces must be a 1x3 RGB vector');
							end
							
						end
						
						
						
						br_plotter.options.colormap = @(num_colors) repmat(br_plotter.options.monocolor_color,num_colors,1);
						
					case 'labels'
						tentative_arg = command_line_options{ii+1};
						
						if ischar(tentative_arg)
							switch tentative_arg
								case {'y','yes','true'}
									br_plotter.options.labels = true;
								case {'n','no','none','false'}
									br_plotter.options.labels = false;
								otherwise
									error('bad option %s for labels',tentative_arg);
							end
							
						else
							if tentative_arg==1
								br_plotter.options.labels = true;
							elseif tentative_arg==0
								br_plotter.options.labels = false;
							else
								error('bad option %f for labels',tentative_arg);
							end
						end
						
					case 'linestyle'
						br_plotter.options.linestyle = command_line_options{ii+1};
						br_plotter.options.use_fixed_linestyle = true;
						
						
					case 'faces'
						
						tentative_arg = command_line_options{ii+1};
						
						if ischar(tentative_arg)
							switch tentative_arg
								case {'y','yes','true'}
									br_plotter.options.render_faces = true;
								case {'n','no','none','false'}
									br_plotter.options.render_faces = false;
								otherwise
									error('bad option %s for faces',tentative_arg);
							end
							
						else
							if tentative_arg==1
								br_plotter.options.render_faces = true;
							elseif tentative_arg==0
								br_plotter.options.render_faces = false;
							else
								error('bad option %f for faces',tentative_arg);
							end
						end
						
												
					otherwise
						error('unexpected option name ''%s''',command_line_options{ii})
				end
			end
		end
		
		
		function load_data(br_plotter)
			
			if isempty(br_plotter.filename)
				prev_filenames = dir('BRinfo*.mat');
				
				if isempty(prev_filenames)
					error('no obvious BRinfo files to load');
				end
				
				max_found = -1;

				for ii = 1:length(prev_filenames)
					curr_name = prev_filenames(ii).name;
					curr_num = str2double(curr_name(7:end-4));
					if max_found < curr_num
						max_found = curr_num;
					end

				end
				br_plotter.filename = ['BRinfo' num2str(max_found) '.mat'];
			end
			
			
			
			if isempty(dir(br_plotter.filename))
				error('nexists file with name ''%s''',br_plotter.filename);
			end
			
			
			file_variables = whos('-file',br_plotter.filename);
			
			if ismember('BRinfo', {file_variables.name})
				temp = load(br_plotter.filename);
				br_plotter.BRinfo = temp.BRinfo;
			else
				error('file ''%s'' does not contain variable ''BRinfo''',br_plotter.filename);
			end
			
			
			
			
			

			[br_plotter.options.containing, br_plotter.options.basename, ~] = fileparts(pwd);

			br_plotter.dimension = br_plotter.BRinfo.dimension;

		end
		
		
		function plot(br_plotter)
			
			
			setupfig(br_plotter);
			
			
			switch br_plotter.dimension
				case 1
					curve_plot(br_plotter);

				case 2	
					surface_plot(br_plotter);

				otherwise

			end

			if ~isempty(br_plotter.legend)
				render_legends(br_plotter);
			end
			
			
			br_plotter.options.plotted = 1;
			
			
			
			button_setup(br_plotter);
			
			controls(br_plotter);
			
			if br_plotter.options.autosave
				try
					save_routine(br_plotter);
				catch exception
					display(exception);
					display('saving render not completed');
				end
			end
			
			
		end
		
		
		function set_label_text_size(br_plotter,~,~,new_size)
			
			f = fieldnames(br_plotter.handles.vertex_text);
			
			for ii = 1:length(f)
				set(br_plotter.handles.vertex_text.(f{ii}),'FontSize',new_size);
			end
			
			
			switch br_plotter.dimension
				case 1
					
					
				case 2
					all_edge_labels = [br_plotter.handles.critcurve_labels;...
								br_plotter.handles.spherecurve_labels;...
								br_plotter.handles.midtext;...
								br_plotter.handles.crittext;...
								br_plotter.handles.singtext];
					
					set(all_edge_labels,'FontSize',new_size);
				otherwise
					
					
			end
			
			br_plotter.options.fontsizes.labels = new_size;
		end
	
		function set_legend_text_size(br_plotter,~,~,new_size)
			set(br_plotter.handles.legend,'FontSize',new_size);
			br_plotter.options.fontsizes.legend = new_size;
		end
		
		function set_axis_text_size(br_plotter,~,~,new_size)
			set(br_plotter.axes.main,'FontSize',new_size);
			br_plotter.options.fontsizes.axis = new_size;
		end
		
		
		function delete_panels(br_plotter)
			
			b = fieldnames(br_plotter.panels);
			for ii = 1:length(b)
				delete(br_plotter.panels.(b{ii}));
			end
			br_plotter.panels = [];
			
		end
		
		
		
		
		
		
		
		% declared headers for functions in other files.
		
		preprocess_data(br_plotter)

		get_indices(br_plotter)
		
		setupfig(br_plotter,varargin)
		
		create_axes(br_plotter)
		label_axes(br_plotter)
		sync_axes(br_plotter)
		
		visibility_setup(br_plotter)
		
		surface_plot(br_plotter)
		handles = plot_surface_edges(br_plotter)
		
		
		curve_plot(br_plotter)
		
		handles = plot_curve_samples(br_plotter,sampler_data,style, color)
		
		% common to all dimensions
		sphere_plot(br_plotter)
		plot_vertices(br_plotter)
		
		plot_edge_points(br_plotter)
		
		
		
		render_legends(br_plotter)
		
		
		
		
		% calls the initial_visibility, visibility_setup, and
		% twiddle_visibility functions
		controls(br_plotter)
		
		camera_setup(br_plotter)
		
		center_camera_on_selected_point(br_plotter,source, event)
		
		save_routine(br_plotter,varargin) % is a callback function
		
		
		% setup the interactive buttons
		button_setup(br_plotter)
		
		
		
		%http://www.mathworks.com/help/matlab/matlab_oop/class-methods-for-graphics-callbacks.html
		
		
		change_alpha(br_plotter,source,event) % is a callback function
		change_text_size(br_plotter,source,event)% is a callback function

		set_initial_visibility(br_plotter)
		twiddle_visibility(br_plotter)
		
		flip_switch(br_plotter,srcHandle,eventData,varargin)
		
		
		
		function resizeui(br_plotter,srcHandle,eventData,varargin)
			br_plotter.window = get(br_plotter.figures.main,'Position');
			w = br_plotter.window;
			if ~isempty(br_plotter.panels)
				
				p = get(br_plotter.panels.buttons,'Position');
				set(br_plotter.panels.buttons,'position',[5    w(4)-p(4)-5     p(3)    p(4)]);
				
				
				c = get(br_plotter.panels.common_visibility,'Position');
				total_vertical_size = c(4);
				if isfield(br_plotter.panels,'vertex')
					v = get(br_plotter.panels.vertex,'Position');
					total_vertical_size = total_vertical_size+v(4);
				end
				
				if isfield(br_plotter.panels,'surface')
					s = get(br_plotter.panels.surface,'Position');
					total_vertical_size = total_vertical_size+s(4);
				end
				
				if isfield(br_plotter.panels,'curve')
					cu = get(br_plotter.panels.curve,'Position');
					total_vertical_size = total_vertical_size+cu(4);
				end
				

				if total_vertical_size+30 > w(4)-p(4)
					P = w(3)-c(3)-5;
				else
					P = 5;
				end
				
				set(br_plotter.panels.common_visibility,'position',[P 5 c(3)    c(4)]);
				Q = 10+c(4);
				if isfield(br_plotter.panels,'vertex')
					set(br_plotter.panels.vertex,'position',[P Q v(3)    v(4)]);
					Q = Q+5+v(4);
				end

				if isfield(br_plotter.panels,'surface')
					set(br_plotter.panels.surface,'position',[P Q s(3)    s(4)]);
					Q = Q+5+s(4);
				end

				if isfield(br_plotter.panels,'curve')
					set(br_plotter.panels.curve,'position',[P Q cu(3)    cu(4)]);
					Q = Q+5+cu(4);
				end				
				
			end
		end	
				
	end%re: methods
	
	

	
end
