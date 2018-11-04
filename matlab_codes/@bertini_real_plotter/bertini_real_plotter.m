%
% a class for plotting the data from a bertini_real run of any computable
% dimension
%
% example invokation: 
% bertini_real_plotter('autosave',false,'vertices',false,'linestyle','-','colormap',@jet,'colorfn',@norm,'num_colors',200,'curves',false)
%
%
% options: 
%	'autosave'          - bool [false]
%	'vertices', 'vert'  - bool [false]  for large samplings, this causes
%							rendering to be very slow.
%	'filename'          - string [BRinfo*.mat]
%	'file'              - struct.  no default. the result of reading
%	                        gathered data
%	'proj'              - handle to function.  no default
%	'mono', 'monocolor' - color or RGB triple.  not on by default
%	'labels'            - bool [true]
%  	'colormap'          - handle to colormap generating function [@jet]
%  	'linestyle'         - string, sets all curves to have this style. 
%                            no default, curves have style by type.
%	'curves', 'curve'   - bool [true]
%	'faces'             - bool [true]
%
%	'indices'		- array [empty]  plot these variables only.  if you use 
%							more than 3, you will generate errors.
%							this accomplishes coordinate projection onto these 
%							coordinates
%							IMPORTANT NOTE:  the indices in this array
%							should be 1-based, not 0-based.  
%
%	'whichfaces'		- array [empty]  If you are plotting faces, and
%							this is empty, all faces will be plotted.  
%							if this is non-empty, only those faces with
%							indices in this array will be rendered.
%							IMPORTANT NOTE:  the indices in this array
%							should be 1-based, not 0-based.  
%	'touchingedgesonly'	- bool [true] Only render edges of curves touching 
%							faces rendered.  Conditional on object being a
%							surface.  For native curves, this will be false
%							because there are no faces to touch
%
%   'colorfn'           - handle to function of x, for generating color
%                            data.  no default value.  if this is not
%                            specified, then the colors for the faces
%                            correspond to the entire face, in the order
%                            computed.  An example of using this colorfn
%                            would be to pass a handle to a function
%                            computing the distance between x and the
%                            origin, perhaps.
%    'colorfnprojdata'       - bool [false]  If true, then the data from
%                            the projection function is used for the color 
%                            function you passed in.  If false, the plain
%                            old coordinates from space will be used for
%                            color, even if you are passing them through a
%                            projection to plot them.
%
%    'num_colors'        - integer [64] the number of colors used in the
%							colormap, particularly used when you use the 
%							'colorfn' option, to
%							specify a function used for coloring the surface.
%							also, 'numcolors'
%
%    'view_through_pi'   - bool [false] use projection through that which 
%                            was used for decomposition as space
%                            coordinates.  can make some very pretty plots.
%                            
%
%
%
%
% danielle brake
% danielleamethystbrake@gmail.com
%
% university of wisconsin -- eau claire
% mathematics
% spring, summer 2018
%
% university of notre dame
% applied and computational mathematics and statistics
% 2014, 2015, 2016, 2017



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
		data = [];
	end
	
	methods
		
		function br_plotter = bertini_real_plotter(varargin)
			initialize(br_plotter);
			set_default_options(br_plotter);
			set_options(br_plotter,varargin);
			
			if needs_data(br_plotter)
				set_filename(br_plotter);
				load_and_process(br_plotter);
			else
				just_process(br_plotter);
			end

			plot(br_plotter);
		end %re: bertini_real_plotter() constructor
		
		function initialize(br_plotter)
			path_sanity_check()
			initialize_handles_surface(br_plotter);
		end
		
		%why is there not a symmetric one for the curve case? -20180614
		initialize_handles_surface(br_plotter);
			
			
		function does_it = needs_data(br_plotter)
			if isempty(br_plotter.BRinfo)
				does_it = true;
			else
				does_it = false;
			end
		end
		
		
		function load_and_process(br_plotter)
			load_data(br_plotter);
			just_process(br_plotter);
		end
		
		
		function just_process(br_plotter)
			sanity_checks(br_plotter)
			adjust_version_numbers(br_plotter)
			version_check(br_plotter);
			set_options_from_BRinfo(br_plotter);
			preprocess_data(br_plotter);
		end
			
		
		function sanity_checks(br_plotter)
			if br_plotter.BRinfo.num_vertices==0
				warning('your decomposition contains 0 vertices.  The real part appears to be empty.  Plotting will now terminate.')
				display(br_plotter.BRinfo);
				return;
			end
		end
		
		
		
		function adjust_version_numbers(br_plotter)
			if ~isfield(br_plotter.BRinfo,'run_metadata')
				br_plotter.BRinfo.run_metadata.version.number = 103; %never change this number.  realistically it should be a commit number for the repo, but... i'm lazy
			end
		end
		
		
		
		function set_options_from_BRinfo(br_plotter)
			
			get_indices(br_plotter);

			
			[br_plotter.options.containing, br_plotter.options.basename, ~] = fileparts(pwd);
			br_plotter.dimension = br_plotter.BRinfo.dimension;

			if (br_plotter.BRinfo.dimension == 2)
				if isempty(br_plotter.options.which_faces)
					br_plotter.options.which_faces = 1:br_plotter.BRinfo.num_faces;
				elseif max(br_plotter.options.which_faces) > br_plotter.BRinfo.num_faces
					error('trying to plot faces which don''t exist. requested index %i > num faces %i',...
						max(br_plotter.options.which_faces),...
						br_plotter.BRinfo.num_faces);
				end
			end

			if (br_plotter.BRinfo.dimension == 1)
				br_plotter.options.touching_edges_only = false;
			elseif (br_plotter.BRinfo.dimension == 2)
				if length(br_plotter.options.which_faces)==br_plotter.BRinfo.num_faces
					br_plotter.options.touching_edges_only = false;
				end
			end
		end
		
		
		
		
		function version_check(br_plotter)
			if isfield(br_plotter.BRinfo.run_metadata.version, 'gather')
				if br_plotter.BRinfo.run_metadata.version.gather < 150
					error('this version of bertini_real_plotter requires data gathered with gather_br_samples at least 150.  please re-gather');
				end
			else
				error('this version of bertini_real_plotter requires data gathered with gather_br_samples at least 150.  please re-gather');
			end
		end
		
		
		
		
		set_default_options(br_plotter)

		%parses the command line options fed into the constructor.
		set_options(br_plotter,command_line_options)
		
		% uses internally set variable 'filename' to load a .mat file
		% containing data gathered previously.
		load_data(br_plotter)
		
		
		function set_filename(br_plotter)
			
			if isempty(br_plotter.filename)
				prev_filenames = dir('BRinfo*.mat');
				
				if isempty(prev_filenames)
					br_plotter.filename = uigetfile();
				else
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
			end
			
		end
		

		
		function load_and_render(br_plotter,source, event)
			
			clear_for_load(br_plotter)
			
			[FileName,PathName,FilterIndex] = uigetfile();
			br_plotter.filename = [PathName FileName];
			
			load_and_process(br_plotter);
			plot(br_plotter);
			
		end
		
		function clear_for_load(br_plotter)
			
			
		end
		
		function plot(br_plotter)
			
			if br_plotter.BRinfo.num_vertices==0
				warning('your decomposition contains 0 vertices.  The real part appears to be empty.')
				return;
			end
				
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
		
		function hide_panels(br_plotter,varargin)
			f = fieldnames(br_plotter.panels);

			for ii = 1:length(f)
				set( findall(br_plotter.panels.(f{ii}), '-property', 'visible'), 'visible', 'off')
				set(br_plotter.panels.(f{ii}),'visible','off');

			end
		end
		
		
		function show_panels(br_plotter,varargin)
			f = fieldnames(br_plotter.panels);

			for ii = 1:length(f)
				set( findall(br_plotter.panels.(f{ii}), '-property', 'visible'), 'visible', 'on')
				set(br_plotter.panels.(f{ii}),'visible','on');
			end
			
			hide_restore_panel(br_plotter)
			cameratoolbar('show')
		end
		
		
		function minimize_panels(br_plotter,varargin)
			f = fieldnames(br_plotter.panels);

			for ii = 1:length(f)
				set( findall(br_plotter.panels.(f{ii}), '-property', 'visible'), 'visible', 'off')
				set(br_plotter.panels.(f{ii}),'visible','off');
			end
			cameratoolbar('hide')
			br_plotter.show_restore_panel()
		end
		
		
		function show_restore_panel(br_plotter,varargin)
			set( findall(br_plotter.panels.restore, '-property', 'visible'), 'visible', 'on')
			set(br_plotter.panels.restore,'visible','on');
		end
		
		
		function hide_restore_panel(br_plotter,varargin)
			set( findall(br_plotter.panels.restore, '-property', 'visible'), 'visible', 'off')
			set(br_plotter.panels.restore,'visible','off');
		end
		
		
		% declared headers for functions in other files.
		
		preprocess_data(br_plotter)

		get_indices(br_plotter)
		
		val = edge_touches_faces(br_plotter, edge_index, curve)
		
		
		%functions specific to surfaces
		surface_plot(br_plotter)
		handles = plot_surface_edges(br_plotter)
		
		
		%functions specific to curves
		curve_plot(br_plotter)
		handles = plot_curve_samples(br_plotter,sampler_data,style, color)
		
		
		% common to all dimensions
		sphere_plot(br_plotter)
		plot_vertices(br_plotter)
		plot_edge_points(br_plotter)
		
		setupfig(br_plotter,varargin)
		
		create_axes(br_plotter)
		label_axes(br_plotter)
		sync_axes(br_plotter)
		
		render_legends(br_plotter)
		visibility_setup(br_plotter)
		
		
		
		% calls the initial_visibility, visibility_setup, and
		% twiddle_visibility functions
		controls(br_plotter)
		camera_setup(br_plotter)
		
		
		
		
		% setup the interactive buttons
		button_setup(br_plotter)
		set_initial_visibility(br_plotter)
		twiddle_visibility(br_plotter)
		
		
		
		
		% callback functions
		
		% for more info on associating callbacks with buttons, see e.g.
		% http://www.mathworks.com/help/matlab/matlab_oop/class-methods-for-graphics-callbacks.html
		change_alpha(br_plotter,source,event) % is a callback function
		change_edge_alpha(br_plotter,source,event) % is a callback function
		change_edge_width(br_plotter,source,event) % is a callback function
        change_line_width(br_plotter,source,event) % is a callback function
		change_text_size(br_plotter,source,event)% is a callback function
		center_camera_on_selected_point(br_plotter,source, event)
		save_routine(br_plotter,varargin) % is a callback function
		flip_switch(br_plotter,srcHandle,eventData,varargin)
		resizeui(br_plotter,srcHandle,eventData,varargin)
		
		
	
	end%re: methods
	
	methods(Static)
		y = project_through_pi(x,pi)
	end%re: static methods

	
end
