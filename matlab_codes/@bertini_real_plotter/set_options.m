%parses the command line options fed into the constructor.
function set_options(br_plotter,command_line_options)
			
			
if mod(length(command_line_options),2)~=0
	error('must have option-value pairs');
end


for ii = 1:2:length(command_line_options)-1

	switch lower(command_line_options{ii})
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
						error('bad option %s for ''vertices''',tentative_arg);
				end

			else
				if tentative_arg==1
					br_plotter.options.render_vertices = true;
				elseif tentative_arg==0
					br_plotter.options.render_vertices = false;
				else
					error('bad option %f for ''vertices''',tentative_arg);
				end
			end


		case {'filename'}
			if ~ischar(br_plotter.filename)
				error('''filename'' argument must be a string filename')
			end
			
			br_plotter.filename = command_line_options{ii+1};
			

		case {'file'}
			if ~isstruct(command_line_options{ii+1})
				error('''file'' needs to be a BRinfo struct.');
			end
			
			br_plotter.BRinfo = command_line_options{ii+1};
			
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

		case {'colorfn'}
			tmp = command_line_options{ii+1};
			if isa(tmp,'function_handle')
				br_plotter.options.use_colorfn = true;
				br_plotter.options.colorfn = tmp;
			else
				error('value for ''colorfn'' must be a handle to a function accepting a real matrix in which the rows are points, and returning a real vector');
			end

		case {'colorfn_uses_raw'}
			tentative_arg = command_line_options{ii+1};
			if ischar(tentative_arg)
				switch tentative_arg
					case {'y','yes','true'}
						br_plotter.options.colorfn_uses_raw = true;
					case {'n','no','none','false'}
						br_plotter.options.colorfn_uses_raw = false;
					otherwise
						error('bad option %s for ''colorfn_uses_raw''',tentative_arg);
				end

			else
				if tentative_arg==1
					br_plotter.options.colorfn_uses_raw = true;
				elseif tentative_arg==0
					br_plotter.options.colorfn_uses_raw = false;
				else
					error('bad option %f for ''colorfn_uses_raw''',tentative_arg);
				end
			end


		case {'num_colors','numcolors'}
			tmp = command_line_options{ii+1};
			if ~isint(tmp)
				error('value for ''num_colors'' must be in integer');
			end

			br_plotter.options.num_colors = tmp;

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
						error('input color string for mono must be one of r g b m c y k.  you can also specify a 1x3 RGB color vector');
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

		case 'linewidth'
			br_plotter.options.linewidth = command_line_options{ii+1};

		case 'edgecolor'
			br_plotter.options.linewidth = command_line_options{ii+1};
			
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

		case 'whichfaces'

			tentative_arg = command_line_options{ii+1};
			if ~isnumeric(tentative_arg)
				error('argument for ''whichfaces'' must be integer array');
			end

			br_plotter.options.which_faces = tentative_arg;

		case 'touchingedgesonly'
			tentative_arg = command_line_options{ii+1};
			br_plotter.options.touching_edges_only = tentative_arg;

		otherwise
			error('unexpected option name ''%s''',command_line_options{ii})
	end
end
end