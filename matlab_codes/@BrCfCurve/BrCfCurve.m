% a class for interfacing between bertini_real curves, and chebfun
% 
% ---
%
% the intended functionality is to provide you with a way to optimize over
% algebraic curves in arbitrary dimensions.
%
% this class provides overloads for functions like min, max, roots, plot, etc
%
% ---
%
% how to use: 
% 1. do a bertini_real decomposition of a curve
% 2. in matlab, navigate to that folder, and `gather_br_samples`
% 3. make an empty BrCfCurve object.  `C = BrCfCurve()`
% 4. compute the chebfuns.  `C.ComputeChebfuns()`
% 5. do what you want with it.  perhaps
% 6.a plot a function: `plot(C,'func',f)` where `f` is a handle to an n-ary
%       function
% 6.b find the mins or maxes (or heck, both) of a function: `C.max(f)` where `f` is a handle to an n-ary
%       vectored function (it should operate on columns of data, as those are the variables during evaluation)
% 
% ---
%
% this class has code dependencies:  
% 1. the code in github.com/ofloveandhate/brakelab/bertini1, it *must*
%     be on your matlab path
% 2. a plotting function called `color_line`.  how about the one at https://www.mathworks.com/matlabcentral/fileexchange/19476-colored-line-or-scatter-plot
%
% ---
%
% the following chebfun arguments are passed through to chebfun:
% 'fixedlength'
% 'minsamples'
% 'maxlength'
% 'sampletest'
% 'splitting'
% 
% see http://www.chebfun.org/docs/guide/guide08.html for the official
% documentation
% 
% ---
% 
% we have provided the following arguments for the constructor of this class:
% 
% * 'indices' - restrict the computation of chebfuns for these variables only.
% * 'tol' - for min/max/minandmax, how close to the boundary or midpoint of and edge a value
% must be to consider for discarding due to spuriousness.  
% this is in terms of the t \in [0,1] parameterization using chebfuns, not a 
% euclidean space distance in the variables.  1e-14 by default
% * `echo` - whether to duplicate log contents to screen.  false by default
% * `filename` - the filename for the .mat file for the gathered BRinfo file
% that contains the decomposition of the curve
% * `autoload` - whether to automatically load the data from the BRinfo.mat file
% when you create the BrCfCurve.  default: true
% * `autosetup` - whether you want the BrCfCurve to automatically set up things
% to do some tracking.  I have no idea why you would want to not do this, except
% that you want to change into a temporary directory or something before
% starting the process.  
%
% ---
%
% the computed chebfuns are such that each edge is broken into two pieces, at
% the midpoint of the edge.  on each half-edge, a chebfun is constructed for eac
% variable.  the interval used is t \in [0,1], with t=0 being the midpoint of the
% edge, and t=1 being the boundary.
%
% ---
%
% important notes:
%
% this class is a handle class.  This means that "copying" it does not
% deep-copy, but you get a shallow reference-like copy.  Beware.
%
% ---
% 
% this code contains known flaws:
% * it does not compute mins/maxes/roots on isolated points of a curve 
% in the optimization routines.  if you need to be able to do this, 
% contact the author, Silviana.  She'd be happy to hear from you.
%
% ---
%
%
% this code comes with absolutely no warranty for any purpose whatsoever.  it
% would be amazing if it worked at all for another human aside from the author.  
% if it does, or doesn't, drop her a line will ya?
%
% silviana amethyst
% university of notre dame, university of wisconsin -- eau claire
% fall 2016 - fall 2018 - fall 2021
% amethyst@uwec.edu
% 


classdef BrCfCurve < handle


	properties
		br_filename;
		BRinfo = [];  % the output from Bertini_real
		b_input = bertini_input();  %parse dat input file, yo.  from bertini_tropical
		options = [];
		
		logfile = [];
		
		edge_chebfuns = [];
		have_chebfuns = false;
		
		func = [];
		cache = struct();
		
		edge_computed_points; % a map of computed points, so can look up rather than compute fresh
		current_working_edge = 0; % default is impossible value, since matlab uses 1-based indexing
	end % re: properties
	
	methods (Access = public)
		
		% constructor for BrCfCurves.
		% 
		% loads the .mat file, produced from gather_br_samples
		function this = BrCfCurve(varargin)

			this.NewLogfile();

			SetDefaultOptions(this);
			
			this.Log('creating BrCfCurve object');
			
			SetOptions(this,varargin);

			if this.options.autoload
				if isempty(this.br_filename)
					this.SetDefaultFilename();
				end
				LoadData(this);
				
			end
			
			if this.options.autosetup
				this.InitialSetup();
			end
			
		end % re: constructor
			
		

		
		
		
		function BertiniSetup(this)
			this.b_input = bertini_input(); % clear the input file.  it's a handle, so persists in unexpected places
			this.SetupInput();
			this.SetupSliceSystem();
		end
		
		function NewLogfile(this)
			this.logfile.name = sprintf('BrCfLogfile_%s.txt',datetime('now','TimeZone','local','Format','yyyy.MM.dd.HH:mm:ss'));
			this.logfile.fid = fopen(this.logfile.name,'w');
		end
		
		
		
		
		
		
		%looks for the highest-numbered BRinfo.mat file in the current
		%folder, and sets that to be the current filename.
		function SetDefaultFilename(this)
			prev_filenames = dir('BRinfo*.mat');

			if isempty(prev_filenames)
				error('no obvious BRinfo files to load, I think you need to `gather_br_samples`');
			end

			max_found = -1;

			for ii = 1:length(prev_filenames)
				curr_name = prev_filenames(ii).name;
				curr_num = str2double(curr_name(7:end-4));
				if max_found < curr_num
					max_found = curr_num;
				end

			end
			this.br_filename = ['BRinfo' num2str(max_found) '.mat'];
		end
		
		
		
		
		
		
		
		
		% sets the current function, which is used for `min`, `max`,
		% `minandmax`, `plot`, and `roots`
		%
		% the reason this exists is because this class uses memoization to help
		% reduce the amount of computation needed, by cacheing the function
		% chebfuns and re-using them when possible
		function set_function(this, func)
			this.ValidateHandle(func);
			this.func = func;
			this.cache.have_functions = false;
		end
		
		% clears the current function from the BrCfCurve.  
		function clear_function(this)
			this.func = [];
			this.cache.have_functions = false;
		end
		
		% errors if don't have a function, or the function is incompatible
		function require_have_function(this)
			if isempty(this.func)
				error('function empty in BrCfCurve.  set the function with `set_function`, using an n-ary function handle');
			end
			
			this.ValidateHandle(this.func);
		end
		
		
		function set_indices(this, ind)
			this.options.plot.indices = ind;
		end
		
		% load the BRinfo file, with name stored at this.br_filename
		% as a data member of the BrCfCurve
		% 
		% will cause errors if you cause a surface to be loaded, e.g.
		function LoadData(this)
			
			if isempty(dir(this.br_filename))
				error('nexists file with name ''%s''',this.br_filename);
			end
			
			
			file_variables = whos('-file',this.br_filename);
			
			if ismember('BRinfo', {file_variables.name})
				temp = load(this.br_filename);
				this.ValidateBRinfo(temp.BRinfo);
				this.BRinfo = temp.BRinfo;
			else
				error('file ''%s'' does not contain variable ''BRinfo''',this.br_filename);
			end

			[this.options.containing, this.options.basename, ~] = fileparts(pwd);
	
			this.Log('loaded bertinireal matlab data file %s\n', [this.options.containing '/' this.options.basename '/' this.br_filename]);
			
			
			
			if isnan(this.options.plot.indices)
				this.options.plot.indices = 1:this.BRinfo.num_variables-1;
			end
		end % re: load data
		
		
		
		
		
		
		%
		% this -- the curve over which to compute the max
		%
		% the function is set via C.set_function(func), not passed into `max`
		% 
		%
		% extra arguments are passed along to Chebfun
		function [varargout] = max(this,varargin)
			[varargout{1:nargout}] = min_or_max(this,'max',varargin{:});
		end
		
		%
		% this -- the curve over which to compute the min
		%
		% the function is set via C.set_function(func), not passed into `min`
		% 
		% extra arguments are passed along to Chebfun
		function [varargout] = min(this,varargin)
			[varargout{1:nargout}] = min_or_max(this,'min',varargin{:});
		end
		
		%
		% this -- the curve over which to compute the min and max
		% 
		% the function is set via C.set_function(func), not passed into
		% `minandmax`
		% 
		% extra arguments are passed along to Chebfun
		function [varargout] = minandmax(this,varargin)

			[varargout{1:nargout}] = min_or_max(this,'min',varargin{:});
			vararg_min = varargout; %copy for later concatenate
			[varargout{1:nargout}] = min_or_max(this,'max',varargin{:});
			vararg_max = varargout; %copy for later concatenate
			
			% merge
			varargout = cell(1,nargout);
			for ii = 1:nargout
				varargout{ii} = [vararg_min{ii}; vararg_max{ii}];
			end
% 			positions = [m_positions;M_positions];
		end
		
		
		
		

		% computes the roots of a n-ary function over the curve.
		%
		% returns with multiplicity ~~2x what it should be 
		% (1x on the boundary, 2x not on the boundary)
		%
		% the function is set previously with `set_function`, not passed into
		% a call to `roots`.
		%
		% the output is [positions, function_values].
		% positions is a row-vector array.  that is, rows are points, and
		% columns are variables.  
		%
		% `function_values` is probably useless, but it lets you check the
		% function values at the found roots.
		function [positions, function_values] = roots(this, varargin)
			
			this.require_have_function();
			
			f = this.func;
			this.Log('running algorithm `%s` for function %s\n','roots',func2str(f));
			
			if nargout == 2
				need_function_values = true;
			else 
				need_function_values = false;
			end
			
			if (~this.have_chebfuns)
				this.ComputeChebfuns()
			else
				this.Log('reusing stored chebfuns');
			end
			
			function_values = [];
			positions = [];
			function_values_by_edge = cell(this.BRinfo.num_edges,2);
			positions_by_edge = function_values_by_edge;
			ind = 1:this.BRinfo.num_variables-1; %stupid stupid stupid -1 for homogenizing var

			for ii = 1:this.BRinfo.num_edges
				
				this.Log('getting %s for edge %i\n','roots',ii);
				
				if is_degenerate(this.BRinfo.edges(ii,:))
					this.Log('skipping degenerate edge %i\n',ii);
					continue;
				end
				
				this.SetupEdgeForTracking(ii);
				
				evalme = 'useme = f(';
				for jj = ind 
					evalme = sprintf('%s this.edge_chebfuns(%i,%i).left',evalme,ii,jj);
					if jj ~= ind(end)
						evalme = sprintf('%s, ',evalme);
					else
						evalme = sprintf('%s);',evalme);
					end
				end
				eval(evalme);
				
				pos = roots(useme,varargin{:});
			
				callme = MakeBrCfHandle(this,ii,ind,'left');
				pts_on_curve = callme(pos);
				positions_by_edge{ii,1} = pts_on_curve;
				positions = [positions;pts_on_curve]; %#ok<AGROW>
					
				
				if need_function_values
					val = useme(pos);
					function_values_by_edge{ii,1} = val;
					function_values = [function_values;val]; %#ok<AGROW>
				end
				
				
				
				
				
				% this block constructs a function handle, to be able to call
				% the desired function.  this should be abstracted somehow.
				evalme = 'useme = f(';
				ind = 1:this.BRinfo.num_variables-1;%stupid stupid stupid -1 for homogenizing var
				for jj = ind  %stupid stupid stupid -1 for homogenizing var
					evalme = sprintf('%s this.edge_chebfuns(%i,%i).right',evalme,ii,jj);
					if jj ~= ind(end)
						evalme = sprintf('%s, ',evalme);
					else
						evalme = sprintf('%s);',evalme);
					end
				end
				eval(evalme);
				
				pos = roots(useme,varargin{:});
			
				callme = MakeBrCfHandle(this,ii,ind,'right');
				pts_on_curve = callme(pos);
				positions_by_edge{ii,2} = pts_on_curve;
				positions = [positions;pts_on_curve]; %#ok<AGROW>
					
				
				if need_function_values
					val = useme(pos);
					function_values_by_edge{ii,2} = val;
					function_values = [function_values;val]; %#ok<AGROW>
				end
			end
		end % re: roots
		
		
		
		
		
		
		
		
		
		
		% overloads the plot commands for BrCfCurves
		%
		% one awesome option to this is `func` which lets you plot a given
		% function over the curve.  cool!  make sure the handle you pass is
		% n-ary, not unary
		% 
		% the not-used-here arguments are passed onto the plot call, things like
		% the LineWidth, LineStyle, etc.
		function h = plot(this, varargin)
			
			assumed_plot_args = {};
			
			

			ii=1;
			while ii <= length(varargin)
				switch varargin{ii}
					case 'color_use_indices'
						this.options.plot.color.use_indexed_vars = true;
					case 'func'
						this.set_function(varargin{ii+1});
						ii = ii+1;
					otherwise 
						assumed_plot_args{end+1} = varargin{ii}; %#ok<AGROW>
				end
				ii = ii+1;
			end
			
			have_custom_colorfun = ~isempty(this.func);
			
			if have_custom_colorfun
				titlestr = func2str(this.func);
				
				if ~this.cache.have_functions
					compute_function_per_edge(this);
				end
				colorfnhandle = this.make_color_fn_handle(this.func);
			else
				titlestr = '';
			end
			
			
						
			
			
			h = [];
			
			prev_hold = ishold();
			
			ind = this.options.plot.indices;
			num_plotvars = length(ind);
			
			switch num_plotvars
				case 2
					if have_custom_colorfun
						plotfn = @(cfs, discr, r_l) color_line(cfs{ind(1)}(discr),cfs{ind(2)}(discr),colorfnhandle(cfs,discr),assumed_plot_args{:});
						colorbar;
					else
						plotfn = @(cfs, discr, r_l) plot(cfs{ind(1)}(discr),cfs{ind(2)}(discr),assumed_plot_args{:});
					end
				case 3
					if have_custom_colorfun
						plotfn = @(cfs, discr, r_l) clinep(cfs{ind(1)}(discr),cfs{ind(2)}(discr),cfs{ind(3)}(discr),colorfnhandle(cfs, discr));
						colorbar;
					else
						plotfn = @(cfs, discr, r_l) plot3(cfs{ind(1)}(discr),cfs{ind(2)}(discr),cfs{ind(3)}(discr),assumed_plot_args{:});
					end
				otherwise
					error('plotting not enable for this many variables (%i).  please select the number of variables at construct-time by indicating the `indices` you want to plot', num_plotvars);
			end
					
			if (~this.have_chebfuns)
				this.ComputeChebfuns();
			end
			
			for ii = 1:this.BRinfo.num_edges
		
				if is_degenerate(this.BRinfo.edges(ii,:))
					this.Log('skipping degenerate edge %i\n',ii);
					continue;
				else
					this.Log('plotting edge %i (%i %i %i)\n',ii,this.BRinfo.edges(ii,1),...
																		 this.BRinfo.edges(ii,2),...
																		 this.BRinfo.edges(ii,3));
				end
				
				m = zeros(num_plotvars,2);
				for jj = 1:num_plotvars
					m(jj,1) = length(this.edge_chebfuns(ii,ind(jj)).left);
					m(jj,2) = length(this.edge_chebfuns(ii,ind(jj)).right);
				end
				
				num_left = max(m(:,1));
				num_right = max(m(:,2));
				
				left_discr = linspace(0,1,num_left);
				right_discr = linspace(0,1,num_right);
				
				h(end+1) = plotfn(sided_chebfuns(this, ii, 'left'), left_discr); %#ok<AGROW>
				hold on
				h(end+1) = plotfn(sided_chebfuns(this, ii, 'right'), right_discr); %#ok<AGROW>
				
							
				hold on
			end
			
			
			title(titlestr);
			
			if prev_hold
				hold on
			else
				hold off
			end
		end
		
		% clears the computed points, and resets to correct size
		function ResetComputedPoints(this)
			
			this.edge_computed_points = cell(this.BRinfo.num_edges,2); 
			for ii = 1:this.BRinfo.num_edges
				
				for jj = 1:2
					this.edge_computed_points{ii,jj} = containers.Map('KeyType','double','ValueType','any');
				end
				% add points for 0 and 1
				% left
				e = this.BRinfo.edges(ii,:);
				this.edge_computed_points{ii,1}(0.0) = real(this.BRinfo.vertices(e(2)).point);
				this.edge_computed_points{ii,1}(1.0) = real(this.BRinfo.vertices(e(1)).point);
				
				% right
				this.edge_computed_points{ii,2}(0.0) = real(this.BRinfo.vertices(e(2)).point);
				this.edge_computed_points{ii,2}(1.0) = real(this.BRinfo.vertices(e(3)).point);
			end
		end
		
		% clears the chebfuns, and resets to all empty
		function ResetChebfuns(this)
			
			this.edge_chebfuns = repmat(...
				repmat(...
					struct('left',chebfun(),'right',chebfun()),[this.BRinfo.num_edges 1]),...
					[1 this.BRinfo.num_variables-1]);
			this.have_chebfuns = 0;
		end
		
		% sets up the stored chebfuns
		function ComputeChebfuns(this)
			this.ResetChebfuns();
			this.ReComputeEmptyChebfuns();
		end
		
		
		% sets up the stored chebfuns
		function ReComputeEmptyChebfuns(this)
			
			SetCfArgs(this);
			
			for ii = 1:this.BRinfo.num_edges
			
				if is_degenerate(this.BRinfo.edges(ii,:))
					this.Log('skipping degenerate edge %i\n',ii);
					continue;
				else
					this.Log('computing chebfuns for edge %i (%i %i %i)\n',ii,this.BRinfo.edges(ii,1),...
																		 this.BRinfo.edges(ii,2),...
																		 this.BRinfo.edges(ii,3));
				end
				
				this.SetupEdgeForTracking(ii);
				
				ind = 1:this.BRinfo.num_variables-1; %eww, this off by one because of the homogenizing variable thing sucks.
				for jj = ind
					if isempty(this.edge_chebfuns(ii,jj).left)
						this.Log('\tcomputing chebfun for variable %i, left\n',jj);
						try
							this.edge_chebfuns(ii,jj).left = chebfun(MakeBrCfHandle(this,ii,jj,'left'),this.options.chebfunargs{:});
						catch ER
							warning(getReport(ER))
							this.Log('%s',warning(getReport(ER)));
							warning('chebfun left for edge %i, variable %i failed to recompute\n',ii,jj);
							this.Log('chebfun left for edge %i, variable %i failed to recompute\n',ii,jj);
						end
					else
						this.Log('not recomputing left');
					end
					
					if isempty(this.edge_chebfuns(ii,jj).right)
						this.Log('\tcomputing chebfun for variable %i, right\n',jj);
						try
							this.edge_chebfuns(ii,jj).right = chebfun(MakeBrCfHandle(this,ii,jj,'right'),this.options.chebfunargs{:});
						catch ER
							warning(getReport(ER))
							this.Log('%s',warning(getReport(ER)));
							warning('chebfun right for edge %i, variable %i failed to recompute\n',ii,jj);
							this.Log('chebfun right for edge %i, variable %i failed to recompute\n',ii,jj);
							pause
						end
					else
						this.Log('not computing right');
					end
				end % re: jj
				
			end
			this.have_chebfuns = true;
			this.Log('done computing chebfuns');
		end
		
		
		
		% logs to the log file, the content.  
		% vargargin is expanded into an sprintf statement.
		function this = Log(this, varargin)
			
			if length(varargin)==1
				fprintf(this.logfile.fid,'%s\n',varargin{1});
			else
				fprintf(this.logfile.fid,varargin{1},varargin{2:end});
			end
			
			if this.options.echo
				if length(varargin)==1
					fprintf('%s\n',varargin{1});
				else
					fprintf(varargin{1},varargin{2:end});
				end
			end
			
		end
		
		
		function s = saveobj(obj)
			s.br_filename = obj.br_filename;
			s.BRinfo = obj.BRinfo;
			s.b_input = copy(obj.b_input);
			s.options = obj.options;
			s.logfile = obj.logfile;
			s.edge_chebfuns = obj.edge_chebfuns;
			s.have_chebfuns = obj.have_chebfuns;
			s.edge_computed_points = obj.edge_computed_points;
			s.current_working_edge = obj.current_working_edge;
			s.cache = obj.cache;
			s.func = obj.func;
		end

		
	end % re: public member methods
	
	methods (Static)
		
		function ValidateBRinfo(BRinfo)
			if BRinfo.dimension ~= 1
				error('loaded BRinfo file does not contain a curve.');
			end

			if isempty(BRinfo.cycle_numbers)
				error('loaded curve does not have cycle numbers stored\nplease run a version of Bertini_real which computes cycle numbers.  This is available starting in 1.2.0');
			end
		end
		
		
		function newObj = loadobj(s)
			if isstruct(s)
				newObj = BrCfCurve('autoload',false,'autosetup',false); 

				newObj.br_filename = s.br_filename;
				newObj.BRinfo = s.BRinfo;
				newObj.b_input = copy(s.b_input);
				newObj.options = s.options;
				newObj.edge_chebfuns = s.edge_chebfuns;
				newObj.have_chebfuns = s.have_chebfuns;
				newObj.edge_computed_points = s.edge_computed_points;
				newObj.current_working_edge = s.current_working_edge;
				
				newObj.cache = s.cache;
				newObj.func = s.func;
			else
				newObj = s;
			end
			
			newObj.NewLogfile();
		end
	end % re: public static methods
   
	
	%%%%%%%%%%%%%%%%%%%   PRIVATE METHODS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	methods (Access = private)
		
		
		%intializes options for BrCfCurve
		function SetDefaultOptions(this)
			
			this.logfile.level = 1;
			
			
			this.options.tol = 1e-14; %todo document this setting
			this.options.cf.fixedlength = nan; % cf setting
			this.options.cf.minsamples = 9; % cf setting
			this.options.cf.maxlength = 129; % cf setting 2^16+1 is the chebfun default, which is really large
			this.options.cf.sampletest = 1; % cf setting
			this.options.cf.splitting = 'off'; % cf setting
			
			this.options.plot.indices = nan; %nan indicates want all indices
			this.options.plot.color.use_indexed_vars = false;
			
			this.options.autoload = true;
			this.options.autosetup = true;
			
			this.options.echo = false;
		end
		
		
		
		function InitialSetup(this)
			this.ResetComputedPoints();
			this.BertiniSetup();
		end
		
		
		
		%parses varargin for the constructor of BrCfCurve
		function SetOptions(this,command_line_options)
			
			if mod(length(command_line_options),2)~=0
				error('must have option-value pairs.  parity incorrect');
			end
			
			
			for ii = 1:2:length(command_line_options)-1
				switch lower(command_line_options{ii})
					case 'indices'
						this.options.plot.indices = command_line_options{ii+1};
					case 'filename'
						this.br_filename = command_line_options{ii+1};
						
					case 'fixedlength'
						this.options.cf.fixedlength = command_line_options{ii+1};
					case 'minsamples'
						this.options.cf.minsamples = command_line_options{ii+1};
					case 'maxlength'
						this.options.cf.maxlength = command_line_options{ii+1};
					case 'sampletest'
						this.options.cf.sampletest = command_line_options{ii+1};
					case 'splitting'
						this.options.cf.splitting = command_line_options{ii+1};
					case 'autoload'
						this.options.autoload = command_line_options{ii+1};
					case 'autosetup'
						this.options.autosetup = command_line_options{ii+1};
					case 'echo'
						this.options.echo = command_line_options{ii+1};
					case 'func'
						this.set_function(command_line_options{ii+1});
					otherwise
						error('bad option name for BrCfCurve ''%s''',command_line_options{ii});
				end
			end	
		end
		
		
		% translate internally stored options for chebfun into a cell array, to pass
		% on to chebfun when it is called.
		function SetCfArgs(this)
			this.options.chebfunargs = {...
				'domain',[0 1],...
				'splitting',this.options.cf.splitting,...
				'sampleTest',this.options.cf.sampletest,...
				'minsamples',this.options.cf.minsamples,'maxLength',this.options.cf.maxlength,...
				'fixedLength',this.options.cf.fixedlength};
		end
				
		% creates a function handle to a given function.  this is necessary to
		% deal with n-ary functions, since they are required to work with this
		% class type.
		function handle = SubsFuncOnSide(this, edge_ind, func, side) %#ok<INUSL,STOUT>
			evalme = 'handle = func(';
			for jj = 1:this.BRinfo.num_variables-1  %stupid stupid stupid -1 for homogenizing var
				evalme = sprintf('%s this.edge_chebfuns(%i,%i).%s',evalme,edge_ind,jj,side);
				if jj ~= this.BRinfo.num_variables-1
					evalme = sprintf('%s, ',evalme);
				else
					evalme = sprintf('%s);',evalme);
				end
			end
			eval(evalme);
		end
		
		
		
		% generates errors if the handle is not compatible.  the function must
		% be a function_handle, and accept n arguments -- that is, f must be
		% n-ary, where n is the number of variables
		function ValidateHandle(this, func)
			
			if ~or( ishandle(func), isa(func,'function_handle')) %really matlab, anonymous functions aren't handles, but they are function_handles?  really?
				error('argument must be a function handle');
			end
			
			if (nargin(func) ~= this.BRinfo.num_variables-1)
				error('function handle does not accept vector arguments of length %i',this.BRinfo.num_variables-1);
			end
			
		end
		
		
		
		% sets up a function handle which Chebfun will use to compute
		% points on the curve.
		%
		% also can be used to compute points on the curve without chebfun
		% driving.
		function h = MakeBrCfHandle(this,edge_index, coord_indices, side)
	
			b = this.BRinfo;
			ed = b.edges(edge_index,:);
			if strcmp(side,'left')
				pi_out = b.vertices(ed(1)).projection_value;
				c = b.cycle_numbers(edge_index,1);
				s = 1;
			elseif strcmp(side,'right')
				pi_out = b.vertices(ed(3)).projection_value;
				c = b.cycle_numbers(edge_index,2);
				s = 2;
			else
				error('bad value %s for ''side''',side);
			end
			pi_mid = b.vertices(ed(2)).projection_value;
			
			if this.logfile.level >= 1
				this.Log('\t\tmaking ChebFun handle for edge %i with cycle number %i\n',edge_index,c);
			end
			
			h = @(x) call_bertini_chebfun(x,...
					pi_out,pi_mid,...
					c,...
					coord_indices,...
					this.BRinfo.num_variables-1,...
					this.edge_computed_points{edge_index,s},...
					this.logfile.level, this.logfile.fid);
		end
		
		
		
		
		%writes the bertini input files and parameter files for doing
		%tracking on an edge of the curve.
		% also stores the edge index into the curve's state.
		function SetupEdgeForTracking(this, edge_index)
			b = this.b_input;
			
			midpoint_index = this.BRinfo.edges(edge_index,2);

			p_start = this.BRinfo.vertices(midpoint_index).projection_value;
			
			write_generic_solns(p_start,'start_parameters');
			write_generic_solns(p_start,'final_parameters');
			write_generic_solns(real(this.BRinfo.vertices(midpoint_index).point),'start');
			
			b.config.deletetempfiles = 0;
			
			b.config.needtodiff = 1;
			write_bertini_input_file(b.variable_group, b.functions, 'filename', 'input_BrCf_slice','options',b.config,'constants',b.constant,'subfunctions',b.subfunction,'parameters',b.parameter);
			bertini('filename','input_BrCf_slice','stifle');
			
			b.config.needtodiff = 2;
			write_bertini_input_file(b.variable_group, b.functions, 'filename', 'input_BrCf_slice','options',b.config,'constants',b.constant,'subfunctions',b.subfunction,'parameters',b.parameter);
			
			this.current_working_edge = edge_index;
		end
		
		
				% sets up an input file, for bertini to use for tracking to the
		% points that chebfun tells it to.
		function [this] = SetupInput(this)
			this.b_input.parse_from_string(this.BRinfo.input);
			
			this.b_input.config.parameterhomotopy = 2;
			this.b_input.config.tracktype = 0;
			if isfield(this.b_input.config,'userhomotopy')
				this.b_input.config = rmfield(this.b_input.config,'userhomotopy');
			end
			
			this.b_input.config.randomseed = 1;
			this.b_input.config.mptype = 2;
			this.b_input.config.tracktolbeforeeg = 1e-7;
			this.b_input.config.tracktolduringeg = 1e-7;
			this.b_input.config.finaltol = 1e-11;
			this.b_input.config.odepredictor = 8;
			this.b_input.config.endgamenum = 1;
			this.b_input.config.endgamebdry = 0.0001;
			this.b_input.config.numsamplepoints = 5;
			this.b_input.config.maxstepsbeforenewton = 0;
			this.b_input.config.maxnewtonits = 1;
			this.b_input.config.sharpendigits = 20;
			this.b_input.config.condnumthreshold = 1e30;
		end
		
		
		% sets up the bertini_input object in this curve to be ready to
		% slice at a given projection value.
		function [] = SetupSliceSystem(this)
			vars = this.b_input.variable_group; % this is a cell array
			p = this.BRinfo.pi;
			b = this.b_input;
			
			for ii = 1:size(vars)
				b.declare_and_define(sprintf('proj_val%i',ii),p(ii),'constant');
			end
			
			b.declare_symbols({'target_proj_val'},'parameter');
			
			
			s = sprintf('-target_proj_val');
			for ii = 1:size(vars)
				s = sprintf('%s + %s*proj_val%i',s,vars{ii},ii');
			end
			
			b.declare_and_define('br_cf_curve_slice',s,'function');
			write_bertini_input_file(b.variable_group, b.functions, 'filename', 'canonical_input_BrCf_slice','options',b.config,'constants',b.constant,'subfunctions',b.subfunction,'parameters',b.parameter);
		end
		
		function [] = compute_function_per_edge(this)
			
			this.require_have_function();
			
			this.Log('computing function chebfuns for each half-edge');
			this.cache.have_functions = false;
			this.cache.functions = cell(this.BRinfo.num_edges,2);
			
			if (~this.have_chebfuns)
				this.ComputeChebfuns();
			else
				this.Log('reusing stored chebfuns');
			end
			
			edge_range = 1:this.BRinfo.num_edges; 
			for ii = edge_range
				
				if is_degenerate(this.BRinfo.edges(ii,:))
					this.Log('skipping degenerate edge %i\n',ii);
					continue;
				else
					this.Log('computing functions for edge %i\n',ii);
				end
				
				% run the min/max algorithm
				this.Log('left half of edge');
				this.cache.functions{ii,1} = this.SubsFuncOnSide(ii,this.func,'left');
				
				this.Log('right half of edge');
				this.cache.functions{ii,2} = this.SubsFuncOnSide(ii,this.func,'right');

			end %re: edge loop
			this.cache.have_functions = true;
			
		end
		
		
		
		% a generic function, intended to be called by min, max, or
		% minandmax.  wraps around chebfun stuffs, and returns the function
		% values, the natural coordinates on the curve of the max and min.
		% also produces edgewise data if requested.
		function [varargout] = min_or_max(this, mode, varargin)
			
			switch mode
				case 'min'
					m_or_m = @min;
				case 'max'
					m_or_m = @max;
				otherwise
					error('bad option %s for mode', mode);
			end
			
			this.require_have_function();
			
			if ~this.cache.have_functions
				compute_function_per_edge(this);
			end

			f = this.func;
			
			
			% initialize
			
			this.Log('running algorithm %s for function %s\n',mode,func2str(f));
			
			
			need_coords = (nargout > 1);
			
			values = [];
			value_edge_index = [];
			coordinates = [];
			values_by_edge = cell(this.BRinfo.num_edges,2);
			coordinates_by_edge = values_by_edge;
			positions_by_edge = values_by_edge; % used to hold the chebfun values of the min/max
			
			% loop over each edge.  
			% skip degenerate ones, 
			% but this is wrong because there may be isolated 
			% singular points with degenerate edges
			
			edge_range = 1:this.BRinfo.num_edges; 
			for ii = edge_range
				
				if is_degenerate(this.BRinfo.edges(ii,:))
% 					this.Log('skipping degenerate edge %i\n',ii);
					continue;
				else
					this.Log('getting %s for edge %i\n',mode,ii);
				end
				
				% run the min/max algorithm
				this.Log('left half of edge');
				[val_left,pos_left]   = m_or_m(this.cache.functions{ii,1},varargin{:});
				this.Log('right half of edge');
				[val_right,pos_right] = m_or_m(this.cache.functions{ii,2},varargin{:});
				
				
				% now check for mid-edge nonsense.
				tol = this.options.tol;
				if any(abs(pos_left)<tol) % then the midpt is labeled as min or max
					this.Log('midpoint of edge is a min/max, from the left');
					this.Log('removing left val %f, edge %i\n',val_left(abs(pos_left)<tol),ii);
					ind = abs(pos_left)<tol;
					val_left(ind) = [];
					pos_left(ind) = [];
				end
				
				if any(abs(pos_right)<tol) % then the midpt is labeled as min or max
					this.Log('midpoint of edge is a min/max, from the right');
					this.Log('removing right val %f, edge %i\n',val_right(abs(pos_right)<tol),ii);
					ind = abs(pos_right)<tol;
					val_right(ind) = [];
					pos_right(ind) = [];
				end
				
				values_by_edge{ii,1} = val_left;
				positions_by_edge{ii,1} = pos_left;
				values_by_edge{ii,2} = val_right;
				positions_by_edge{ii,2} = pos_right;
				
			end %re: computing m/m edge loop
			
			
			
			this.Log('\n\tentering post processing\n\n',[]);
			% some post-processing
			% we must remove the edge boundary false min/max points
			% we know they are false if not all edges with this point as a
			% boundary point have it as a min/max.
			% we have to check ALL (nondegenerate) edges.
			for ii = edge_range
				if is_degenerate(this.BRinfo.edges(ii,:))
					continue;
				end
				this.Log('\nedge %i\n',ii);
				pos_left = positions_by_edge{ii,1};
				val_left = values_by_edge{ii,1};
				
				pos_right = positions_by_edge{ii,2};
				val_right = values_by_edge{ii,2};
				
				
				tol = this.options.tol;
				
				
				%%% remove spurious m/m on the left side
				
				%%%%%%%%%
				if any(abs(pos_left-1)<tol) % then the endpt of left edge is labeled as min or max
					other_edge_has_this_mm_at_bdry = zeros(this.BRinfo.num_edges,1); %preallocate
					
					non_degenerate = ~is_degenerate(this.BRinfo.edges); %get a boolean vector indicating the nondegenerate edges
					
					has_as_bound = this.BRinfo.edges(:,1) == this.BRinfo.edges(ii,1); % boolean vector indicating those edges with the left bdry pt on its left bdry
					consider_me_l = all([has_as_bound non_degenerate],2);
					for jj = 1:this.BRinfo.num_edges
						if or(jj==ii, ~consider_me_l(jj))
							continue;
						end
						
						if any(abs(positions_by_edge{jj,1}-1)<tol) %check the left
							other_edge_has_this_mm_at_bdry(jj) = true;
						end
					end

					has_as_bound = this.BRinfo.edges(:,3) == this.BRinfo.edges(ii,1); % boolean vector indicating those edges with the left bdry pt on its left bdry
					consider_me_r = all([has_as_bound non_degenerate],2);
					for jj = 1:this.BRinfo.num_edges
						if or(jj==ii, ~consider_me_r(jj))
							continue;
						end
						
						if any(abs(positions_by_edge{jj,2}-1)<tol) %check the right half-edge
							other_edge_has_this_mm_at_bdry(jj) = true;
						end
					end
					
					if and(~any(other_edge_has_this_mm_at_bdry),sum(or(consider_me_l,consider_me_r))>1) %if the left edge doesnt appear as a max
							% then we can eliminate it as a min/max.  

							ind = abs(pos_left-1)<tol;

							this.Log('removing m/m at left boundary with val %f, edge %i\n',...
								val_left(ind),ii);

							pos_left(ind) = [];
							val_left(ind) = [];
					end
				end % if any(pos_left-1)...
				positions_by_edge{ii,1} = pos_left;
				values_by_edge{ii,1} = val_left;
				%%%%%%%%%%%%%%%%
				
				
				%%% now to remove spurious m/m on the right side
				
				if any(abs(pos_right-1)<tol) % then the endpt of left edge is labeled as min or max
					other_edge_has_this_mm_at_bdry = zeros(this.BRinfo.num_edges,1); %preallocate
					
					non_degenerate = ~is_degenerate(this.BRinfo.edges); %get a boolean vector indicating the nondegenerate edges
					
					has_as_bound = this.BRinfo.edges(:,1) == this.BRinfo.edges(ii,3); % boolean vector indicating those edges with the left bdry pt on its left bdry
					consider_me_l = all([has_as_bound non_degenerate],2);
					for jj = 1:this.BRinfo.num_edges
						if or(jj==ii, ~consider_me_l(jj))
							continue;
						end
						
						if any(abs(positions_by_edge{jj,1}-1)<tol) %check the left
							other_edge_has_this_mm_at_bdry(jj) = true;
						end
					end

					has_as_bound = this.BRinfo.edges(:,3) == this.BRinfo.edges(ii,3); % boolean vector indicating those edges with the left bdry pt on its left bdry
					consider_me_r = all([has_as_bound non_degenerate],2);
					for jj = 1:this.BRinfo.num_edges
						if or(jj==ii, ~consider_me_r(jj))
							continue;
						end
						
						if any(abs(positions_by_edge{jj,2}-1)<tol) %check the right half-edge
							other_edge_has_this_mm_at_bdry(jj) = true;
						end
					end
					
					
					if and(~any(other_edge_has_this_mm_at_bdry),sum(or(consider_me_l,consider_me_r))>1) %if the left edge doesnt appear as a max
							% then we can eliminate it as a min/max.  

							ind = abs(pos_right-1)<tol;

							this.Log('removing m/m at right boundary with val %f, edge %i\n',...
								val_right(ind),ii);

							pos_right(ind) = [];
							val_right(ind) = [];
					end
				end % if any(pos_left-1)...
				positions_by_edge{ii,2} = pos_right;
				values_by_edge{ii,2} = val_right;
			end % postprocessing for edges

			
			
			
			
			% final assembly of the data, and computation of the points on 
			% the curve if needed
			
			% flattens the data into a rectangular array, so we can return it
			if need_coords
				this.Log('computing actual coordinates of computed optima');
			end
			
			for ii = edge_range
				if is_degenerate(this.BRinfo.edges(ii,:))
					continue;
				end
				
				value_edge_index = [value_edge_index;ii*ones(length(values_by_edge{ii,1}),1)]; %#ok<AGROW>
				value_edge_index = [value_edge_index;ii*ones(length(values_by_edge{ii,2}),1)]; %#ok<AGROW>
				values = [values;values_by_edge{ii,1}]; %#ok<AGROW>
				values = [values;values_by_edge{ii,2}]; %#ok<AGROW>
				
				if need_coords
					
					num_left_optima = length(positions_by_edge{ii,1});
					num_right_optima = length(positions_by_edge{ii,2});
					
					if or(num_left_optima>0, num_right_optima>0)
						this.SetupEdgeForTracking(ii);
					end
					
					if num_left_optima>0
						callme = MakeBrCfHandle(this,ii,1:this.BRinfo.num_variables-1,'left');
						this.Log('\tedge %i, left side, with %i optima on it\n',ii,num_left_optima);
						pts_on_curve = callme(positions_by_edge{ii,1});
						coordinates_by_edge{ii,1} = pts_on_curve; % 1 is left
					else
						coordinates_by_edge{ii,1} = [];
					end
					
					if num_right_optima>0
						callme = MakeBrCfHandle(this,ii,1:this.BRinfo.num_variables-1,'right');
						this.Log('\tedge %i, right side, with %i optima on it\n',ii,num_right_optima);
						pts_on_curve = callme(positions_by_edge{ii,2});
						coordinates_by_edge{ii,2} = pts_on_curve; % 2 is right
					else
						coordinates_by_edge{ii,2} = [];
					end
					
					coordinates = [coordinates;coordinates_by_edge{ii,1}]; %#ok<AGROW>
					coordinates = [coordinates;coordinates_by_edge{ii,2}]; %#ok<AGROW>
				end
			end % for ii=edge_range, final assembly for edges
			
			
			% if things aren't local...
			% then we just grab the min or max value and its location
			if ~any(ismember(varargin,'local'))
				[m,ind] = m_or_m(values);
				values = m;
				value_edge_index(ind)
				if need_coords
					coordinates = coordinates(ind,:);
				end
				
			end % if ~any
			
			switch nargout
				case 1
					varargout = {values};
				case 2
					varargout = {values, coordinates};
				case 3
					varargout = {values, coordinates, values_by_edge};
				case 4
					varargout = {values, coordinates, values_by_edge, coordinates_by_edge};
			end
		end % re: function min_or_max

	
		% takes a function handle and makes it callable from a 1xn cell array of
		% chebfuns
		function h = make_color_fn_handle(this, colorfn) %#ok<INUSD>  
			h = ''; % to shut up the linter
			evalme = 'h = @(cfs, discr) colorfn('; %cfs is for chebfuns
			
			if this.options.plot.color.use_indexed_vars
				ind = this.options.plot.indices;
			else
				ind = 1:this.BRinfo.num_variables-1; %-1 for homogenizing coordinate.  yep.  it's annoying.  regrets are strong.
			end
			

			for ii = 1:length(ind)
				marker = ',';
				if ii == length(ind)
					marker = ')';
				end  %cfs is for chebfuns.  discr is for discretization
				evalme = sprintf('%scfs{%i}(discr)%s',evalme,ind(ii),marker);
			end
			eval(evalme);
		end

		% extracts into a 1xn cell array the chebfuns for half of an edge.
		%
		% `side` should be 'right' or 'left'
		function cfs = sided_chebfuns(this, edge, side)
			cfs = cell(1,this.BRinfo.num_variables-1);
			for ii = 1:this.BRinfo.num_variables-1
				cfs{ii} = this.edge_chebfuns(edge,ii).(side)(:);
			end
		end
	end % re: private methods






end
