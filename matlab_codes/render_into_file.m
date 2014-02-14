% render_into_file(basename, plot_params)
%
% a utility for saving plots to a file from matlab.
% intended for use with the paramotopy suite of functions.
%
% input:
%  plotmode -- a string which forms the basis of a filename

%  plot_params -- a structure with data members for setting up the save.
%      members:  fontsize, an integer
%                window, a 4-element vector in standard matlab format for
%                     specifying the window into which to plot.
%                format, a string for the format to save.  see the print
%                         help.
%                format_flag, again a string for the save driver to use.
%
%  outputs a file of the appropriate format to the current working dir.
%
%
%

% daniel brake
% colorado state university
% mathematics
% 2013
% danielthebrake@gmail.com


function plot_params = render_into_file(varargin)
%


if isempty(varargin)
	
	plot_params.fontsize = 16;
	plot_params.window = get(gcf,'Position');
	plot_params.format = 'eps';
	plot_params.format_flag = 'epsc2';
	plot_params.basename = 'default_filename';
elseif and(ischar(varargin{1} ), length(varargin)==1)
	
	plot_params.fontsize = 16;
	plot_params.window = get(gcf,'Position');
	plot_params.format = 'eps';
	plot_params.format_flag = 'epsc2';
	
	plot_params.basename = varargin{1};
	
elseif isstruct(varargin{1})
	plot_params = varargin{1};
else
	plot_params.fontsize = 16;
	plot_params.window = get(gcf,'Position');
	plot_params.format = 'eps';
	plot_params.format_flag = 'epsc2';
	plot_params.basename = 'default_filename';
end



fig1 = gcf;

set(fig1,'PaperPositionMode','auto');
set(gca,'FontSize',plot_params.fontsize-2);

% labelText('FontName','Times New Roman');

% basename = sprintf('%s_%s',plotmode,dataname(1:end-4));

currname = increment_name(plot_params.basename);
nameforfile = sprintf('%s.%s',currname,plot_params.format);

evalme = sprintf('print(fig1,''%s'',''-d%s''',nameforfile,plot_params.format_flag);
for ii = 2:length(varargin)
    evalme = sprintf('%s,''%s''',evalme,varargin{ii});
end
evalme = sprintf('%s);',evalme);
evalme
eval(evalme);
end

