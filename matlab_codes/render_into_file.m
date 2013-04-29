% render_into_file(basename, plot_params)
%
% a utility for saving plots to a file from matlab.
% intended for use with the paramotopy suite of functions.
%
% input:
%  plotmode -- a string which forms the basis of a filename
%  dataname -- a string ending with .mat which further specifies the
%                filename
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


function plot_params = render_into_file(basename, plot_params)
%

if nargin<2
	fontsize = 16;
	window = get(gcf,'Position');
	format = 'eps';
	format_flag = 'epsc2';
else
	fontsize  = plot_params.fontsize;
	window   = plot_params.window;
	format   = plot_params.format;
	format_flag   = plot_params.format_flag;
end


if nargin < 1
	basename = 'default_filename';
end

fig1 = gcf;

set(fig1,'Position',window,'PaperPositionMode','auto');
set(gca,'FontSize',fontsize-2);

% basename = sprintf('%s_%s',plotmode,dataname(1:end-4));

currname = increment_name(basename);
nameforfile = sprintf('%s.%s',currname,format);
print(fig1,nameforfile,sprintf('-d%s',format_flag));
% 
% if input('save another view of this plot? \n if so, set the window, and press 1 then enter.');
% 	plot_params.fontsize = fontsize;
% 	plot_params.window = window;
% 	plot_params.format = format;
% 	plot_params.format_flag = format_flag;
% 	
% 	render_into_file(basename,plot_params);
% 	
% end
plot_params.fontsize = fontsize;
plot_params.window = window;
plot_params.format = format;
plot_params.format_flag = format_flag;

end