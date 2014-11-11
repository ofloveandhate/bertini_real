% render the current scene to a file, from a bertini_real_plotter object.
%
% daniel brake
% north carolina state university
% mathematics
% 2014
% danielthebrake@gmail.com


function render_into_file(br_plotter)



fig1 = br_plotter.figures.main;



currname = increment_name(br_plotter.options.basename);
nameforfile = sprintf('%s.%s',currname,br_plotter.format);

try
	print(fig1,nameforfile,['-d' br_plotter.format_flag]);
	display(sprintf('saved scene with filename: %s',nameforfile));
catch
	display('rendering the picture failed for some reason');
end


end

