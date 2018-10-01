function br_plotter = save_routine(br_plotter,varargin) % additional arguments go afer the first two, which must be put here ,even if they are not used.
%perhaps get info from other calls here?

hide_panels(br_plotter)

try
	if and(br_plotter.options.render_faces,or(br_plotter.switches.display_faces == 1,br_plotter.switches.display_face_samples == 1))
		br_plotter.options.figure_saving.format = 'png';
		br_plotter.options.figure_saving.format_flag = 'png';
		br_plotter.options.figure_saving.resolution = 600;
		render_into_file(br_plotter.options.figure_saving);
	else
		br_plotter.options.figure_saving.format = 'png';
		br_plotter.options.figure_saving.format_flag = 'png';
		br_plotter.options.figure_saving.resolution = 600;
		render_into_file(br_plotter.options.figure_saving);
	end
catch e
	msg = getReport(e);
	warning('unable to complete saving for this reason:\n\n%s',msg);
	
	fprintf('you probably just need to clone and add to the path https://github.com/ofloveandhate/brakelab');
end

show_panels(br_plotter)

end
