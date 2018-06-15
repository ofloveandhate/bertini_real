function br_plotter = save_routine(br_plotter,varargin) % additional arguments go afer the first two, which must be put here ,even if they are not used.
%perhaps get info from other calls here?

hide_panels(br_plotter)

try
	if and(br_plotter.options.render_faces,or(br_plotter.switches.display_faces == 1,br_plotter.switches.display_face_samples == 1))
		br_plotter.options.format = 'png';
		br_plotter.options.format_flag = 'png';
		br_plotter.options.resolution = 600;
		render_into_file(br_plotter.options);
	else
		br_plotter.options.format = 'png';
		br_plotter.options.format_flag = 'png';
		br_plotter.options.resolution = 600;
		render_into_file(br_plotter.options);
	end
catch e
	msg = getReport(e);
	warning('unable to complete saving for this reason:\n\n%s',msg);
end

show_panels(br_plotter)

end
