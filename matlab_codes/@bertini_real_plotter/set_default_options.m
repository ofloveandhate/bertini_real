function set_default_options(br_plotter)
	br_plotter.options.use_custom_projection = false;
	br_plotter.options.markersize = 10;
	br_plotter.options.sample_alpha = 0.5;
	br_plotter.options.face_alpha = 0.5;
	br_plotter.options.raw_triangulation_alpha = 0.4;
	br_plotter.options.sample_triangulation_alpha = 0.4;
	br_plotter.options.fontsizes.legend = 12;
	br_plotter.options.fontsizes.labels = 16;
	br_plotter.options.fontsizes.axis = 20;
	br_plotter.options.linewidth = 2;
	br_plotter.options.edge_width = 1;
	br_plotter.options.autosave = false;

	br_plotter.options.labels = false;
	br_plotter.options.monocolor = false;

	br_plotter.options.render_vertices = false;
	br_plotter.options.render_curves = true;
	br_plotter.options.render_faces = true;
	br_plotter.options.render_samples = true;
	br_plotter.options.which_faces = [];
	br_plotter.options.user_indices = [];
	
	br_plotter.options.touching_edges_only = true;


	br_plotter.options.use_colorfn = false;
	br_plotter.options.colorfn_uses_raw = false;

	br_plotter.options.num_colors = 256;
	if isempty(which('parula'))
		br_plotter.options.colormap = @jet;
	else
		br_plotter.options.colormap = @parula;
	end
	
	br_plotter.options.figure_saving = render_into_file('generate');
end