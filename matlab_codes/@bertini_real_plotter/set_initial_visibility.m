

%labels and set visibility

function br_plotter = set_initial_visibility(br_plotter)

br_plotter.switches.legend = 0;

br_plotter.switches.label_faces = 0;
br_plotter.switches.label_spherecurve = 0;
br_plotter.switches.label_critcurve = 0;
br_plotter.switches.label_critedges = 0;
br_plotter.switches.label_midedges = 0;
br_plotter.switches.label_singular = 0;
br_plotter.switches.show_edges = 0;
br_plotter.switches.show_curve_samples = 0;

br_plotter.switches.display_vertices = 0;
br_plotter.switches.label_vertices = 0;

if br_plotter.options.render_vertices
	f = br_plotter.legend.vertices.types;

	for ii = 1:length(f)
		if strcmp('CRITICAL',f{ii})
			br_plotter.switches.vertex_set.(f{ii}) = 1;
		else
			br_plotter.switches.vertex_set.(f{ii}) = 0;
		end
	end
end

switch br_plotter.dimension
	case 2
	%         plot_params.switches.display_faces = 1;
	%         plot_params.switches.display_face_samples = 0;
		if isempty(br_plotter.handles.surface_samples)
			br_plotter.switches.display_faces = 1;
			br_plotter.switches.display_face_samples = 0; 
			
			br_plotter.switches.curve_refinements.main = 0;
			br_plotter.switches.curve_refinements.show_midslices = 0;
			br_plotter.switches.curve_refinements.show_critslices = 0;
			br_plotter.switches.curve_refinements.show_spherecurve = 0;
			br_plotter.switches.curve_refinements.show_critcurve = 0;
			br_plotter.switches.curve_refinements.show_singularcurve = 0;


		else
			br_plotter.switches.display_faces = 0;
			br_plotter.switches.display_face_samples = 1;
			
			br_plotter.switches.curve_refinements.main = 1;
			br_plotter.switches.curve_refinements.show_midslices = 0;
			br_plotter.switches.curve_refinements.show_critslices = 0;
			br_plotter.switches.curve_refinements.show_spherecurve = 1;
			br_plotter.switches.curve_refinements.show_critcurve = 1;
			br_plotter.switches.curve_refinements.show_singularcurve = 1;
			
			
		end
		
		br_plotter.switches.raw_curves_main = 0;
		br_plotter.switches.show_critcurve = 0;
		br_plotter.switches.show_spherecurve = 0;
		br_plotter.switches.show_critslices = 0;
		br_plotter.switches.show_midslices = 0;
		br_plotter.switches.show_singular = 1;
	
	case 1

		if isempty(br_plotter.handles.sample_edges)
			br_plotter.switches.show_edges = 1;
			br_plotter.switches.show_curve_samples = 0; 
		else
			br_plotter.switches.show_edges = 0;
			br_plotter.switches.show_curve_samples = 1;
		end

		br_plotter.switches.display_faces = 0;
		br_plotter.switches.display_face_samples = 0;

	otherwise
		
end


common_init(br_plotter);



end



function common_init(br_plotter)

br_plotter.switches.display_projection = 0;

br_plotter.switches.main_axes = 1;

br_plotter.switches.show_sphere = 0;


end









