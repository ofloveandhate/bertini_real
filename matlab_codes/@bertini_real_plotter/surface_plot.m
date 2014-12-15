%  surface specific code



function br_plotter = surface_plot(br_plotter)


intialize_handles(br_plotter)

create_axes(br_plotter);

label_axes(br_plotter);

sphere_plot(br_plotter);

plot_vertices(br_plotter);

adjust_axes(br_plotter);


if br_plotter.options.render_curves
	plot_surface_edges(br_plotter);
end

if br_plotter.options.render_faces
	plot_faces(br_plotter);
	plot_surface_samples(br_plotter);
end

plot_projection(br_plotter);




sync_axes(br_plotter);


end







function intialize_handles(br_plotter)

br_plotter.handles.faces = [];
br_plotter.handles.face_labels = [];

br_plotter.handles.surface_samples = [];


br_plotter.handles.critcurve_labels = [];
br_plotter.handles.spherecurve_labels = [];
br_plotter.handles.refinements.critcurve = [];
br_plotter.handles.critcurve = [];


br_plotter.handles.crittext = [];
br_plotter.handles.critslices = [];
br_plotter.handles.refinements.critslice = [];

br_plotter.handles.refinements.spherecurve = [];
br_plotter.handles.spherecurve = [];


br_plotter.handles.midtext = [];
br_plotter.handles.midslices = [];
br_plotter.handles.refinements.midslice = [];




br_plotter.handles.singtext = [];
br_plotter.handles.singular_curves = [];
br_plotter.handles.refinements.singularcurve = [];


end










