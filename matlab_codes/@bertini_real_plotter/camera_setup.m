function br_plotter = camera_setup(br_plotter)

cameratoolbar

curr_axes = br_plotter.axes.main;


set(curr_axes,'CameraTargetMode','manual');
br_plotter.scene.target = real(br_plotter.BRinfo.center(br_plotter.indices));
set(curr_axes,'CameraTarget',br_plotter.scene.target);


% init_campos = real(BRinfo.center(plot_params.ind));
% init_campos = init_campos + 9*BRinfo.radius;

br_plotter.scene.campos = br_plotter.cam_pos;

set(curr_axes,'CameraPosition',br_plotter.scene.campos);
set(curr_axes,'CameraViewAngle',70);

rotate3d off

end