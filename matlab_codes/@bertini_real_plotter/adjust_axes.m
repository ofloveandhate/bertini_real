function br_plotter = adjust_axes(br_plotter)

curr_axes = br_plotter.axes.main;


lower = min(br_plotter.fv.vertices,[],1);
upper = max(br_plotter.fv.vertices,[],1);


ten_percent = 0.10 *(upper-lower);

ten_percent(ten_percent<0.01) = 0.01;


lower = lower - ten_percent;
upper = upper + ten_percent;

set(curr_axes,'XLim',[lower(1) upper(1)]);
set(curr_axes,'YLim',[lower(2) upper(2)]);
br_plotter.cam_pos = [upper(1) lower(2)];%(upper(2)+lower(2))/2];
% BRinfo.center(plot_params.ind(2))-1.1*BRinfo.radius BRinfo.center(plot_params.ind(2))+1.1*BRinfo.radius
if size(br_plotter.fv.vertices,2)==3
    set(curr_axes,'ZLim',[lower(3) upper(3)]);
	br_plotter.cam_pos = [br_plotter.cam_pos upper(3)];
end



end
