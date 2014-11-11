function init_cam_pos = adjust_axes_br(vertices,curr_axes)

lower = min(vertices);
upper = max(vertices);


ten_percent = 0.10 *(upper-lower);

lower = lower - ten_percent;
upper = upper + ten_percent;

set(curr_axes,'XLim',[lower(1) upper(1)]);
set(curr_axes,'YLim',[lower(2) upper(2)]);
init_cam_pos = [upper(1) lower(2)];%(upper(2)+lower(2))/2];
% BRinfo.center(plot_params.ind(2))-1.1*BRinfo.radius BRinfo.center(plot_params.ind(2))+1.1*BRinfo.radius
if size(vertices,2)==3
    set(curr_axes,'ZLim',[lower(3) upper(3)]);
	init_cam_pos = [init_cam_pos upper(3)];
end



end