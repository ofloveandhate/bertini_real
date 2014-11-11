
function create_axes_br()
global plot_params


curr_axes = gca;

set(curr_axes,'visible','on');
set(curr_axes,'FontSize',plot_params.fontsize);

hold(curr_axes,'on');



set(curr_axes,'dataaspectratio',[1 1 1]);

set(curr_axes,'Tag','bertini_real');



set(curr_axes,'Position', [0.1 0.1 0.8 0.8]);
plot_params.axes.main = curr_axes;



end
