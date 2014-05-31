function sphere_plot(BRinfo)
global plot_params


center = real(BRinfo.center(plot_params.ind));
radius = real(BRinfo.radius);

switch length(plot_params.ind)

	case 2
		num_samp = 100;
		to_2_pi = linspace(0,2*pi,num_samp);
		x = radius*cos(to_2_pi) + center(1);
		y = radius*sin(to_2_pi) + center(2);
		plot_params.handles.sphere = plot(x,y,':k','LineWidth',2);
		
	case 3
		num_samp = 20;
		[x,y,z] = sphere(num_samp);
		h = patch(surf2patch(radius*x+center(1),radius*y+center(2),radius*z+center(3)));
		set(h,'FaceAlpha',0.03,'EdgeColor','none','FaceColor','k');
		plot_params.handles.sphere = h;
	otherwise
		display('bad length(ind) in sphere_plot.')
		return;
end

switch BRinfo.dimension
	case 1
		% i have no idea right now how to tell from the data i already
		% collect whether this is bounded or not.
			plot_params.is_bounded = 0;
		
	case 2
		if all( and(BRinfo.sphere_curve.edges(:,1)==BRinfo.sphere_curve.edges(:,2),BRinfo.sphere_curve.edges(:,3)==BRinfo.sphere_curve.edges(:,2)))
			plot_params.is_bounded = 1;
		else
			plot_params.is_bounded = 0;
		end
	otherwise
		
end



end