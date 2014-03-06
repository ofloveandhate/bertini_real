
function k = plot_projection(BRinfo,ind)
%
global plot_params

curr_axes = plot_params.axes.main;

pi = BRinfo.pi;

labeltext = {};
for ii = 1:BRinfo.dimension
	switch length(ind)
		case {2}
			k = plot([0 pi(ind(1),ii)],[0 pi(ind(2),ii)],'--k','Parent',curr_axes);
			
		case {3}
			k = plot3([0 pi(ind(1),ii)],[0 pi(ind(2),ii)],[0 pi(ind(3),ii)],'--k','Parent',curr_axes);
	end
	plot_params.handles.projection(ii) = k;
	labeltext{ii} = ['\pi_' num2str(ii-1)];
end

switch length(ind)
	case {2}
		texthandles = text(pi(ind(1),:),pi(ind(2),:),labeltext,'HorizontalAlignment','right','FontSize',plot_params.fontsize-2,...
			'Parent',curr_axes,'Color','r');
	case {3}
		texthandles = text(pi(ind(1),:),pi(ind(2),:),pi(ind(3),:),labeltext,'HorizontalAlignment','right','FontSize',plot_params.fontsize-2,...
			'Parent',curr_axes,'Color','r');
end

plot_params.handles.projection = [plot_params.handles.projection'; texthandles];
end



