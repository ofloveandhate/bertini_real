
function k = plot_projection(BRinfo,ind)
	%
	global plot_params

	curr_axis = plot_params.axes.projection;

	pi = BRinfo.pi;

    labeltext = {};
	for ii = 1:BRinfo.dimension
        switch length(ind)
            case {2}
                k = plot([0 pi(ind(1),ii)],[0 pi(ind(2),ii)],'--k','Parent',curr_axis);

            case {3}
                k = plot3([0 pi(ind(1),ii)],[0 pi(ind(2),ii)],[0 pi(ind(3),ii)],'--k','Parent',curr_axis);
        end
        labeltext{ii} = ['\pi_' num2str(ii-1)];
	end
    
    switch BRinfo.dimension
        case {2}
            plot_params.handles.projection_text = text(pi(ind(1),:),pi(ind(2),:),labeltext,'HorizontalAlignment','right','FontSize',plot_params.fontsize-2,...
			'Parent',plot_params.axes.projection,'Color','r');
        case {3}
            plot_params.handles.projection_text = text(pi(ind(1),:),pi(ind(2),:),pi(ind(3),:),labeltext,'HorizontalAlignment','right','FontSize',plot_params.fontsize-2,...
			'Parent',plot_params.axes.projection,'Color','r');
    end
    
end



