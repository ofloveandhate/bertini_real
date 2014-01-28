%curve specific code

function curve_plot(sampler_data, BRinfo,plot_indices)
global plot_params
%

num_non_degen = 0;
for ii = 1:BRinfo.num_edges
	
	left_plot = BRinfo.vertices(BRinfo.edges(ii,1)).point;
	mid_plot = BRinfo.vertices(BRinfo.edges(ii,2)).point;
	right_plot = BRinfo.vertices(BRinfo.edges(ii,3)).point;
	
    if and(left_plot==mid_plot, mid_plot==right_plot)
       continue; 
    else
        num_non_degen = num_non_degen+1;
    end
    
end

colors = jet(num_non_degen);








create_axes_curve(BRinfo);

label_axes(plot_indices,BRinfo,plot_params.axes.vertices);

% actually plot
plot_projection(BRinfo,plot_indices);



plot_vertices(plot_indices, BRinfo);

if ~isempty(sampler_data)
	plot_sampler_data(plot_indices, BRinfo,sampler_data,colors);
else
	plot_edge_points(plot_indices,BRinfo,0.8*colors);
end




sync_axes();

% colorbar
make_discrete_colorbar(1,num_non_degen,plot_params.fontsize,'edge number');
end













function plot_edge_points(ind, BRinfo,colors)
global plot_params

%
curr_axis = plot_params.axes.edges;

nondegen_edge_ind = 1;
for ii = 1:BRinfo.num_edges
	
	left_plot = BRinfo.vertices(BRinfo.edges(ii,1)).point;
	mid_plot = BRinfo.vertices(BRinfo.edges(ii,2)).point;
	right_plot = BRinfo.vertices(BRinfo.edges(ii,3)).point;
	
    if and(left_plot==mid_plot, mid_plot==right_plot)
       continue; 
    end
	switch length(ind)
		case {2}
			plotme = [left_plot(ind(1)) left_plot(ind(2));...
				mid_plot(ind(1)) mid_plot(ind(2));...
				right_plot(ind(1)) right_plot(ind(2))];
			plotme = real(plotme);
			
			h = main_plot_function(plotme,[1 2], curr_axis);
			
		case {3}
		
			plotme = [left_plot(ind(1)) left_plot(ind(2)) left_plot(ind(3));...
				mid_plot(ind(1)) mid_plot(ind(2)) mid_plot(ind(3));...
				right_plot(ind(1)) right_plot(ind(2)) right_plot(ind(3))];
			plotme = real(plotme);
			h = main_plot_function(plotme,[1 2 3], curr_axis);
% 			h = plot3(plotme(:,1),plotme(:,2),plotme(:,3),'--');
	end
	set(h,'Color',colors(nondegen_edge_ind,:),'LineStyle','--','LineWidth',2);
    nondegen_edge_ind = nondegen_edge_ind+1;
end


end

function plot_sampler_data(ind, BRinfo,sampler_data,colors)
%
global plot_params
curr_axis = plot_params.axes.sampler;

nondegen_edge_ind = 1;

for ii = 1:BRinfo.num_edges
    
    left_plot = BRinfo.vertices(BRinfo.edges(ii,1)).point;
	mid_plot = BRinfo.vertices(BRinfo.edges(ii,2)).point;
	right_plot = BRinfo.vertices(BRinfo.edges(ii,3)).point;
    if and(left_plot==mid_plot, mid_plot==right_plot)
       continue; 
    end
    
	plotme = zeros(sampler_data.sample_sizes(ii),BRinfo.num_variables-1);
	
	for jj = 1:sampler_data.sample_sizes(ii)
		plotme(jj,:) = sampler_data.edge(ii).samples(jj).soln;
	end
	plotme = real(plotme);
	h = main_plot_function(plotme,ind,curr_axis);
	set(h,'Color',colors(nondegen_edge_ind,:));
	set(h,'LineWidth',2);
    nondegen_edge_ind = nondegen_edge_ind+1;
end

end






function create_axes_curve(BRinfo)
global plot_params



new_axes = axes('Parent',plot_params.figures.main_plot);
delete( get(new_axes,'Children') );
set(new_axes,'visible','on');
hold(new_axes,'on');

set(new_axes,'XLim',[BRinfo.center(plot_params.ind(1))-1.1*BRinfo.radius BRinfo.center(plot_params.ind(1))+1.1*BRinfo.radius]);
set(new_axes,'YLim',[BRinfo.center(plot_params.ind(2))-1.1*BRinfo.radius BRinfo.center(plot_params.ind(2))+1.1*BRinfo.radius]);

if length(plot_params.ind)==3
    set(new_axes,'ZLim',[BRinfo.center(plot_params.ind(3))-1.1*BRinfo.radius BRinfo.center(plot_params.ind(3))+1.1*BRinfo.radius]);
end

plot_params.axes.vertices = new_axes;


new_axes = copyobj(plot_params.axes.vertices,plot_params.figures.main_plot);
delete( get(new_axes,'Children') );
set(new_axes,'visible','off');
hold(new_axes,'on')
plot_params.axes.projection = new_axes;


new_axes = copyobj(plot_params.axes.vertices,plot_params.figures.main_plot);
delete( get(new_axes,'Children') );
set(new_axes,'visible','off');
hold(new_axes,'on');
plot_params.axes.edges = new_axes;

new_axes = copyobj(plot_params.axes.vertices,plot_params.figures.main_plot);
delete( get(new_axes,'Children') );
set(new_axes,'visible','off');
hold(new_axes,'on');
plot_params.axes.sampler = new_axes;

end





