function handles = plot_sampler_data(ind, vertices,sampler_data,colors)
%
global plot_params


curr_axis = plot_params.axes.main;
nondegen_edge_ind = 1;

handles = [];
% handles = zeros(sum( length(sampler_data.sample_sizes>3)),1);
for ii = 1:length(sampler_data.edge)
    if sampler_data.sample_sizes(ii)<=3
       continue; 
    end
    
	plotme = zeros(sampler_data.sample_sizes(ii),length(ind));
	
	for jj = 1:sampler_data.sample_sizes(ii)
		plotme(jj,:) = vertices(sampler_data.edge(ii).samples(jj)+1).point(ind);
	end
	plotme = real(plotme);
	h = main_plot_function(plotme,1:length(ind),curr_axis);

	colors
	nondegen_edge_ind
	
	set(h,'Color',colors(nondegen_edge_ind,:));
	set(h,'LineWidth',3);
    
    
    handles(nondegen_edge_ind) = h;
	
	
	nondegen_edge_ind = nondegen_edge_ind+1;
	
end

end