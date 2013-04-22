function plot_br_samples(edges,BRinfo)

if nargin==0
	load('edges.mat');
end

close all

colors = jet(BRinfo.num_edges);

for ii = 1:BRinfo.num_edges
	plotme = zeros(BRinfo.sample_sizes(ii),BRinfo.num_variables-1);
	
	
	for jj = 1:BRinfo.sample_sizes(ii)
		plotme(jj,:) = edges(ii).samples(jj).soln;
	end
	
	h = plot(plotme(:,1),plotme(:,2));
	set(h,'Color',colors(ii,:));
	hold on
end
	axis square;


render_into_file('br_samples');


end%re: function