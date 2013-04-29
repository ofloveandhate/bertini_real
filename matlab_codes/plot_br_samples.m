function plot_br_samples(edges,BRinfo)

%% initialize
	fontsize = 18;
	window = [20 20 1024 768];
	format = 'eps';
	format_flag = 'epsc2';
	
	
	
	
	
if nargin==0
	load('edges.mat');
end

close all

colors = jet(BRinfo.num_edges);


plot_params.fontsize = fontsize;
plot_params.window = window;
plot_params.format = format;
plot_params.format_flag = format_flag;




%% set up the plotting indices -- which data to plot
if (BRinfo.num_variables-1)==2
	plot_indices = [1 2];
elseif (BRinfo.num_variables-1)==3
	plot_indices = [1 2 3];
else
		%first scan for zero columns.
	tmpdata = zeros(sum(BRinfo.sample_sizes),BRinfo.num_variables-1);
	counter = 1;
	for ii = 1:BRinfo.num_edges
		for jj = 1:BRinfo.sample_sizes(ii)
			tmpdata(counter,:) = edges(ii).samples(jj).soln;
			counter = counter+1;
		end
	end

	indices_of_nonzero_cols = find_zero_vars(tmpdata);

	if length(indices_of_nonzero_cols)<4
		plot_indices = indices_of_nonzero_cols;
	else
		plot_indices = get_user_indices(indices_of_nonzero_cols,BRinfo);
	end

end


%% actually plot
for ii = 1:BRinfo.num_edges
	plotme = zeros(BRinfo.sample_sizes(ii),BRinfo.num_variables-1);
	
	for jj = 1:BRinfo.sample_sizes(ii)
		plotme(jj,:) = edges(ii).samples(jj).soln;
	end
	plotme = real(plotme);
	

	h = main_plot_function(plotme,plot_indices,BRinfo,plot_params);
	
	set(h,'Color',colors(ii,:));
	set(h,'LineWidth',2);
	hold on
end
% 	axis square;



render_into_file('br_samples',plot_params);


end%re: function


function h = main_plot_function(data,ind,infostruct,plot_params)

switch length(ind)
	case {2}
		h = plot(data(:,ind(1)),data(:,ind(2)),'-');
		xlabel(infostruct.var_names{ind(1)}, 'interpreter', 'none','FontSize',plot_params.fontsize)
		ylabel(infostruct.var_names{ind(2)}, 'interpreter', 'none','FontSize',plot_params.fontsize)
		
	case {3}
		h = plot3(data(:,ind(1)),data(:,ind(2)),data(:,ind(3)),'-');
		xlabel(infostruct.var_names{ind(1)}, 'interpreter', 'none','FontSize',plot_params.fontsize)
		ylabel(infostruct.var_names{ind(2)}, 'interpreter', 'none','FontSize',plot_params.fontsize)
		zlabel(infostruct.var_names{ind(3)}, 'interpreter', 'none','FontSize',plot_params.fontsize)
	otherwise
			
end
	
end





function ind = get_user_indices(suggested_indices, infostruct)

'getting user indices'

'code me!'

ind = [];

end

function ind = find_zero_vars(data,dim)
zerothresh = 1e-8;
if nargin < 2 
	dim = 1;
end


biggies = max(abs(data),[],dim);



ind = find(biggies>zerothresh);


	

end