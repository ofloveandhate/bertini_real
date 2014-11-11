function br_trig

global plot_params


prev_filenames = dir('BRinfo*.mat');
max_found = -1;

for ii = 1:length(prev_filenames)
	curr_name = prev_filenames(ii).name;
	curr_num = str2num(curr_name(7:end-4));
	if max_found < curr_num
		max_found = curr_num;
	end

end
filename = ['BRinfo' num2str(max_found) '.mat'];
load(filename);
[containing, here, ~] = fileparts(pwd);
plot_params.basename = here;


if ~isvarname('sampler_data')
	return;
end


colors = hsv(BRinfo.num_edges);
for ii = 1:BRinfo.num_edges
	curr_vars = zeros(sampler_data.sample_sizes(ii),(BRinfo.num_variables-1)/2);
	for jj = 1:sampler_data.sample_sizes(ii)
		for kk = 1:2:BRinfo.num_variables-1
			
			curr_vars(jj,(kk+1)/2) = atan2( real(sampler_data.edge(ii).samples(jj).soln(kk)), real(sampler_data.edge(ii).samples(jj).soln(kk+1)) );
		end
	end
% 	h = main_plot_function(curr_vars,[1 2],gca);
	h = main_plot_function(mod(curr_vars,2*pi),[1 2],gca);
	set(h,'Color',colors(ii,:));
	%plot here
	hold on
end


axis([0 2*pi 0 2*pi]);
hold off;

end