function gather_br_samples()


% 19
% output_comp0_curve
% 2

fid = fopen('Dir_Name','r');
nchars = fscanf(fid,'%i\n',[1 1]);
dirname = fscanf(fid,'%s\n',[1 1])
nchars = fscanf(fid,'%i\n',[1 1]);

fclose(fid);




fid = fopen([dirname '/' 'E.edge'],'r');
num_variables = fscanf(fid,'%i\n',[1 1]);
fclose(fid);


filename = [dirname '/' 'samp.dat'];




fid = fopen(filename,'r');

num_edges = fscanf(fid,'%i\n\n',[1 1]);


% sample = zeros(num_variables,1);
% prototype = struct('projval',[],'samples',sample);


dehom_me = zeros(num_variables,1);
for ii = 1:num_edges
	num_samples = fscanf(fid,'%i\n\n',[1 1]);
	sample_sizes(ii) = num_samples;
	for jj=1:num_samples
		tmpsoln = zeros(num_variables-1,1);
		
		
		tmp_proj = fscanf(fid,'%f %f\n',[1 2]);
		for kk = 1:num_variables
		
			tmp = fscanf(fid,'%f %f\n',[1 2]);
			dehom_me(kk) = tmp(1)+1i*tmp(2);
		end
		
		for kk = 1:num_variables-1
			tmpsoln(kk) = dehom_me(kk+1)/dehom_me(1);
		end
		edges(ii).samples(jj).soln = tmpsoln;
		edges(ii).samples(jj).proj = tmp_proj;
	end
		
	
end


fclose(fid);


BRinfo.num_edges = num_edges;
BRinfo.num_variables = num_variables;
BRinfo.dir_name = dirname;
BRinfo.sample_sizes = sample_sizes;



save('edges.mat','edges','BRinfo');
end%re: function


