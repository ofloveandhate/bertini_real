%
%get_path_linprod()
%
%reads in the file 'pathtrack_linprod', output by a suitably compiled
%version of bertini_real, and stores the path data in a .mat file.
%
% daniel brake
% colorado state university
% 2013
% danielthebrake@gmail.com
%


function get_path_solver(solver)


if nargin<1
	solver = 'linprod';
end

fid = fopen(sprintf('pathtrack_%s',solver),'r');



%read in the top two lines of header info.

tempstr = fgetl(fid);
parsed=strread(tempstr,'%s','delimiter',' ');
num_variables = str2num(parsed{1});

variable_names = cell(num_variables,1);
for ii = 1:num_variables
	variable_names{ii} = parsed{ii+1};
end

tempstr = fgetl(fid);
parsed=strread(tempstr,'%s','delimiter',' ');
num_pts = str2num(parsed{1});
MPType = str2num(parsed{2});
odepredictor = str2num(parsed{3});
endgame = str2num(parsed{4});

 path_ = zeros(1,2);
 
 tmp = struct('path',path_,'num_steps',0);
 data = repmat(tmp,[1,num_pts]);
 data(1).num_pts = num_pts;
 data(1).num_variables = num_variables-1;
 data(1).variable_names = variable_names;
 data(1).MPType = MPType;
 data(1).odepredictor = odepredictor;
 data(1).endgame = endgame;
 data(1).solver = solver;
 
for ii = 1:num_pts
	display(sprintf('reading path %i',ii));
	data(ii).path = zeros(3,2*num_variables);
	counter = 0;
	while 1
		counter = counter+1;
		tmp = fscanf(fid,'%f %f',[1 2*num_variables]);

		if tmp(1)==-100 % the break flag
			data(ii).retval = fscanf(fid,'%i\n\n',[1 1]);
			data(ii).num_steps = tmp(2)+1;
			break
		else %store the data
			data(ii).path(counter,:) = tmp;
		end
	end
	
end

fclose(fid);



save(sprintf('pathdata_%s_ode%i_eg%i_mp%i.mat',solver,odepredictor,endgame,MPType),'data','-v6');

end%re: function