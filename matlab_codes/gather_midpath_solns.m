%  solns = gather_midpath_solns(filename,num_vars,varargin)
%
%  reads solutions from a bertini 1 midpath output file.
%
%
% copyright 2017 Dani Brake
% University of Notre Dame
% Applied and Computational Mathematics and Statistics
% danielthebrake@gmail.com
%
%  Bertini (TM) is a registered trademark.
%

function solns = gather_midpath_solns(filename,num_vars,varargin)

if ~ischar(filename)
	display('input not a string');
	return;
elseif isempty(dir(filename))
	error('folder has no file with the name %s',filename);
end

use_vpa = false;
use_str = false;

struct_output = false;

converter = @str2double;

for ii = 1:length(varargin)
	switch(varargin{ii})
		case 'struct'
			if varargin{ii+1}
				struct_output = true;
			end
			ii = ii+1;
		case 'vpa'
			use_vpa = true;
			use_str = false;
			converter = @sym;
		case 'str'
			use_str = true;
			use_vpa = false;
			converter = @(x) x;
		otherwise
			error('unknown option %s',varargin{ii});
	end
end


fid = fopen(filename,'r');	



num_solns = 1;

if use_vpa
	tempsoln = sym(zeros(num_vars,1));
	if struct_output
		solns = repmat(struct('soln',[]),[1,num_solns]);
	else
		solns = vpa(zeros(num_vars,num_solns));
	end
elseif use_str
	tempsoln = cell(num_vars,1);
	solns = cell(1,num_solns);
else
	tempsoln = zeros(num_vars,1);
	if struct_output
		solns = repmat(struct('soln',[]),[1,num_solns]);
	else
		solns = zeros(num_vars,num_solns);
	end
end

max_digits = 0;
ii = 0;
while ~feof(fid)
	ii = ii+1;
	pathid = fscanf(fid,'%i\n',[1 1])+1;
	for jj = 1:num_vars
		tempstr = fgetl(fid);
		parsed=strread(tempstr,'%s','delimiter',' ');
		if length(parsed{1})>max_digits
			max_digits = length(parsed{1});
		end
		if use_str
			tempsoln{jj} = [parsed{1} ' ' parsed{2}];
		else
			tempsoln(jj) = converter(parsed{1})+1i*converter(parsed{2});
		end
	end
	
	if struct_output
		solns(pathid).soln = tempsoln;
	else
		
		if ~use_str
			solns(:,pathid) = tempsoln;
		else
			solns{pathid} = tempsoln;
		end
		
	end
	fgetl(fid); %burn a line
end


fclose(fid);
if use_vpa
	digits(max_digits);
end

end
