%  solns = get_generic_solns(filename,num_vars,varargin)
%
%  reads solutions from a bertini 1 standard output file.
%
%
% copyright 2015, 2016 Daniel Brake
% University of Notre Dame
% Applied and Computational Mathematics and Statistics
% danielthebrake@gmail.com
%
%  Bertini (TM) is a registered trademark.
%
%     This file is part of Bertini_tropical.
%
%     Bertini_tropical is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     Bertini_tropical is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


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
	pathid = fscanf(fid,'%i\n');
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
		solns(ii).soln = tempsoln;
	else
		
		if ~use_str
			solns(:,ii) = tempsoln;
		else
			solns{ii} = tempsoln;
		end
		
	end
	fgetl(fid); %burn a line
end


fclose(fid);
if use_vpa
	digits(max_digits);
end

end
