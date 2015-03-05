function solns = get_generic_solns(filename,num_vars,varargin)

if ~ischar(filename)
	display('input not a string');
	return;
elseif isempty(dir(filename))
	error(sprintf('folder has no file with the name %s',filename));
end


struct_output = false;

for ii = 1:length(varargin)
	switch(varargin{ii})
		case 'struct'
			if varargin{ii+1}
				struct_output = true;
			end
			
	end
end

fid = fopen(filename,'r');	

tempstr = fgetl(fid);
parsed=strread(tempstr,'%s','delimiter',' ');
num_solns = str2num(parsed{1});

fgetl(fid); %burn a line

if struct_output
	solns = repmat(struct('soln',[]),[1,num_solns]);
else
	solns = zeros(num_vars,num_solns);
end


for ii = 1:num_solns
	tempsoln = zeros(num_vars,1);
	for jj = 1:num_vars
		tempstr = fgetl(fid);
		parsed=strread(tempstr,'%s','delimiter',' ');
		tempsoln(jj) = str2num(parsed{1})+1i*str2num(parsed{2});
	end
	
	if struct_output
		solns(ii).soln = tempsoln;
	else
		solns(:,ii) = tempsoln;
	end
	fgetl(fid); %burn a line
end


fclose(fid);


end
