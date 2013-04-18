function solns = get_generic_solns(filename,num_vars)

if ~ischar(filename)
	display('input not a string');
	return;
else if isempty(dir(filename))
	display(sprintf('folder has no file with the name %s',filename));
	return;
	end


fid = fopen(filename,'r');	

tempstr = fgetl(fid);
parsed=strread(tempstr,'%s','delimiter',' ');
num_solns = str2num(parsed{1});

fgetl(fid); %burn a line
solns = repmat(struct('soln',[]),[1,num_solns]);

for ii = 1:num_solns
	tempsoln = zeros(num_vars,1);
	for jj = 1:num_vars
		tempstr = fgetl(fid);
		parsed=strread(tempstr,'%s','delimiter',' ');
		tempsoln(jj) = str2num(parsed{1})+1i*str2num(parsed{2});
	end
	solns(ii).soln = tempsoln;
	fgetl(fid); %burn a line
end


fclose(fid);


end