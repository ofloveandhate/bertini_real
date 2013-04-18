function r = read_random_matrix()

fid = fopen('Rand_matrix','r');

rows = fscanf(fid,'%i\n',[1 1]);
cols = fscanf(fid,'%i\n',[1 1]);
precision = fscanf(fid,'%i\n',[1 1]);

r = zeros(rows,cols);

for ii = 1:rows
	tempstr = fgetl(fid);
	parsed=strread(tempstr,'%s','delimiter',' ');
	for jj = 1:cols
		a = parsed{2*jj-1};
		b = parsed{2*jj};
		r(ii,jj) = str2num(a(2:end))+1i*str2num(b(1:end-1));
	end
end

fclose(fid);

end