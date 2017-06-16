%this function generates a filename, based on highest number already in
%directory.  

%could potentially be used elsewhere for autogeneration of
%names.

% daniel brake
% colorado state university, north carolina state university, notre dame, 
% university of wisconsin eau claire
% mathematics and applied mathematics
% 2013-17
%
% danielthebrake@gmail.com
% brakeda@uwec.edu



function newname = increment_name(basename)

file_list = dir([basename '_*']);
if ~isempty(file_list)
	filenumbers = zeros(1,length(file_list));
	for ii = 1:length(file_list)
		filenumbers(ii) = extract_number(file_list(ii).name, basename);
	end
else
	filenumbers = 0;
end

newname = sprintf('%s_%i',basename,max(filenumbers)+1);


end

function n = extract_number(name, basename)
	
if has_extension(name)
	r = [length(basename)+2 length(name)-4];
else
	r = [length(basename)+2 length(name)];
end
n = str2double(name(r(1):r(2)));

if isnan(n)
	n = -1;
end


end


function result =  has_extension(name)

if strcmp((name(end-3)), '.')
	result = true;
else
	result = false;
end

end

