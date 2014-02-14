%this function generates a filename, based on highest number already in
%directory.  

%could potentially be used elsewhere for autogeneration of
%names.

% daniel brake
% colorado state university
% mathematics
% 2013
% danielthebrake@gmail.com



function newname = increment_name(picname)

file_list = dir([picname '_*']);
if ~isempty(file_list)
	filenumbers = zeros(1,length(file_list));
	for ii = 1:length(file_list)

		filenumbers(ii) = str2double(file_list(ii).name(length(picname)+2:end-4));%subtract 4 to remove extension
	end
else
	filenumbers = 0;
end

newname = sprintf('%s_%i',picname,max(filenumbers)+1);

display(sprintf('saving with filename: %s',newname));
end