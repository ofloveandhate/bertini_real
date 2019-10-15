% [] = save_all_fvs(fvs,root_filename)
%
% saves a cell array of faces-vertices structs to stl's.
%
% requires two arguments.  
%
% fvs is a cell array of fv
% root_filename is a string


function save_all_fvs(fvs,root_filename)

options.remove_duplicates = false;
options.align = false;
options.save_mat = false;

for ii = 1:length(fvs)
	options.filename = sprintf('%s_%i.stl',root_filename, ii);
	fv2stl(fvs{ii}, options);
end

end