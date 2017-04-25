function save_all_fvs(fvs,root_filename)

options.remove_duplicates = false;
options.align = false;
for ii = 1:length(fvs)
	options.filename = sprintf('%s_%i.stl',root_filename, ii);
	fv_to_stl(fvs{ii}, options);
end

end