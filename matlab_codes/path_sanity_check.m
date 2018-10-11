function path_sanity_check()


loc.gather = which('gather_br_samples');

if isempty(loc.gather)
	error('something crazy is going on, you don''t even have gather_br_samples on your path.  please add `path/to/bertini_real/matlab_codes` to your matlab path and try again')
end

[codepath,~] = fileparts(loc.gather);

brakelabstuff = dir(sprintf('%s/%s',codepath,'brakelab'));
if isempty(brakelabstuff)
	warning(failure_to_find_message());
end


loc.render = which('render_into_file');
if isempty(loc.render)
	warning('please add `%s/brakelab/rendering` to your matlab path for the `save` button in the plotter to work',codepath)
end

loc.dehom = which('dehomogenize');
if isempty(loc.dehom)
	warning('please add `%s/brakelab/bertini1` to your matlab path for `gather_br_samples` to work',codepath)
end





end


function msg = failure_to_find_message()

msg = 'when doing a sanity check on the matlab code, `brakelab` was not found, but its a submodule of bertini_real... from the bertini_real code folder, simply run `git submodule update --init`.  then, add the three subfolders to the matlab path';
end