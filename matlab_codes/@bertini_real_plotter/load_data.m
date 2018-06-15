% uses internally set variable 'filename' to load a .mat file
% containing data gathered previously.
function load_data(br_plotter)
	if isempty(br_plotter.filename)
		error('unset filename in br_plotter object');
	end

	if isempty(dir(br_plotter.filename))
		if isempty(dir([br_plotter.filename '.mat']))
			error('nexists file with name ''%s''',br_plotter.filename);
		else
			br_plotter.filename = [br_plotter.filename '.mat'];
		end
	end

	file_variables = whos('-file',br_plotter.filename);

	if ismember('BRinfo', {file_variables.name})
		temp = load(br_plotter.filename);
		br_plotter.BRinfo = temp.BRinfo;
	else
		error('file ''%s'' does not contain variable ''BRinfo''',br_plotter.filename);
	end

end


