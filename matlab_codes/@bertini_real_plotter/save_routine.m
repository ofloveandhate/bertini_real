function br_plotter = save_routine(br_plotter,varargin) % additional arguments go afer the first two, which must be put here ,even if they are not used.
%perhaps get info from other calls here?

f = fieldnames(br_plotter.panels);

for ii = 1:length(f)
	set( findall(br_plotter.panels.(f{ii}), '-property', 'visible'), 'visible', 'off')
	set(br_plotter.panels.(f{ii}),'visible','off');
	
end


if or(br_plotter.switches.display_faces == 1,br_plotter.switches.display_face_samples == 1)
    br_plotter.format = 'png';
    br_plotter.format_flag = 'png';
else
    br_plotter.format = 'eps';
    br_plotter.format_flag = 'epsc2';
end

render_into_file(br_plotter);

for ii = 1:length(f)
	set(br_plotter.panels.(f{ii}),'visible','on');
	set( findall(br_plotter.panels.(f{ii}), '-property', 'visible'), 'visible', 'on')
end

end