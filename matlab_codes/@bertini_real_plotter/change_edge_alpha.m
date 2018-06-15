% a callback function to be linked to a button in the bertini_real_plotter
% interface


% this is sooo repetitive with several other similar functions.  begs for
% abstraction
function br_plotter = change_edge_alpha(br_plotter,source,event)
			

[alpha_options, user_cancelled] = getopt(br_plotter);
			
if user_cancelled
	return;
end

new_raw_alpha = alpha_options.blocky;
new_sample_alpha = alpha_options.refinement;

if br_plotter.options.raw_triangulation_alpha ~= new_raw_alpha
	set(br_plotter.handles.faces(:),'EdgeAlpha',new_raw_alpha);
	br_plotter.options.raw_triangulation_alpha = new_raw_alpha;
end

if br_plotter.options.sample_triangulation_alpha ~= new_sample_alpha
set(br_plotter.handles.surface_samples(:),'EdgeAlpha',new_sample_alpha);
br_plotter.options.sample_triangulation_alpha = new_sample_alpha;
end
	
end









%http://www.mathworks.com/matlabcentral/fileexchange/25862-inputsdlg--enhanced-input-dialog-box--v2-0-5-


function [answer, cancelled] = getopt(br_plotter)

Title = 'Bertini_real Edge Alpha options';

%%%% SETTING DIALOG OPTIONS
% Options.WindowStyle = 'modal';
Options.Resize = 'on';
Options.Interpreter = 'none';
Options.CancelButton = 'on';
Options.ApplyButton = 'off';
Options.ButtonNames = {'Continue','Cancel'}; %<- default names, included here just for illustration


Prompt = {};
Formats = {};
default_answers = struct([]);


Prompt(1,:) = {'Raw edge alpha       ', 'blocky',[]};
Formats(1,1).type = 'edit';
Formats(1,1).format = 'float';
Formats(1,1).limits = [0 1]; % 9-digits (positive #)
Formats(1,1).size = 60;
default_answers(1).blocky = br_plotter.options.raw_triangulation_alpha;


Prompt(2,:) = {'Refined edge alpha                ', 'refinement',[]};
Formats(2,1).type = 'edit';
Formats(2,1).format = 'float';
Formats(2,1).limits = [0 1]; % 9-digits (positive #)
Formats(2,1).size = 60;
default_answers(1).refinement = br_plotter.options.sample_triangulation_alpha;


[answer,cancelled] = inputsdlg(Prompt,Title,Formats,default_answers,Options);


end

