% a callback function to be linked to a button in the bertini_real_plotter
% interface

function br_plotter = change_alpha(br_plotter,source,event)
			



[alpha_options, user_cancelled] = getopt(br_plotter);
			

if user_cancelled
	return;
end

new_face_alpha = alpha_options.blocky;
new_sample_alpha = alpha_options.refinement;

if br_plotter.options.face_alpha ~= new_face_alpha
set(br_plotter.handles.faces(:),'FaceAlpha',new_face_alpha);
br_plotter.options.face_alpha = new_face_alpha;
end

if br_plotter.options.sample_alpha ~= new_sample_alpha
set(br_plotter.handles.surface_samples(:),'FaceAlpha',new_sample_alpha);
br_plotter.options.sample_alpha = new_sample_alpha;
end

set(br_plotter.buttons.alpha,'String',sprintf('FaceAlpha'));
% set(br_plotter.buttons.alpha,'String',sprintf('FaceAlpha (%1.2f / %1.2f)',br_plotter.options.face_alpha,br_plotter.options.sample_alpha));


			
end









%http://www.mathworks.com/matlabcentral/fileexchange/25862-inputsdlg--enhanced-input-dialog-box--v2-0-5-


function [answer, cancelled] = getopt(br_plotter)

Title = 'Bertini_real FontSize options';

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


Prompt(1,:) = {'Raw face alpha       ', 'blocky',[]};
Formats(1,1).type = 'edit';
Formats(1,1).format = 'float';
Formats(1,1).limits = [0 1]; % 9-digits (positive #)
Formats(1,1).size = 60;
default_answers(1).blocky = br_plotter.options.face_alpha;


Prompt(2,:) = {'Refined face alpha                ', 'refinement',[]};
Formats(2,1).type = 'edit';
Formats(2,1).format = 'float';
Formats(2,1).limits = [0 1]; % 9-digits (positive #)
Formats(2,1).size = 60;
default_answers(1).refinement = br_plotter.options.sample_alpha;


[answer,cancelled] = inputsdlg(Prompt,Title,Formats,default_answers,Options);


end

