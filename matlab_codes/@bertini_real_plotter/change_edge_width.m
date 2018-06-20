% a callback function to be linked to a button in the bertini_real_plotter
% interface


% this is sooo repetitive with several other similar functions.  begs for
% abstraction
function br_plotter = change_edge_width(br_plotter,source,event)
			

[alpha_options, user_cancelled] = getopt(br_plotter);
			
if user_cancelled
	return;
end

new_val = alpha_options.all;

if br_plotter.options.edge_width ~= new_val
	set(br_plotter.handles.faces(:),'LineWidth',new_val); 
	set(br_plotter.handles.faces.samples,'LineWidth',new_val);
	br_plotter.options.edge_width = new_val;
end
	
end









%http://www.mathworks.com/matlabcentral/fileexchange/25862-inputsdlg--enhanced-input-dialog-box--v2-0-5-


function [answer, cancelled] = getopt(br_plotter)

Title = 'Bertini_real Edge Width options';

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


Prompt(1,:) = {'Edge width       ', 'all',[]};
Formats(1,1).type = 'edit';
Formats(1,1).format = 'float';
Formats(1,1).limits = [0 20]; % 9-digits (positive #)
Formats(1,1).size = 60;
default_answers(1).all = br_plotter.options.edge_width;





[answer,cancelled] = inputsdlg(Prompt,Title,Formats,default_answers,Options);


end

