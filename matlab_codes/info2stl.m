% [] = info2stl(BRinfo, options)
%
% a function for extracting the faces and vertices from
% a BRinfo object, and saving the result as an stl
% possibly for 3d printing, or whatever; i don't get to tell
% you what to do with your stl's
% 
% there's an internal bit you can flip to save the raw faces...

function [] = info2stl(BRinfo, varargin)


options = process_options(BRinfo, varargin);

fv = info2fv(BRinfo, options.which_faces, options.smooth_faces);
fv2stl(fv,options);


end % the function



function options = process_options(BRinfo, varargin)

[options.containing, options.filename, ~] = fileparts(pwd);

options.which_faces = [];
if isempty(options.which_faces)
	options.which_faces = 1:BRinfo.num_faces;
elseif max(options.which_faces) > BRinfo.num_faces
	error('trying to plot faces which don''t exist. requested index %i > num faces %i',...
		max(options.which_faces),...
		BRinfo.num_faces);
end
				

options.remove_duplicates = 0;
options.align = 0;
options.save_mat = 0;
options.smooth_faces = 0;
end
