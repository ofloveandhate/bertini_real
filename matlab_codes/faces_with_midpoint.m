% function [indices] = faces_with_midpoint(BRinfo)
% 
% computes which faces on a surface have midpoints satisfying a
% predicate


function [indices] = faces_with_midpoint(BRsurf, predicate)

check_data(BRsurf);
check_pred(predicate);

indices = [];

for ii=1:BRsurf.num_faces
	if predicate(BRsurf.vertices(BRsurf.faces(ii).midpoint).point)
		indices(end+1) = ii;
	end
end


end



function check_data(BRsurf)
	if BRsurf.dimension ~= 2
		error('this function designed to work on surfaces decomposed with bertini_real.  your object has dimension %i', BRsurf.dimension);
	end
end


function check_pred(predicate)
	if ~isa(predicate,'function_handle')
		error('second input argument must be an evaluable handle returning a bool');
	end
end