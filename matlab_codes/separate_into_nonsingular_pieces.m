function [indices] = separate_into_nonsingular_pieces(BRsurf)

check_data(BRsurf);

indices = cell(0);

connected = [];

unconnected_this = 1;

while ~isempty(unconnected_this)
	seed = unconnected_this(1);
	[connected_this, unconnected_this] = faces_nonsingularly_connected(BRsurf, seed);
	indices{end+1} = connected_this;
	connected = [connected connected_this];
	unconnected_this = setdiff(unconnected_this, connected);

end

end




function check_data(BRsurf)
	if BRsurf.dimension ~= 2
		error('this function designed to work on surfaces decomposed with bertini_real.  your object has dimension %i', BRsurf.dimension);
	end
end