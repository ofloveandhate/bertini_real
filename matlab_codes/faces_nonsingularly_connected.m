% gets the indices of all faces nonsingularly connected to a gived face.
% 1-based indices is the convention here
%
% the computed set will include the seed index
%


function [connected, unconnected] = faces_nonsingularly_connected(BRsurf, seed_index)

	check_data(BRsurf);



	new_indices = seed_index;
	connected = [];

	while ~isempty(new_indices)
		connected = [connected new_indices];
		[new_indices] = find_connected_faces(BRsurf, connected);
	end
	connected = sort(connected);
	unconnected = setdiff(1:BRsurf.num_faces, connected);
end


function [new_indices] = find_connected_faces(BRsurf, current)
	% perform a double loop over the faces
	new_indices = [];

	unexamined_indices = 1:BRsurf.num_faces;
	unexamined_indices(current) = []; % delete the faces we already know connect. 
	
	for ii=1:length(current)
		c = current(ii);
		f = BRsurf.faces(c); % unpack the current face
		deleteme = [];
		for jj = 1:length(unexamined_indices)
			d = unexamined_indices(jj);
			g = BRsurf.faces(d); % unpack the examined face
			if faces_nonsingularly_connect(BRsurf, f, g)
				new_indices(end+1) = d;
				deleteme(end+1) = d;
			end
		end
		unexamined_indices = setdiff(unexamined_indices,deleteme);
	end
end


%queries whether faces f and g nonsingularly connect
function val = faces_nonsingularly_connect(BRsurf, f, g)
	val = false; %assume no

	if cannot_possibly_meet(BRsurf, f, g);
		return;
	elseif meet_at_left(BRsurf,f,g)
		val = true;
	elseif meet_at_right(BRsurf,f,g)
		val = true;
	elseif meet_at_top(BRsurf,f,g)
		val = true;
	elseif meet_at_bottom(BRsurf,f,g)
		val = true;
	end

end

function val = cannot_possibly_meet(BRsurf, f,g)
	val = false;
	if abs(f.midslice_index - g.midslice_index)>=2
		val = true;
	end
end

function val = meet_at_left(BRsurf, f, g)
	val = false;
	for ii = 1:f.num_left
		e = BRsurf.critpoint_slices{f.midslice_index}.edges(f.left(ii),:);
		a = e(2);
		for jj = 1:g.num_left
			E = BRsurf.critpoint_slices{g.midslice_index}.edges(g.left(jj),:);
			b = E(2);
			if a == b && ~is_degenerate (e) && ~is_degenerate (E)
				val = true;
				return
			end
		end

		for jj = 1:g.num_right
			E = BRsurf.critpoint_slices{g.midslice_index+1}.edges(g.right(jj),:);
			b = E(2);
			if a == b && ~is_degenerate (e) && ~is_degenerate (E)
				val = true;
				return
			end
		end
	end
end


function val = meet_at_right(BRsurf,f, g)
	val = false;
	for ii = 1:f.num_right
		e = BRsurf.critpoint_slices{f.midslice_index+1}.edges(f.right(ii),:);
		a = e(2);
		for jj = 1:g.num_left
			E = BRsurf.critpoint_slices{g.midslice_index}.edges(g.left(jj),:);
			b = E(2);
			if a == b && ~is_degenerate (e) && ~is_degenerate (E)
				val = true;
				return
			end
		end

		for jj = 1:g.num_right
			E = BRsurf.critpoint_slices{g.midslice_index+1}.edges(g.right(jj),:);
			b = E(2);
			if a == b && ~is_degenerate (e) && ~is_degenerate (E)
				val = true;
				return
			end
		end
	end
end

function val = meet_at_top(BRsurf,f, g)
	val = false;
	if strcmp(f.system_top(1:15),'input_singcurve')
		return % cannot meet singularly, because edge is singular
	elseif f.midslice_index~=g.midslice_index % at least they are in the same interval
		return
	end

	if strcmp(f.system_top, g.system_top)
		if strcmp(BRsurf.crit_curve.inputfilename, f.system_top);
			if f.top==g.top
				val = true;
				return;
			end
		end
	end

	if strcmp(f.system_top, g.system_bottom)
		if strcmp(BRsurf.crit_curve.inputfilename, f.system_top);
			if f.top==g.bottom
				val = true;
				return;
			end
		end
	end
end



function val = meet_at_bottom(BRsurf,f, g)
	val = false;
	if strcmp(f.system_bottom(1:15),'input_singcurve')
		return % cannot meet singularly, because edge is singular
	elseif f.midslice_index~=g.midslice_index % at least they are in the same interval
		return
	end

	if strcmp(f.system_bottom, g.system_top)
		if strcmp(BRsurf.crit_curve.inputfilename, f.system_bottom);
			if f.bottom==g.top
				val = true;
				return;
			end
		end
	end

	if strcmp(f.system_bottom, g.system_bottom)
		if strcmp(BRsurf.crit_curve.inputfilename, f.system_bottom);
			if f.bottom==g.top
				val = true;
				return;
			end
		end
	end
end

function val = is_degenerate(e)
	val = or(e(1)==e(2), e(2)==e(3));
end

function check_data(BRsurf)
	if BRsurf.dimension ~= 2
		error('this function designed to work on surfaces decomposed with bertini_real.  your object has dimension %i', BRsurf.dimension);
	end
end