% check whether an edge touches the rendered faces
function val = edge_touches_faces(br_plotter, edge_index, curve)
	
	val = false;
	e = curve.edges(edge_index,:); % the current edge
	m = e(2); % the midpoint index
	
	for ii = 1:length(br_plotter.options.which_faces)
		face_ind = br_plotter.options.which_faces(ii);

		this_face_touches = false;

		
		f = br_plotter.BRinfo.faces(face_ind); % unpack the current face
		
		
		if m == f.midpoint 
			this_face_touches = true;
			where = 'midpoint of face';
		end

		left_critslice = br_plotter.BRinfo.critpoint_slices{f.midslice_index};
		for jj = 1:f.num_left
			if m == left_critslice.edges(f.left(jj),2)
				this_face_touches = true;
				where = 'left critslice of face';
			end
		end

		right_critslice = br_plotter.BRinfo.critpoint_slices{f.midslice_index+1};
		for jj = 1:f.num_right
			if m == right_critslice.edges(f.right(jj),2)
				this_face_touches = true;
				where = 'right critslice of face';
			end
		end

		if this_face_touches
			val = this_face_touches;
			break;
		end
	end



end
