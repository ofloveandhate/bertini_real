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
		end

		left_critslice = br_plotter.BRinfo.critpoint_slices{f.midslice_index};
		for jj = 1:f.num_left
			if m == left_critslice.edges(f.left(jj),2)
				this_face_touches = true;
			end
		end

		right_critslice = br_plotter.BRinfo.critpoint_slices{f.midslice_index+1};
		for jj = 1:f.num_right
			if m == right_critslice.edges(f.right(jj),2)
				this_face_touches = true;
			end
		end
	
		% check the top and bottom
		
		curve_top = curve_with_inputname(br_plotter.BRinfo, f.system_top);
		if m == curve_top.edges(f.top,2)
			this_face_touches = true;
		end
		
		curve_bottom = curve_with_inputname(br_plotter.BRinfo, f.system_bottom);
		if m == curve_bottom.edges(f.bottom,2)
			this_face_touches = true;
		end
		
		if this_face_touches
			val = this_face_touches;
			break;
		end
	end



end


function curve = curve_with_inputname(BRinfo, name)


if strcmp(BRinfo.crit_curve.inputfilename,name)
	curve = BRinfo.crit_curve;
	return
elseif strcmp(BRinfo.sphere_curve.inputfilename,name)
	curve = BRinfo.sphere_curve;
	return
else
	for ii = 1:BRinfo.num_singular_curves
		if strcmp(BRinfo.singular_curves{ii}.inputfilename,name)
			curve = BRinfo.singular_curves{ii};
			return
		end
	end
end


error('unable to find crit-like curve of name %s in BRinfo',name);

end



