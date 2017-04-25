% fv = info2fv(BRinfo, which_faces)
%
% extract the faces from a bertini_real output.
%
% if which_faces is empty, or missing, all faces will be extracted

function [fv] = info2fv(BRinfo, which_faces)


if nargin<=1
	which_faces = 1:BRinfo.num_faces;
elseif isempty(which_faces)
	which_faces = 1:BRinfo.num_faces;
end

coord_ind = [1 2 3];

fv.vertices = make_vertices(coord_ind, BRinfo);

if ~isempty(BRinfo.sampler_data)
	fv.faces = make_faces_sampled(BRinfo, which_faces);
else
	fv.faces = make_faces_blocky(BRinfo, which_faces);
end

degen = any(diff(fv.faces(:,[1:3 1]),[],2)==0,2);
fv.faces(degen,:) = [];

end




function faces = make_faces_blocky(BRinfo, which_faces)

%
num_total_faces = 0;
for ii = 1:length(which_faces)
	
	curr_face = BRinfo.faces(which_faces(ii));
	num_total_faces = num_total_faces + curr_face.num_left + curr_face.num_right + curr_face.top>=0 + curr_face.bottom>=0;
end
num_total_faces = num_total_faces*2;
faces = zeros(num_total_faces,3);
curr_face_index = 1; %this will be incremented in the loop over the faces




for ii = 1:length(which_faces)
	
	curr_face = BRinfo.faces(which_faces(ii));
	
	if curr_face.midslice_index == -1
		continue
	end
	

	
	pass = 1; left_edge_counter = 1;  right_edge_counter = 1;
	
	
	
	while 1
		
		switch pass
			case 1  %the top edge
				pass = pass+1;
				
				
				if curr_face.top<0
					continue;
				end
				
				curr_edge = -10;
				if strcmp(curr_face.system_top,'input_critical_curve')
					curr_edge = BRinfo.crit_curve.edges(curr_face.top,:);
				elseif strcmp(curr_face.system_top,'input_surf_sphere')
					curr_edge = BRinfo.sphere_curve.edges(curr_face.top,:);
				else
					%do a lookup
					for zz = 1:length(BRinfo.singular_curves)
						if strcmp(BRinfo.singular_names{zz},curr_face.system_top)
							curr_edge = BRinfo.singular_curves{zz}.edges(curr_face.top,:);
						end
					end
					
				end
				
				if curr_edge<0
					continue;
				end
				
				curr_edge = curr_edge([3 2 1]);
				
			case 2  %the bottom edge
				pass = pass+1;
				
				if curr_face.bottom<0
					continue;
				end
				
				curr_edge = -10;
				if strcmp(curr_face.system_bottom,'input_critical_curve')
					curr_edge = BRinfo.crit_curve.edges(curr_face.bottom,:);
				elseif strcmp(curr_face.system_bottom,'input_surf_sphere')
					curr_edge = BRinfo.sphere_curve.edges(curr_face.bottom,:);
				else
					%do a lookup
					for zz = 1:length(BRinfo.singular_curves)
						if strcmp(BRinfo.singular_names{zz},curr_face.system_bottom)
							curr_edge = BRinfo.singular_curves{zz}.edges(curr_face.bottom,:);
						end
					end
					
				end
				
				
				if curr_edge<0
					continue;
				end
				
			case 3
				if left_edge_counter <= curr_face.num_left
					if curr_face.left(left_edge_counter)<0 %an error check
						continue;
					end
					
					slice_ind = curr_face.midslice_index; %offset by 1.
					edge_ind = curr_face.left(left_edge_counter); %offset by 1.
					
					curr_edge = BRinfo.critpoint_slices{slice_ind}.edges(edge_ind,:);
					left_edge_counter = left_edge_counter +1; %increment
					
					
				else
					pass = pass+1;
					continue;
				end
			case 4
				if right_edge_counter <= curr_face.num_right
					
					if curr_face.right(right_edge_counter)<0
						continue;
					end
					
					slice_ind = curr_face.midslice_index+1;
					edge_ind = curr_face.right(right_edge_counter);
					curr_edge = BRinfo.critpoint_slices{slice_ind}.edges(edge_ind,:);
					right_edge_counter = right_edge_counter +1;
					
					curr_edge = curr_edge([3 2 1]);
					
				else
					pass = pass+1;
					continue;
				end
			case 5
				break;
				
		end
		
		
		faces(curr_face_index,:) = [curr_edge(1) curr_edge(2) curr_face.midpoint];
		faces(curr_face_index+1,:) = [curr_edge(2) curr_edge(3) curr_face.midpoint];
		curr_face_index = curr_face_index+2;
		
		

	end

	
end


end















function [faces] = make_faces_sampled(BRinfo, which_faces)


faces = [];

if ~isempty(BRinfo.sampler_data)
	for ii = 1:length(which_faces)
		
		curr_face = which_faces(ii);
		
		f = BRinfo.sampler_data{curr_face};
		f(any(f<=0,2),:) = []; % omit problematic faces.
		faces = [faces;f];
	end
end


end




















function vertices = make_vertices(ind, BRinfo)



vertices = zeros(BRinfo.num_vertices,length(ind));


for ii=1:BRinfo.num_vertices
	vertices(ii,:) = real(transpose(BRinfo.vertices(ii).point(ind)));
end




end %re: function vertices












