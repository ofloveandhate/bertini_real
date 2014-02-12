function gather_br_samples()



dirname = parse_dirname;

BRinfo = parse_decomp(dirname);

BRinfo.dirname = dirname;




BRinfo = gather_vertices(BRinfo);

tmpnames = get_names(BRinfo.num_variables);
BRinfo.var_names = tmpnames(2:end);

switch BRinfo.dimension
	
	case 1
		
		[BRinfo.num_edges, BRinfo.edges, BRinfo.sampler_data] = gather_curve(BRinfo.dirname);

	case 2
		
		[BRinfo.num_faces, BRinfo.faces, BRinfo.midpoint_slices, BRinfo.critpoint_slices, BRinfo.crit_curve, BRinfo.sphere_curve,BRinfo.sampler_data] = gather_surface(BRinfo.dirname);
		
	otherwise
		
		
	
end

display('done gathering decomposition.');






prev_filenames = dir('BRinfo*.mat');
max_found = -1;

for ii = 1:length(prev_filenames)
	curr_name = prev_filenames(ii).name;
	curr_num = str2num(curr_name(7:end-4));
	if max_found < curr_num
		max_found = curr_num;
	end
	
end

filename = ['BRinfo' num2str(max_found+1) '.mat'];

display('writing data to file.');
save(filename,'BRinfo');

display(['file saved to ' filename]);
end%re: function








function [num_faces, faces, midpoint_slices, critpoint_slices, crit_curve, sphere_curve,sampler_data] = gather_surface(dirname)



fin = fopen([dirname '/S.surf']);

[num_faces] = fscanf(fin,'%i',[1 1]);
[num_edges] = fscanf(fin,'%i',[1 1]);
[num_midpoint_slices] = fscanf(fin,'%i',[1 1]);
[num_critpoint_slices] = fscanf(fin,'%i',[1 1]);
midpoint_slices = [];
critpoint_slices = [];
crit_curve = [];

fclose(fin);



[faces] = gather_faces(dirname);


for ii =1:num_midpoint_slices
	[midpoint_slices(ii).num_edges, midpoint_slices(ii).edges, midpoint_slices(ii).sampler_data] = gather_curve([dirname '/curve_midslice_' num2str(ii-1) ]);
end

for ii =1:num_critpoint_slices
	[critpoint_slices(ii).num_edges, critpoint_slices(ii).edges, critpoint_slices(ii).sampler_data] = gather_curve([dirname '/curve_critslice_' num2str(ii-1) ]);
end

[crit_curve.num_edges, crit_curve.edges,crit_curve.sampler_data] = gather_curve([dirname '/curve_crit']);

[sphere_curve.num_edges, sphere_curve.edges, sphere_curve.sampler_data] = gather_curve([dirname '/curve_sphere']);

sampler_data = gather_surface_samples(dirname);
end



function [faces] = gather_faces(dirname)

faces = [];
fid = fopen([dirname '/F.faces'],'r');

num_faces = fscanf(fid,'%i\n',[1 1]);

for ii = 1:num_faces
	[faces(ii).midpoint] = fscanf(fid,'%i', [1 1]);
	faces(ii).midslice_index = fscanf(fid,'%i\n', [1 1]);
	
	faces(ii).top = fscanf(fid,'%i',[1 1]);
	faces(ii).bottom = fscanf(fid,'%i\n',[1 1]);
	
	faces(ii).system_top = fscanf(fid,'%i',[1 1]);
	faces(ii).system_bottom = fscanf(fid,'%i\n',[1 1]);
	
	
	faces(ii).num_left = fscanf(fid,'%i\n',[1 1]); 
	faces(ii).left = zeros(1,faces(ii).num_left); % preallocate
	for jj = 1:faces(ii).num_left
		faces(ii).left(jj) = fscanf(fid,'%i\n',[1 1]);
	end
	
	[faces(ii).num_right] = fscanf(fid,'%i\n',[1 1]);
	faces(ii).right = zeros(1,faces(ii).num_right);
	for jj = 1:faces(ii).num_right
		faces(ii).right(jj) = fscanf(fid,'%i\n',[1 1]);
	end
	
% 	faces(ii).num_bottom
% 	faces(ii).num_top

end

fclose(fid);

end


function BRinfo = parse_decomp(dirname)




fid = fopen([dirname '/' 'decomp'],'r');

% nchars = fscanf(fid,'%i\n',[1 1]);
BRinfo.inputfilename = fscanf(fid,'%s\n',[1 1]);


BRinfo.num_variables = fscanf(fid,'%i',[1 1]);
BRinfo.dimension = fscanf(fid,'%i',[1 1]);

BRinfo.vertex.num_types = fscanf(fid,'%i',[1 1]);


BRinfo.vertex.types = zeros(BRinfo.vertex.num_types,1);
BRinfo.vertex.counters = zeros(BRinfo.vertex.num_types,1);
BRinfo.vertex.indices = zeros(BRinfo.vertex.num_types,1);

for ii =1:BRinfo.vertex.num_types
	BRinfo.vertex.types(ii) = fscanf(fid,'%i',[1 1]);
	BRinfo.vertex.counters(ii) = fscanf(fid,'%i\n',[1 1]);
	for jj = 1:BRinfo.vertex.counters(ii)
		BRinfo.vertex.indices(ii,jj) = fscanf(fid,'%i\n',[1 1]);
	end
	
end


BRinfo.pi = zeros(BRinfo.num_variables,BRinfo.dimension);

for jj = 1:BRinfo.dimension
    tmp1 = fscanf(fid,'%i',[1 1]);
    tmp_num_vars = tmp1;
for ii = 1:tmp_num_vars
	tmp1 = fscanf(fid,'%e',[1 1]);
	tmp2 = fscanf(fid,'%e\n',[1 1]); % this should scan right over the empty line separating them
	BRinfo.pi(ii,jj) = tmp1+1i*tmp2;
end
end


BRinfo.pi=BRinfo.pi(2:end,:);


BRinfo.num_patches = fscanf(fid,'%i',[1 1]);
for jj = 1:BRinfo.num_patches
    BRinfo.patch.sizes(jj) = fscanf(fid,'%i',[1 1]);
    for ii = 1:BRinfo.patch.sizes(jj)
        tmp1 = fscanf(fid,'%e',[1 1]);
        tmp2 = fscanf(fid,'%e\n',[1 1]); % this should scan right over the empty line separating them
        BRinfo.patch.vectors{jj}(ii) = tmp1 + 1i*tmp2;
    end
end

tmp1 = fscanf(fid,'%e',[1 1]);
tmp2 = fscanf(fid,'%e\n',[1 1]); % this should scan right over the empty line separating them
BRinfo.radius = tmp1+1i*tmp2;
    
center_size = fscanf(fid,'%i',[1 1]);
BRinfo.center = zeros(1,center_size);
for ii = 1:center_size
    tmp1 = fscanf(fid,'%e',[1 1]);
    tmp2 = fscanf(fid,'%e\n',[1 1]); % this should scan right over the empty line separating them
    BRinfo.center(ii) = tmp1 + 1i*tmp2;
end
fclose(fid);





end


function dirname = parse_dirname()


if isempty('Dir_Name')
	display('no file ''Dir_Name''.  please run bertini_real');
	return;
end

fid = fopen('Dir_Name','r');
nchars = fscanf(fid,'%i\n',[1 1]);
dirname = fscanf(fid,'%s\n',[1 1]);
nchars = fscanf(fid,'%i\n',[1 1]);

fclose(fid);


end


function [sampler_data] = gather_surface_samples(dirname)
sampler_data = [];

if isempty(dir([dirname '/samp.surfsamp']))
    return;
end

fid = fopen([dirname '/samp.surfsamp']);

num_faces = fscanf(fid,'%i',[1 1]);
sampler_data = cell(num_faces,1);

for ii = 1:num_faces
    
    
    num_triangles_this_edge = fscanf(fid,'%i',[1 1]);
    
    sampler_data{ii} = zeros(num_triangles_this_edge,3);
    for jj = 1:num_triangles_this_edge
        temp_triangle = fscanf(fid,'%i',[1 3]);
        sampler_data{ii}(jj,:) = temp_triangle;
    end
end


fclose(fid);



end

function [num_edges, edges, sampler_data] = gather_curve(dirname)
fid = fopen([dirname '/E.edge'],'r');

num_edges = fscanf(fid,'%i\n',[1 1]);


edges = zeros(num_edges,3);  % 3, as left, mid, right
for ii = 1:num_edges
	tmp = fscanf(fid,'%i',[1 3]);
	edges(ii,:) = tmp+1;
end
	
fclose(fid);



filename = [dirname '/' 'samp.curvesamp'];


if isempty(dir(filename))
	display('no sampler data.  run sampler.');
	sampler_data = [];
else
	fid = fopen(filename,'r');
	num_edges = fscanf(fid,'%i\n\n',[1 1]);
	sample_sizes = zeros(num_edges,1);
	for ii = 1:num_edges
		num_samples = fscanf(fid,'%i\n',[1 1]);
		sample_sizes(ii) = num_samples;
		for jj=1:num_samples
			tmp_ind = fscanf(fid,'%i\n',[1 1]);
			edge(ii).samples(jj) = tmp_ind;
		end
	end
	fclose(fid);
	sampler_data.edge = edge;
	sampler_data.sample_sizes = sample_sizes;
    
    display('done gathering sampler data.');
end

end


function [BRinfo] = gather_vertices(BRinfo)

fname = 'V.vertex';
if ~isempty(dir(sprintf('%s/V_samp.vertex',BRinfo.dirname)))
    fname = 'V_samp.vertex';
end

fid = fopen(sprintf('%s/%s',BRinfo.dirname,fname));
BRinfo.num_vertices = fscanf(fid,'%i',[1 1]);
num_projections = fscanf(fid,'%i',[1 1]);
num_natural_vars = fscanf(fid,'%i',[1 1]);
num_filenames = fscanf(fid,'%i',[1 1]);
BRinfo.vertices = repmat(struct('point',[]),[1 BRinfo.num_vertices]);



for ii = 1:num_projections
	for jj = 1:num_natural_vars
		tmp = fscanf(fid,'%e %e\n',[1 2]);
	end
end

BRinfo.input_filenames = cell(1,num_filenames);
for ii = 1:num_filenames
    nchars = fscanf(fid,'%i\n',[1 1]);
   BRinfo.input_filenames{ii} = fscanf(fid,'%s\n',[1 1]); 
end

for ii = 1:BRinfo.num_vertices
	tmp_num_variables = fscanf(fid,'%i',[1 1]); % number variables
	tmpvertex = zeros(tmp_num_variables-1,1);
	for jj = 1:tmp_num_variables
		tmp = fscanf(fid,'%e %e\n',[1 2]);
		tmpvertex(jj) = tmp(1)+1i*tmp(2);
	end
	
	
	BRinfo.vertices(ii).point = dehomogenize(tmpvertex);
	
	num_proj_vals = fscanf(fid,'%i\n',[1 1]);
	for jj = 1:num_proj_vals
		tmp = fscanf(fid,'%e %e\n',[1 2]);
		BRinfo.vertices(ii).projection_value(jj) = tmp(1)+1i*tmp(2);
	end
    BRinfo.vertices(ii).input_filename_index = fscanf(fid,'%i',[1 1]);
	BRinfo.vertices(ii).type = fscanf(fid,'%i',[1 1]);
end
fclose(fid);


end




