function gather_br_samples()



[BRinfo.dirname,BRinfo.dimension] = parse_dirname;







BRinfo = gather_vertices(BRinfo);



switch BRinfo.dimension
	
	case 1
		
		[BRinfo] = gather_curve(BRinfo.dirname, BRinfo);

	case 2
		
		[BRinfo] = gather_surface(BRinfo.dirname, BRinfo);
		
	otherwise
		
		
	
end

tmpnames = get_names(BRinfo.num_variables);
BRinfo.var_names = tmpnames(2:end);

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




function input = read_input(dirname,decomp)


input = fileread([dirname '/' decomp.inputfilename]);

end



function [BRinfo] = gather_surface(dirname, BRinfo)

BRinfo = parse_decomp(dirname,BRinfo);

% .num_faces, BRinfo.faces, BRinfo.midpoint_slices, BRinfo.critpoint_slices, BRinfo.crit_curve, BRinfo.sphere_curve,BRinfo.sampler_data

fin = fopen([dirname '/S.surf'],'r');

[BRinfo.num_faces] = fscanf(fin,'%i',[1 1]);
[BRinfo.num_edges] = fscanf(fin,'%i',[1 1]);
[num_midpoint_slices] = fscanf(fin,'%i',[1 1]);
[num_critpoint_slices] = fscanf(fin,'%i',[1 1]);
[num_singular_curves] = fscanf(fin,'%i',[1 1]);
BRinfo.singular_curve_multiplicities = [];
for ii = 1:num_singular_curves
    BRinfo.singular_curve_multiplicities(ii) = fscanf(fin,'%i',[1 1]);
end

% BRinfo.midpoint_slices = [];
% BRinfo.critpoint_slices = [];
% BRinfo.crit_curve = [];

fclose(fin);



[BRinfo.faces] = gather_faces(dirname);


for ii =1:num_midpoint_slices
    a = gather_curve([dirname '/curve_midslice_' num2str(ii-1) ], []);
	BRinfo.midpoint_slices(ii) = a;
end

for ii =1:num_critpoint_slices
    a = gather_curve([dirname '/curve_critslice_' num2str(ii-1) ],[]);
	[BRinfo.critpoint_slices(ii)] = a;
end

a = gather_curve([dirname '/curve_crit'],[]);
[BRinfo.crit_curve] = a;


[BRinfo.sphere_curve] = gather_curve([dirname '/curve_sphere'],[]);

% BRinfo.singular_curves = [];
for ii = 1:num_singular_curves
    a = gather_curve([dirname '/curve_singular_mult_' num2str(BRinfo.singular_curve_multiplicities(ii))],[]);
    [BRinfo.singular_curves(ii)] = a;
    BRinfo.singular_names{ii} = a.inputfilename;
end


BRinfo.sampler_data = gather_surface_samples(dirname);


BRinfo.input = read_input(dirname,BRinfo);
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
	
	faces(ii).system_top = fscanf(fid,'%s',[1 1]);
	faces(ii).system_bottom = fscanf(fid,'%s\n',[1 1]);
	
	
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
	

end

fclose(fid);

end


function BRinfo = parse_decomp(dirname,BRinfo)




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


function [dirname,dimension] = parse_dirname()


if isempty('Dir_Name')
	display('no file ''Dir_Name''.  please run bertini_real');
	return;
end

fid = fopen('Dir_Name','r');
dirname = fscanf(fid,'%s',[1 1]);
MPtype = fscanf(fid,'%i',[1 1]);
dimension = fscanf(fid,'%i',[1 1]);

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

function [curve] = gather_curve(dirname, curve)

curve = parse_decomp(dirname,curve);


fid = fopen([dirname '/E.edge'],'r');

curve.num_edges = fscanf(fid,'%i\n',[1 1]);


curve.edges = zeros(curve.num_edges,3);  % 3, as left, mid, right
for ii = 1:curve.num_edges
	tmp = fscanf(fid,'%i',[1 3]);
	curve.edges(ii,:) = tmp+1;
end
	
fclose(fid);



filename = [dirname '/' 'samp.curvesamp'];


if isempty(dir(filename))
	curve.sampler_data = [];
else
	fid = fopen(filename,'r');
	curve.num_edges = fscanf(fid,'%i\n\n',[1 1]);
	curve.sampler_data.sample_sizes = zeros(curve.num_edges,1);
    curve.sampler_data.edge = [];
	for ii = 1:curve.num_edges
		num_samples = fscanf(fid,'%i\n',[1 1]);
		curve.sampler_data.sample_sizes(ii) = num_samples;
		for jj=1:num_samples
			tmp_ind = fscanf(fid,'%i\n',[1 1]);
			curve.sampler_data.edge(ii).samples(jj) = tmp_ind;
		end
	end
	fclose(fid);
    
    display('done gathering sampler data.');
end


curve.input = read_input(dirname,curve);
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




