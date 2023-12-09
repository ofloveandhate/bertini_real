% [] = ribs_for_5axis_cnc(BRinfo)
%
% this function writes to a file the sequence of points and normal vectors used
% for a 5-axis cnc machine.  it was written for one purpose: cutting one lobe of
% the barth sextic on Edmund Harris' machine.  maybe it's useful for other
% things, but probably not.
%
% options you can set:
% * `which_faces`, `whichfaces` -- an array of indices of faces to use
% * `render` -- a boolean, indicates whether we'll render the paths to the
%			screen.  default is true.
%
% * 'normals_mode', 'normals' -- modes: 'outward_from_origin', 'none';
%
% this list of options may be incomplete because i'm lazy.  see `parse_options`
% in this .m file.s
%
% silviana Amethyst, 2019

function [] = ribs_for_5axis_cnc(BRinfo,varargin)

sanity_checks(BRinfo);

options = parse_options(BRinfo,varargin{:});

J = generate_evaluable_jacobian(BRinfo);

ribs = cell(1,length(options.which_faces));
for ii = 1:length(options.which_faces)
	face_index = options.which_faces(ii);
	ribs{ii} = generate_zigzag(BRinfo,face_index,J,options);
end


if options.render
	visualize_ribs(ribs,options);
end

write_ribs(ribs,options);

end


function toolpath = generate_zigzag(BRinfo, face_index, J, options)



ribs = BRinfo.ribs{face_index}; % unpack

num_ribs = length(ribs);
rib_sizes = zeros(1,num_ribs);
for ii = 1:num_ribs
	rib_sizes = length(ribs{ii});
end

alignment_function = make_alignment_handle(options); % ooh, meta
affine_transform = make_affine_transform_handle(options);

toolpath = zeros(6,sum(rib_sizes)); % preallocate

c = 1;
parity = 1;
for ii = 1:num_ribs
	r = ribs{ii};
	if parity == 1
		parity = -1;
	else
		r = r(end:-1:1);
		parity = 1;
	end
	for jj = 1:length(r)
		x = real(BRinfo.vertices(r(jj)).point(1:3));
		n = compute_normal_vector(x,J);

		[x,n] = affine_transform(x,n);
		n = alignment_function(n,x);
		
		toolpath(:,c) = [x;n];
		c = c+1;
	end
end

end 


function n = align_normal_outward(n,x)

if dot(n,x) < 0
	n = -n;
end

end

function n = compute_normal_vector(x,jac)

nullspace = null(jac(x));
n = cross(nullspace(:,1),nullspace(:,2));
end


function [] = visualize_ribs(ribs,options)
hold_cached = ishold;

for ii = 1:length(ribs)
	r = ribs{ii};
	x = r(1,:);
	y = r(2,:);
	z = r(3,:);
	u = r(4,:);
	v = r(5,:);
	w = r(6,:);
	
	quiver3(x,y,z,u,v,w);
	hold on
	plot3(x,y,z);
	hold on
end

daspect([1 1 1 ]);

if hold_cached
	hold on
else
	hold off
end


end

function write_ribs(ribs,options)

base_name = options.base_filename;

fprintf('writing files to files named `%s_*`\n',base_name);

for ii = 1:length(options.which_faces)
	fid = fopen(sprintf('%s_%i',base_name,ii),'w');
	write_rib(fid,ribs{ii});
	fclose(fid);
end

end


function write_rib(fid,rib)
[~,n] = size(rib);
fprintf(fid,'%i\n',n);
for ii = 1:n
	fprintf(fid,'%1.16f ',rib(:,ii));
	fprintf(fid,'\n');
end
end

function [affine_transf] = make_affine_transform_handle(options)

if options.rotate.flag
	R = make_rotation_matrix(options.rotate.axis,options.rotate.angle);
	r = @(x,n) rotate_xn(x,n,R);
else
	r =  @parrot_2;
end


if options.translate.flag
	t = @(x,n) translate_xn(x,n,options.translate.vector);
else
	t = @parrot_2;
end

if options.rotate_first
	first = r;
	second = t;
else
	first = t;
	second = r;
end

affine_transf = @(x,n) double_func_comp(x,n,first,second);

end

function [x,n] = double_func_comp(x,n,first,second)
[x,n] = first(x,n);
[x,n] = second(x,n);
end

function [x,n] = translate_xn(x,n,t)
x = x+t;
end

% see also
% https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
function [x,n] = rotate_xn(x,n,R)
x = R*x;
n = R*n;
end


function R = make_rotation_matrix(axis,angle)

x = axis(1);
y = axis(2);
z = axis(3);

c = cos(angle);
s = sin(angle);
R = [c+x^2*(1-c),    x*y*(1-c)-z*s,  x*z*(1-c)+y*s;...
	 y*x*(1-c)+z*s,  c+y^2*(1-c),    y*z*(1-c)-x*s;...
	 z*x*(1-c)-y*s,  z*y*(1-c)+x*s,  c+z^2*(1-c)];
end


function alignment_function = make_alignment_handle(options)

switch options.normals_mode
	case 'outward_from_origin'
		alignment_function = @align_normal_outward;
	case 'none'
		alignment_function = @(n,x) n;
	otherwise
		error('invalid choice of normals_mode');
end


end




% the un-function.  aaaraaarp!
function [a,b] = parrot_2(a,b,varargin)

end


function options = parse_options(BRinfo,varargin)

options = default_options(BRinfo);
options.base_filename = 'toolpaths';
options.translate.flag = false;
options.rotate.flag = false;
options.rotate_first = true;
options.affine_precedes_ribbing = true;

if mod(length(varargin),2)~=0
	error('options must come in pairs -- name, argument.  you have a mismatch');
end

opt_ind = 1;
while opt_ind < length(varargin)
	arg_ind = opt_ind+1;
	switch varargin{opt_ind}
		case {'which_faces','whichfaces'}
			options.which_faces = varargin{arg_ind};
		case {'render'}
			options.render = varargin{arg_ind};
		case {'normals','normals_mode'}
			options.normals_mode = varargin{arg_ind};
		case {'translate','translate_by'}
			options.translate.flag = true;
			options.translate.vector = reshape(varargin{arg_ind},[3 1]);
		case {'axis','rotation_axis'}
			options.rotate.flag = true;
			options.rotate.axis = varargin{arg_ind}/norm(varargin{arg_ind});
		case {'angle','theta','rotation_angle'}
			options.rotate.flag = true;
			options.rotate.angle = varargin{arg_ind};
		case {'rotate_translate_first'}
			options.rotate_first = varargin{arg_ind};
	end
	opt_ind = opt_ind + 2;
end


end


function options = default_options(BRinfo)

options.which_faces = 1:BRinfo.num_faces;
options.render = true; 
options.normals_mode = 'outward_from_origin';
end



function sanity_checks(BRinfo)
if BRinfo.dimension~=2
	error('sorry, this only works for surfaces');
end
if BRinfo.num_variables~=4 % there is an off-by-one because of the homogenizing variable.  deal with it.
	error('too many variables, you have %i natural vars, need 3',BRinfo.num_variables);
end

if ~isfield(BRinfo,'ribs')
	error('your data doesn''t have the rib data.  bummer.  run sampler with the `-saveribs` option and re-gather');
end


end

