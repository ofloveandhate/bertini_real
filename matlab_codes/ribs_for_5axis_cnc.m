% [] = ribs_for_5axis_cnc(BRinfo)
%
% this function writes to a file the sequence of points and normal vectors used
% for a 5-axis cnc machine.  it was written for one purpose: cutting one lobe of
% the barth sextic on Edmund Harris' machine.  maybe it's useful for other
% things, but probably not.
%
% 
%
% Danielle Amethyst, 2019

function [] = ribs_for_5axis_cnc(BRinfo,varargin)

sanity_checks(BRinfo);

options = parse_options(BRinfo,varargin{:});

J = generate_evaluable_jacobian(BRinfo);

ribs = cell(1,length(options.which_faces));
for ii = 1:length(options.which_faces)
	face_index = options.which_faces(ii);
	ribs{ii} = generate_zigzag(BRinfo,face_index,J);
end

visualize_ribs(ribs,options);
write_ribs(ribs,options);

end


function output = generate_zigzag(BRinfo, face_index, J)

ribs = BRinfo.ribs{face_index}; % unpack

num_ribs = length(ribs);
rib_sizes = zeros(1,num_ribs);
for ii = 1:num_ribs
	rib_sizes = length(ribs{ii});
end

output = zeros(6,sum(rib_sizes)); % preallocate

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
		
		output(:,c) = [x;n];
		c = c+1;
	end
end

end 




function n = compute_normal_vector(x,jac)

nullspace = null(jac(x));
n = cross(nullspace(:,1),nullspace(:,2));

end


function [] = visualize_ribs(ribs,options)
hold off
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
end

daspect([1 1 1 ]);
end

function write_ribs(ribs,options)

base_name = options.base_filename;

for ii = 1:length(options.which_faces)
	face_index = options.which_faces(ii);
	fid = fopen(sprintf('%s_%i',base_name,ii),'w');
	
	fclose(fid);
end

end


function options = parse_options(BRinfo,varargin)

options = default_options(BRinfo);
options.base_filename = 'toolpaths';

if mod(length(varargin),2)~=0
	error('options must come in pairs -- name, argument.  you have a mismatch');
end

opt_ind = 1;
while opt_ind < length(varargin)
	arg_ind = opt_ind+1;
	switch varargin{opt_ind}
		case {'which_faces','whichfaces'}
			options.which_faces = varargin{arg_ind};
	end
	opt_ind = opt_ind + 2;
end


end


function options = default_options(BRinfo)

options.which_faces = 1:BRinfo.num_faces;


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

