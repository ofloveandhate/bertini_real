% fv2stl(src,hndl,fv)
%
% writes an fv to .stl file
%
% with interactive menu if you don't pass options as second argument
%
% example options structs:
%
% options.remove_duplicates = false;
% options.filename = 'asdf.stl'
%

function fv2stl(src,hndl,fv,filename)

if nargin==0
	error('must pass an fv to this function');
end

if nargin<3
	fv = src;
end


if nargin==2
	options = hndl;
end

if nargin==4
    options.filename = filename;
end

%delete the bad degenerate faces.  
degen = any(diff(fv.faces(:,[1:3 1]),[],2)==0,2);
fv.faces(degen,:) = [];

num_degen = sum(degen);
if num_degen>0
	display(sprintf('removed %d degenerate faces',num_degen));
end


options.filename = adjust_name(options.filename);

T = triangulation(fv.faces, fv.vertices);
stlwrite(T,sprintf('%s.stl',options.filename));

disp(sprintf('wrote to file %s',options.filename));
end


function name = adjust_name(old_name)

if strcmp(old_name(end-3:end),'.stl')
	name = old_name(1:end-4);
else
	name = old_name;
end

name = increment_name(name);
end

