% plot some of the set of paths saved in 'paths'
function [handles,paths] = plot_path(n)

% data_set = '_cauchy';
% data_set = '_ps';
data_set = '';
use_text = 1;
noplot = 0;

% if length(n) > 1
% colors = jet(length(n));
% else 
% 	colors = [0 0 1];
% end

paths = cell(1,length(n));
handles = zeros(1,length(n));
for ii = 1:length(n)
	p = n(ii);
	
	[time, path, cond] = get_data(p,data_set);
	
	paths{ii} = struct('time',time,'path',path,'cond',cond,'pathnum',p);


	path = real(dehomogenize(path(:,1:end),2));

	if meh(path)
			continue
	end
	

	if ~noplot
		handles(ii) = path_colored_by_cond(path,cond, use_text,p);
		hold on
		
		

	end
	
end
	

	
if ~noplot
	title('real part of path')

	a = 10000;
	axis([-a a -a a -a a])
	cameratoolbar
end

end


function h = path_colored_by_cond(path,cond, use_text, path_num)

	h = patch(path(:,1),path(:,2),path(:,3),log10(cond)); % ,abs(data(:,8))
	set(h,'facecolor','none')
	set(h, 'edgecolor', 'interp');
	set(h, 'linewidth', 5);
	if use_text
		t = text(path(end,1),path(end,2),path(end,3),sprintf('path %i',path_num));
% 		t.Color = colors(ii,:);
	end
end

function [time, path, cond] = get_data(n,data_set)

path = [];

	fid = fopen(sprintf('paths%s/path_%i',data_set,n),'r');
	while 1
		ell = fgetl(fid);
		temp = str2num(ell);
		if feof(fid)
			break
		end

		path = [path ; temp];
	end
% 	path = [path ;nan(size(path,2))];
	
	fclose(fid);

	cond = path(:,end);
	time = path(:,1)+1i*path(:,2);
	path = path(:,3:2:end-1)+1i*path(:,4:2:end-1);
	
	
end


function rgba = cdata2rgb(cdata)

rgba = colormap('parula'); % take your pick (doc colormap)
rgba = interp1(linspace(min(cdata),max(cdata),length(rgba)),rgba,cdata); % map color to y values
rgba = uint8(rgba'*255); % need a 4xN uint8 array
rgba(4,:) = 255; % last column is transparency

end

function c = meh(path)

tol = 1e-11;
c = false;

if isempty(path)
	warning('empty path %i',p);
	c = true;
end

if abs(imag(path(end-1,1)))> tol
	c = true;
end

end