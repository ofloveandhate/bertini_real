% plot some of the set of paths saved in 'paths'
function [handles,paths] = plot_paths(n)

% data_set = '_cauchy';
data_set = '_ps';
% data_set = '';
use_text = true;
noplot = 0;
method = 'complex'; % or 'complex'

method = 'mono';

% if length(n) > 1
% colors = jet(length(n));
% else 
% 	colors = [0 0 1];
% end

paths = cell(1,length(n));
handles = zeros(1,length(n));

hold off
for ii = 1:length(n)
	p = n(ii);
	
    display(p)
	[time, path, cond] = get_data(p,data_set);
	

	path = dehomogenize(path(:,1:end),2);
    
	if meh(path)
			continue
    end
      
    if strcmp(method,'real')
        [h,path_as_plotted] = path_colored_by_cond(path, time, cond, use_text,p);
    end
    if strcmp(method,'complex')
        [h,path_as_plotted] = path_x_complex(path, time, cond, use_text,p);
    end
    hold on
    
    
    if use_text
        plot_text(path_as_plotted, time)
    end
    paths{ii} = struct('time',time,'path',path,'cond',cond,'pathnum',p,'as_plotted',path_as_plotted, 'handle',h);
	
end
	
hold off
if ~noplot
% 	title('real part of path')
    
    endpoint = paths{end}.as_plotted(end,:);
%     
	a = 1;
	axis([-a+endpoint(1) a+endpoint(1) -a+endpoint(2) a+endpoint(2) -a+endpoint(3) a+endpoint(3)])

    view(2)
	cameratoolbar

end

view(2)
axis off
axis square
end

function h = plot_text(path, time)
    h = [];
    for ii = 0:100
        n = length(path(:,1)) - ii;
        t = time(n);
        if imag(t)==0
            txt = sprintf('      $t = %1.3d$',t);
        else
            txt = sprintf('      $t = %1.3d+%1.3d i$',real(t), imag(t));
        end
        
        h(end+1) = text(path(n,1),path(n,2),path(n,3),txt, interpreter='latex',fontsize=20, margin=10);
    end
end


function h = path_mono(path,color, use_text, path_num)

% 	h = patch(path(:,1),path(:,2),path(:,3),log10(cond)); % ,abs(data(:,8))
% 	set(h,'facecolor','none')
% 	set(h, 'edgecolor', 'interp');
% 	
	h = plot(path(:,1),path(:,2),color);
% 	set(h, 'linewidth', 5);
	if use_text
		t = text(path(end,1),path(end,2),path(end,3),sprintf('path %i',path_num));
% 		t.Color = colors(ii,:);
	end
end


function [h,path] = path_colored_by_cond(path, time, cond, use_text, path_num)


% 	h = patch(path(:,1),path(:,2),path(:,3),log10(cond)); % ,abs(data(:,8))
% 	set(h,'facecolor','none')
% 	set(h, 'edgecolor', 'interp');
% 	
    path = real(path);
	h = color_line(path(:,1),path(:,2),path(:,3),-log10(abs(time)));
    set(h, 'linewidth', 5);
    hold on
    h2 = plot3(path(:,1),path(:,2),path(:,3),'x');
    set(h2, 'MarkerSize', 20);
    
    
	if use_text
		t = text(path(end,1),path(end,2),path(end,3),sprintf('path %i',path_num));
% 		t.Color = colors(ii,:);
	end
end

function [h,path] = path_x_complex(path, time, cond, use_text, path_num)

% 	h = patch(path(:,1),path(:,2),path(:,3),log10(cond)); % ,abs(data(:,8))
% 	set(h,'facecolor','none')
% 	set(h, 'edgecolor', 'interp');
%      
    ind = 1;
    
    path = [real(path(:,ind)) imag(path(:,ind)) zeros(length(path(:,ind)),1)]; 
    
    
    
	h = color_line(path(:,1),path(:,2), path(:,3), -log10(abs(time)));
    set(h, 'linewidth', 5);
    hold on
    scatter3(path(:,1),path(:,2),path(:,3),100)
    
    view(2)
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
    path = path(:,[3,5,7])+1i*path(:,[4,6,8]);
% 	path = path(:,3:2:end-1)+1i*path(:,4:2:end-1);
	
	
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
	warning('empty path');
	c = true;
end

% if abs(imag(path(end-1,1)))> tol
%     warning('something else');
% 	c = true;
% end

end