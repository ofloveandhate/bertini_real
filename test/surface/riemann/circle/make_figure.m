% generates plots for the NAG undergrad book
% silviana amethyst, 2019

function [handles] = make_figure()

save_figs = false;

handles(1) = main_plot('stl/untitled-Xr_Xi_Yr.stl',1,save_figs);
handles(2) = main_plot('stl/untitled-Xr_Xi_Yi.stl',2,save_figs);

end

function h = main_plot(stl_name,fig_ind,save_figs)


f = figure(fig_ind);

tri = stlread(stl_name);

h = trimesh(tri);


set(h,'FaceColor',[1 0.6 0.6]);
set(h,'FaceAlpha',0.9)
set(h,'EdgeColor','none')
set(h,'FaceLighting','gouraud')


camlight
cameratoolbar

% light('Position',[1 0 0])
adjust_figure(f);

t = linspace(0,2*pi,100);
x = cos(t);
y = zeros(size(t));
z = sin(t);

if fig_ind==1
	hold on
	circ = pipe_surface(x,y,z);
	set(circ,'FaceColor','w');
	set(circ,'FaceAlpha',1)
	set(circ,'EdgeColor','none')
	set(circ,'FaceLighting','gouraud')
	zlabel('real(y)')
else
	zlabel('imag(y)')
end


xlabel('real(x)')
ylabel('imag(x)')


view(116.3, 13);
view(43, 17);
view(21, 8);

if save_figs
	w = render_into_file('gendef');
	w.format = 'png';
	w.format_flag = 'png';
	w.resolution = 300;
	w.basename = 'circle_as_complex_curve';

	render_into_file(w)
end

hold off

end

function adjust_figure(f)
f.Units = 'inches';
% f
f.Position = [5 5 4 3];

daspect([1 1 1])
end