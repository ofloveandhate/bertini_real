function h = plot_cond_hist(paths)

h = zeros(1,2);

h(1) = plot1(paths);
h(2) = plot2(paths);
h(3) = plot3(paths);
end

function h = plot1(paths)

	f = figure(1);
	ax = gca;

	cond = flatten(paths,'cond');
max(cond)
	h = histogram(log10(cond));

	set(ax,'FontSize',16);
	set(ax,'FontName','Courier')
	xlabel('$\log_{10}(\textnormal{cond})$','interpreter','latex')
	ylabel('count','interpreter','latex')

	title(sprintf('condition number each step -- %s',get_base_name()));
	figure(f)	
	
end


function h = plot2(paths)

	f = figure(2);
	ax = gca;

	cond = flatten(paths,'cond');
	time = flatten(paths,'time');
	
	h = scatter(real(time),log10(cond));

	set(ax,'FontSize',16);
	set(ax,'FontName','Courier')
	xlabel('real(time)','interpreter','latex')
	ylabel('$\log_{10}(\textnormal{cond})$','interpreter','latex')

	title(sprintf('condnum vs time -- %s',get_base_name()));
	figure(f)		
end


function h = plot3(paths)

	f = figure(3);
	ax = gca;

	cond = flatten(paths,'cond');
	path = flatten(paths,'path',3);
	time = flatten(paths,'time');

	q = max(real(path),[],2);
	whos
	h = scatter(-log10(abs(time)),log10(cond),20,max(real(path),[],2));

	set(ax,'FontSize',16);
	set(ax,'FontName','Courier')
	xlabel('-log10(real(time))','interpreter','latex')
	ylabel('$\log_{10}(\textnormal{cond})$','interpreter','latex')

	title(sprintf('condnum vs logtime -- %s',get_base_name()));
	colormap(jet)
	figure(f)		
end


function data = flatten(paths,fname,ndim)

if nargin <3
	ndim = 1;
end
data = [];
for ii = 1:length(paths)
	data = [data; paths{ii}.(fname)(:,1:ndim)];
end

end

function name = get_base_name()
[~, name, ~] = fileparts(pwd);
name = name(~isspace(name));
end