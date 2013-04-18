function plot_paths(data,bound)

if nargin<2
	bound = 3;
	
	bound = calculate_bound(data);
end

fontsize  =16;
window = [40 40 1024 768];
% window = getmondim;


plot_params.fontsize = fontsize;
plot_params.window = window;
plot_params.format = 'jpg';
plot_params.format_flag = 'jpeg';



plot_imag = 0;


close all

fig1 = figure(1);
set(fig1,'Position',window,'PaperPositionMode','auto');




if data(1).num_variables>3
	num_plot_vars = 3;
	display('please code in a choice for the vars to plot');
	return;
else
	num_plot_vars = data(1).num_variables;
end

colors = jet(data(1).num_pts);

num_failures = 0;
for ii = 1:data(1).num_pts
	if data(ii).retval~=0
		num_failures = num_failures+1;
	end
end

	display(sprintf('this set has %d path failures\n',num_failures));
	
	

% first two are time
for ii = 1:data(1).num_pts
	if data(ii).retval==0
		endpt_marker = 'o';
		linestyle = '-';
	else
		endpt_marker = 'x';
		linestyle = '--';
	end
	
	if plot_imag
		subplot(1,2,1);
	end
	
	hold on
	
	if num_plot_vars==3
		g = plot3(data(ii).path(1,3),data(ii).path(1,5),data(ii).path(1,7),'Marker','p','MarkerSize',10);
		h = plot3(data(ii).path(:,3),data(ii).path(:,5),data(ii).path(:,7),'LineStyle',linestyle);
		k = plot3(data(ii).path(end,3),data(ii).path(end,5),data(ii).path(end,7),'Marker',endpt_marker,'MarkerSize',10);

	elseif num_plot_vars==2
		g = plot(data(ii).path(1,3),data(ii).path(1,5),'Marker','p','MarkerSize',10);	
		h = plot(data(ii).path(:,3),data(ii).path(:,5),'LineStyle',linestyle);
		k = plot(data(ii).path(end,3),data(ii).path(end,5),'Marker',endpt_marker,'MarkerSize',10);
	end
	set(g,'Color',colors(ii,:));
	set(h,'Color',colors(ii,:));
	
	
	set(gca,'FontSize',fontsize-2)
	
	if plot_imag
		subplot(1,2,2)
		hold on
		if num_plot_vars==3
			g = plot3(data(ii).path(1,4),data(ii).path(1,6),data(ii).path(1,8),'Marker','p','MarkerSize',10);
			h = plot3(data(ii).path(:,4),data(ii).path(:,6),data(ii).path(:,8),'LineStyle',linestyle);
			k = plot3(data(ii).path(end,4),data(ii).path(end,6),data(ii).path(end,8),'Marker',endpt_marker,'MarkerSize',10);
		elseif num_plot_vars==2
			g = plot(data(ii).path(1,4),data(ii).path(1,6),'Marker','p','MarkerSize',10);	
			h = plot(data(ii).path(:,4),data(ii).path(:,6),'LineStyle',linestyle);
			k = plot(data(ii).path(end,4),data(ii).path(end,6),'Marker',endpt_marker,'MarkerSize',10);
		end
		set(g,'Color',colors(ii,:));
		set(h,'Color',colors(ii,:));

	end
end

if plot_imag==1
	num_plots = 2;
else
	num_plots = 1;
end

for ii = 1:num_plots
	subplot(1,num_plots,ii);
	
	if num_plot_vars==3 
		view(3)
		plot3(0,0,0,'Marker','x','Color','k');
		xlabel(data(1).variable_names{2},'FontSize',fontsize);
		ylabel(data(1).variable_names{3},'FontSize',fontsize);
		zlabel(data(1).variable_names{4},'FontSize',fontsize);
		axis([-bound bound -bound bound -bound bound]);
	else
		view(2)
		plot(0,0,'Marker','x','Color','k','MarkerSize',10);
		xlabel(data(1).variable_names{2},'FontSize',fontsize);
		ylabel(data(1).variable_names{3},'FontSize',fontsize);
		axis([-bound bound -bound bound]);
	end
	set(gca,'FontSize',fontsize-2)

	if ii==1
		title(sprintf('Real; ODE %i EG %i MP %i',data(1).odepredictor, data(1).endgame, data(1).MPType),'FontSize',fontsize);
	else
		title(sprintf('Imaginary Part'),'FontSize',fontsize);
	end
	
end

render_into_file(sprintf('paths_linprod_ode%i_eg%i_mp%i',data(1).odepredictor,data(1).endgame,data(1).MPType),plot_params);

end






function bound = calculate_bound(data)


start_pts_real = zeros(data(1).num_pts,data(1).num_variables)
start_pts_imag = start_pts_real;

for ii = 1:data(1).num_pts
	start_pts_real(ii,:) = data(ii).path(1,3:2:end);
	start_pts_imag(ii,:) = data(ii).path(1,4:2:end);
end


start_pts_real


bound = 1.25*max_recursive(abs(start_pts_real))




end