function setupfig(br_plotter,varargin)
			
	if mod(length(varargin),2)~=0
		error('setupfig takes only pairs of arguments');
	end

	if ~isempty(which('getmondim'))
		mondim = getmondim;
		br_plotter.window = [20 20 mondim(3)-60 mondim(4)-130];
	end

	fig = figure();

	set(fig,'PaperPositionMode','auto');



	set(fig,'Name',sprintf('Bertini_real -- %s',br_plotter.filename),...
								'NumberTitle','off',...
								'Position',br_plotter.window);

	set(fig,'color','w');
	
	br_plotter.figures.main = fig;
end %re: setupfig
