% make_nsolns_colorbar(maxiii,fontsize)
%
% makes a colorbar for a data set, and sets a colormap to correspond.
% the data set must be integer, and be positive (0 ok).
%
% input: integer maxiii, the upper limit for the colorbar.  must be an
%            integer.
%        integer fontsize, an integer for the size of the labels, etc.
%
% output: none
%
% default fontsize is 16. 
%
% daniel brake
% colorado state university
% mathematics
% 2013
% danielthebrake@gmail.com


function cbar = make_discrete_colorbar(miniii,maxiii,fontsize,labeltxt)

max_num_bins = 7;

if nargin<3
	fontsize = 16;
	labeltxt = '';
elseif nargin<4
	labeltxt = '';
end





% if (or(and(miniii>0,maxiii>0),and(miniii<0,maxiii<0)))
% 	range = maxiii-miniii;
% else
	range = maxiii-miniii+1; %include a color for 0
% end


% if maxiii ==0
% 	maxiii = 1;
% end




J = jet(range);

colormap(J);
set(gca,'CLim',[0 1]);

cbar = colorbar;
f = factor(range);



if range<=max_num_bins
	labels = miniii:maxiii;
elseif min(f)>max_num_bins
	num_bins = max_num_bins;
	labels = miniii:range/num_bins:maxiii;
	labels = ceil(labels);
	if labels(end)>((max_num_bins-2)/max_num_bins*maxiii)
		labels = [labels(1:end-1) maxiii];
	else
		labels = [labels(1:end) maxiii];
		
	end

else
	
	small_factors = f(f<max_num_bins);
	
	if prod(small_factors)>max_num_bins
		if max(small_factors)>4
			num_bins = max(small_factors);
		else
			num_bins = 1;
			counter = 1;
			while (num_bins*f(counter)<max_num_bins)
				num_bins= num_bins*f(counter);
				counter = counter+1;
			end
		end
	else
		num_bins = prod(small_factors);
	end
	

	labels = miniii:range/num_bins:maxiii; 
end

placeme = maxiii*(labels-min(labels))/(range-1);
placements = placeme/(maxiii)*((2*(maxiii+1)-1)/(2*(maxiii+1)) - 1/(2*(maxiii+1))) + 1/(2*(maxiii+1)); % +
%rescale labels

if range==1
	placements = 0.5;
end

set(cbar,'Ytick',placements,'YTicklabel',labels);
set(gca,'FontSize',fontsize-2);

if ~isempty(labeltxt)
	set(get(cbar,'ylabel'),'String', labeltxt,'FontSize',fontsize-2);
end


end