function result = dehomogenize(point)

if length(point)<2
	error('point too short to dehomogenize');
end

result = point(2:end)/point(1);

end
