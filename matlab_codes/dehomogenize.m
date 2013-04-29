function result = dehomogenize(point)

result =[];

if length(point)<2
	display('point too short to dehomogenize');
	return;
end

result = point(2:end)/point(1);

end