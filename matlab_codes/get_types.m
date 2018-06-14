% unpacks all types of a vertex stored in a BRinfo

function t = get_types(index, BRinfo)

v = BRinfo.vertices(index);

has_type = @(x,y) bitand(x.type,y)>0;

nums = BRinfo.vertex_types.nums;
names = BRinfo.vertex_types.names;

counter = 0;
for ii = 1:length(nums)
	if has_type(v, nums(ii))
		counter = counter+1;
		t{counter} = names{ii};
	end
end

end