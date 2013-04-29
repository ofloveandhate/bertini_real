function data = min_recursive(data)

level = ndims(data);

if level == 2
    if length(data)>1
    data = min(data);
    data = min_recursive(data);
    else
        data = min(data);
    end
else
    data = min(data,[],level);
    data = min_recursive(data);
end

end