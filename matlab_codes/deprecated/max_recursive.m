function data = max_recursive(data)

level = ndims(data);

if level == 2
    if length(data)>1
    data = max(data);
    data = max_recursive(data);
    else
        data = max(data);
    end
else
    data = max(data,[],level);
    data = max_recursive(data);
end

end