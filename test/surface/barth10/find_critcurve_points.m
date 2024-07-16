function [p,q] = find_critcurve_points(BRinfo)

p = {};
for ii=1:BRinfo.num_vertices
    if and(BRinfo.vertices(ii).input_filename_index == 0, BRinfo.vertices(ii).type==513)
        p{end+1} = BRinfo.vertices(ii);
    end
end

q = zeros(3,0);
for ii=1:length(p)
    q(:,end+1) = real(p{ii}.point(1:3));
end


scatter3(q(1,:),q(2,:),q(3,:))

end





