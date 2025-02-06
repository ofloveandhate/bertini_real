function indicator = is_degenerate(obj)

if isstruct(obj) % then is face.  this is not programmed yet
	indicator = asdl;
else % is edge
	indicator = any([obj(:,1)==obj(:,2) obj(:,2)==obj(:,3) obj(:,1)==obj(:,3)],2);
end
	
end