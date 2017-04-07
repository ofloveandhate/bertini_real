% produces a cell array of fv structs, one for each 
% critically separated section of a decomposed surface


function fvs = extract_noncritical_fv(BRsurf)

inds = separate_into_noncritical_pieces(BRsurf);
num_blobs = length(inds);

fvs = cell(num_blobs,1);

for ii = 1:num_blobs
	fvs{ii} = info2fv(BRsurf, inds{ii});
end

end