function plot_separated_pieces(BRinfo)
    
pieces = separate_into_noncritical_pieces(BRinfo);


colors = jet(length(pieces));

for ii=1:length(pieces)
    
bertini_real_plotter('file',BRinfo,'whichfaces',pieces{ii},'mono',colors(ii,:),'touchingedgesonly',true);


fv = info2fv(BRinfo, pieces{ii}, true);

options.whichfaces=pieces{ii};
options.filename=sprintf('critpieces.stl');
fv2stl(fv,options);

end

options.whichfaces= 1:BRinfo.num_faces;
options.filename=sprintf('complete.stl');

fv = info2fv(BRinfo, options.whichfaces, true);
fv2stl(fv,options);
end