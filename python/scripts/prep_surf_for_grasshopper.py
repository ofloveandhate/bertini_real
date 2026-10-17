

import bertini_real as br
import json
import os
import numpy as np


class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                            np.int16, np.int32, np.int64, np.uint8,
                            np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float16, np.float32,
                              np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

surf = br.data.read_most_recent()

surface_name = 'br_piece'  # os.getcwd().split('/')[-1]
surf.write_piece_data()

pieces = surf.separate_into_nonsingular_pieces()

for p in pieces:
    p.export_smooth() # basename=f"{surface_name}_smooth"
    p.export_raw() # basename=f"{surface_name}_raw"
    p.solidify_smooth(0.02, basename=f"{surface_name}_solidified_smooth_0.02") # basename=f"{surface_name}_solidified_smooth_0.02"
    p.solidify_raw(0.02, basename=f"{surface_name}_solidified_smooth_0.02") # basename=f"{surface_name}_solidified_raw_0.02"

    edge_pieces = p.edge_pieces()

    edge_pieces_as_points = []
    for ii,edge_piece in enumerate(edge_pieces):
        edge_pieces_as_points.append( (edge_piece.inputfilename()+"_"+str(ii),edge_piece.to_points()) ) 

    with open(p.generate_filename_no_ext(basename="touching_curves")+'.json','w') as f:

        # print(edge_pieces_as_points)
        json.dump(edge_pieces_as_points, f,cls=NumpyEncoder, indent=4) 