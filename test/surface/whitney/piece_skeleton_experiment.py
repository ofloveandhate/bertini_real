import bertini_real as br

import json
import numpy as np

br.data.gather_and_save()
surface = br.data.read_most_recent()
surface_pieces = surface.separate_into_nonsingular_pieces()




class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                            np.int16, np.int32, np.int64, np.uint8,
                            np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32,
                              np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


import os
surface_name = os.getcwd().split('/')[-1]

for surface_piece in surface_pieces:
    surface_piece.export_smooth(basename="_thin_smooth")
    surface_piece.solidify_smooth(0.08, basename="_solidified_smooth_0.08")
    edge_pieces = surface_piece.edge_pieces()

    edge_pieces_as_points = []
    for ii,edge_piece in enumerate(edge_pieces):
        edge_pieces_as_points.append( (edge_piece.inputfilename()+"_"+str(ii),edge_piece.to_points()) ) 

    with open(surface_piece.generate_filename_no_ext(basename="touching_curves")+'.json','w') as f:

        # print(edge_pieces_as_points)
        json.dump(edge_pieces_as_points, f,cls=NumpyEncoder, indent=4)  