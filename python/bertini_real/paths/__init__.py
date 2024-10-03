
"""
    :platform: Unix, Windows, MacOS
    :synopsis: Utilities for plotting paths computed by bertini_real
"""


import numpy as np

class Path(object):
    """docstring for Path"""
    def __init__(self, filename):

        super(Path, self).__init__()
        self.path = None
        self.num = -1

        self.time = None
        self.space = None
        self.condition_number = None



        self.num = int(filename.split('_')[1])
        self._read(filename)


    def _read(self, filename):

        with open(filename,'r') as file:
            raw = file.read()

        as_lines = raw.split('\n')

        path_as_numbers = []
        for line in as_lines:
            if line:
                path_as_numbers.append( np.array( [float(n) for n in line.strip().split()] ) )

        all_together = np.array(path_as_numbers)

        self.time = all_together[:,0] + 1j*all_together[:,1]

        self.condition_number = all_together[:,-1]
        self.space = all_together[:,2:-1:2] + 1j*all_together[:,3:-1:2]

        h = self.space[:,0]

        self.space = self.space[:,1:] / np.tile(  np.expand_dims(h,axis=1),  (1,self.space.shape[1]-1)  )




        