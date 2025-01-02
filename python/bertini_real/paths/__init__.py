
"""
    :platform: Unix, Windows, MacOS
    :synopsis: Utilities for plotting paths computed by bertini_real
"""


import numpy as np

def construct_matrix(lines):

    # read the size of the matrix
    m,n = [int(q) for q in lines.pop(0).split()]

    # preallocate
    mat = np.empty((m,n), dtype = np.complex64)

    # read one entry at a time
    for ii in range(m):
        for jj in range(n):
            a,b = lines.pop(0).split()
            c = float(a) + 1j*float(b)

            mat[ii,jj] = c

    return mat

class Path(object):
    """docstring for Path"""
    def __init__(self, filename, dehomogenize = False):

        super(Path, self).__init__()
        self.path = None
        self.num = -1

        self.time = None
        self.space = None

        self.space_der = None
        self.time_der = None

        self.dx_dt = None

        self.condition_number = None

        self.post_process = None

        self.num = int(filename.split('_')[1])
        self.is_dehomogenized = dehomogenize


        self._read(filename, dehomogenize)


    def _read(self, filename, dehomogenize):

        with open(filename,'r') as file:
            raw = file.read()

        the_path, the_postprocess = raw.split('-------------------------')
        
        self._process_the_path(the_path, dehomogenize)
        self._process_the_postprocess(the_postprocess)



    def _process_the_postprocess(self,the_postprocess):
        from collections import namedtuple

        PostProcess = namedtuple('PostProcess', ['condition_number', 
                                                 'function_residual', 
                                                 'newton_residual',
                                                 'time_at_final_sample_point',
                                                 'max_precision_used',
                                                 'time_first_precision_increase',
                                                 'accuracy_estimate_user_coord',
                                                 'accuracy_estimate_internal_coord',
                                                 'cycle_num',
                                                 'solution'])
        cond, res_f, res_n, t, max_p, time_first_p, acc_u, acc_i, cycle_num, soln = [0]*10

        as_lines = the_postprocess.strip().split('\n')

        for line in as_lines:
            if line.startswith('Cycle'):
                cycle_num = int(line.split(':')[1])

        self.post_process = PostProcess(cond, res_f, res_n, t, max_p, time_first_p, acc_u, acc_i, cycle_num, soln)


    def _process_the_path(self, the_path, dehomogenize):


        as_lines = the_path.strip().split('\n')

        path_num = int(as_lines.pop(0))
        assert path_num == self.num

        self.space_der = []
        self.time_der = []
        self.dx_dt = []


        path_as_numbers = []

        
        as_lines.pop(0) # it's an empty line

        while as_lines:

                path_as_numbers.append( np.array( [float(n) for n in as_lines.pop(0).strip().split()] ) )

                Jv = construct_matrix(as_lines)# these modify `as_lines` by popping
                Jp = construct_matrix(as_lines)# these modify `as_lines` by popping

                self.space_der.append(Jv) 
                self.time_der.append(Jp)

                from scipy.linalg import lu_factor, lu_solve
                
                lu, piv = lu_factor(Jv)
                dx_dt = -lu_solve((lu, piv), Jp).flatten()

                # dx_dt = -np.matmul(np.linalg.inv(Jv), Jp).flatten()  # see line 461, powerseries.hpp in Bertini 2

                self.dx_dt.append(dx_dt)


        # join all time values in one array
        all_together = np.array(path_as_numbers)

        # the time is the first (two) coordinate.  remember, time is complex in bertini
        self.time = all_together[:,0] + 1j*all_together[:,1]

        # i saved the condition number at every step
        self.condition_number = all_together[:,-1]
        self.space = all_together[:,2:-1:2] + 1j*all_together[:,3:-1:2]

        if dehomogenize:
            h = self.space[:,0] # the dehomogenizing coordinate

            # divide
            self.space = self.space[:,1:] / np.tile(  np.expand_dims(h,axis=1),  (1,self.space.shape[1]-1)  )






        