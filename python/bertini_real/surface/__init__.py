import bertini_real.parse

from bertini_real.decomposition import Decomposition
from bertini_real.curve import Curve


class Surface(Decomposition):

    def __init__(self, directory):
        self.directory = directory
        self.inputfilename = None
        self.num_variables = 0
        self.dimension = 9
        self.pi = []
        self.num_patches = 0
        self.patch = {}
        self.radius = 0
        self.center_size = 0
        self.center = []
        self.num_faces = 0
        self.num_midpoint_slices = 0
        self.num_critical_slices = 0
        self.num_singular_curves = 0
        self.singular_curve_multiplicities = []
        self.faces = {}    # stores all data from F.faces file
        self.midpoint_slices = []
        self.critical_point_slices = []
        self.critical_curve = []
        self.sphere_curve = []
        self.singular_curves = []
        self.singular_names = []
        self.surface_sampler_data = []   # store all surface_sampler data
        self.input = None

        # automatically parse data files to gather curve data
        self.parse_decomp(self.directory)
        self.parse_surf(self.directory)
        self.gather_faces(self.directory)
        self.gather_curves(self.directory)
        try:
            self.gather_surface_samples(self.directory)
        except:
            print("no samples found")

        self.read_input(self.directory)

    def __str__(self):
        result = "surface with:\n"
        result += "{} faces".format(self.num_faces)
        return result

    def parse_surf(self, directory):
        surf_data = bertini_real.parse.parse_Surf(directory)
        self.num_faces = surf_data[0]
        self.num_edges = surf_data[1]
        self.num_midpoint_slices = surf_data[2]
        self.num_critical_slices = surf_data[3]
        self.num_singular_curves = surf_data[4]
        self.singular_curve_multiplicities = surf_data[5]

    def gather_faces(self, directory):
        self.faces = bertini_real.parse.parse_Faces(directory)

    def gather_curves(self, directory):
        for ii in range(self.num_midpoint_slices):
            new_curve = Curve(directory + '/curve_midslice_' + str(ii))
            self.midpoint_slices.append(new_curve)
        for ii in range(self.num_critical_slices):
            new_curve = Curve(directory + '/curve_critslice_' + str(ii))
            self.critical_point_slices.append(new_curve)

        self.critical_curve = Curve(directory + '/curve_crit')
        self.sphere_curve = Curve(directory + '/curve_sphere')

        for ii in range(self.num_singular_curves):
            filename = directory + '/curve_singular_mult_' + \
                str(self.singular_curve_multiplicities[ii][0]) + '_' + str(
                    self.singular_curve_multiplicities[ii][1])
            new_curve = Curve(filename)
            self.singular_curves.append(new_curve)
            self.singular_names.append(new_curve.inputfilename)

    def gather_surface_samples(self, directory):
        self.surface_sampler_data = bertini_real.parse.parse_surface_Samples(
            directory)

    def check_data(self):
        try:
            if self.dimension != 2:
                print('This function designed to work on surfaces decomposed with bertini_real.  your object has dimension ' + self.dimension)

        except:
            return

    def faces_nonsingularly_connected(self, seed_index):
        self.check_data()

        new_indices = seed_index
        # connected = []
        connected = []

        while not(new_indices == []):
            # connected = concat([connected, new_indices])
            connected = [connected, seed_index]
            # print(connected)
            [new_indices] = self.find_connected_faces(connected)

        connected = connected.sort()
        set_num_faces = range(1, self.num_faces)

        # unconnected = setdiff(arange(1, self.num_faces), connected)
        unconnected = list(set(set_num_faces) - set(connected))

        return connected, unconnected

    def find_connected_faces(self, current):
        new_indices = []
        # unexamined_indices = arrange(1, self.num_faces)
        unexamined_indices = list(range(self.num_faces))

        # print(current)
        # print(unexamined_indices)

        # delete the faces we already know connect.
        # unexamined_indices[current] = []
        unexamined_indices.pop(0)
        current.pop(0)
        print(current)
        # print(unexamined_indices)

        for ii in range(len(current)):
            c = current[ii]
            f = self.faces[c]  # unpack the current face
            deleteme = []

            for jj in range(len(unexamined_indices)):
                d = unexamined_indices[jj]
                g = self.faces[d]  # unpack the examined face

                if self.faces_nonsingularly_connect(f, g):
                    new_indices.append(d)
                    deleteme.append(d)

            unexamined_indices = list(set(unexamined_indices) - set(deleteme))

        return new_indices

    def faces_nonsingularly_connect(self, f, g):
        val = False

        if self.cannot_possibly_meet(f, g):
            return  # assume no
        else:
            if self.meet_at_left(f, g):
                val = True
            else:
                if self.meet_at_right(f, g):
                    val = True
                else:
                    if self.meet_at_top(f, g):
                        val = True
                    else:
                        if self.meet_at_bottom(f, g):
                            val = True
        return val

    def cannot_possibly_meet(self, f, g):
        val = False

        if abs(f['middle slice index'] - g['middle slice index']) >= 2:
            val = True

        return val

    def meet_at_left(self, f, g):
        val = False

        for ii in range(f['num left']):
            e = self.critical_point_slices[
                f['middle slice index']].edges[f['left'][0]]
            a = e[2]

            for jj in range(g['num left']):
                E = self.critical_point_slices[
                    g['middle slice index']].edges[g['left'][0]]
                b = E[2]

                if a == b and not(self.is_degenerate(e)) and not(self.is_degenerate(E)):
                    val = True
                    return

            for jj in range(g['num right']):
                E = self.critical_point_slices[
                    g['middle slice index'] + 1].edges[g['right'][0]]
                b = E[2]

                if a == b and not(self.is_degenerate(e)) and not(self.is_degenerate(E)):
                    val = True
                    return
        return val

    def meet_at_right(self, f, g):
        val = False

        for ii in range(f['num right']):
            e = self.critical_point_slices[
                f['middle slice index'] + 1].edges[f['right'][0]]
            a = e[2]

            for jj in range(g['num left']):
                E = self.critical_point_slices[
                    g['middle slice index']].edges[g['left'][0]]
                b = E[2]

                if a == b and not(self.is_degenerate(e)) and not(self.is_degenerate(E)):
                    val = True
                    return

            for jj in range(g['num right']):
                E = self.critical_point_slices[
                    g['middle slice index'] + 1].edges[g['right'][0]]
                b = E[2]

                if a == b and not(self.is_degenerate(e)) and not(self.is_degenerate(E)):
                    val = True
                    return
        return val

    def meet_at_top(self, f, g):
        val = False

        if(f['system top'][1:15] == 'input_singcurve'):
            return  # cannot meet singularly, because edge is singular
        else:
            # at least they are in the same interval
            if f['middle slice index'] != g['middle slice index']:
                return

        if (f['system top'] == g['system top']):
            if (self.critical_curve.inputfilename == f['system top']):
                if (f['top'] == g['top']):
                    val = True
                    return

        if (f['system top'] == g['system bottom']):
            if (self.critical_curve.inputfilename == f['system top']):
                if (f['top'] == g['bottom']):
                    val = True
                    return

        return val

    def meet_at_bottom(self, f, g):
        val = False

        if(f['system bottom'][1:15] == 'input_singcurve'):
            return  # cannot meet singularly, because edge is singular
        else:
            # at least they are in the same interval
            if f['middle slice index'] != g['middle slice index']:
                return

        if (f['system bottom'] == g['system top']):
            if (self.critical_curve.inputfilename == f['system bottom']):
                if (f['bottom'] == g['top']):
                    val = True
                    return

        if (f['system bottom'] == g['system top']):
            if (self.critical_curve.inputfilename == f['system bottom']):
                if (f['bottom'] == g['top']):
                    val = True
                    return

        return val

    def is_degenerate(self, e):
        val = (e[0] == e[1]) or (e[1] == e[2])
        return val

    def separate_into_nonsingular_pieces(self):
        self.check_data()

        indices = []
        connected = []
        unconnected_this = [1]

        while not(unconnected_this == []):
            seed = unconnected_this[0]
            [connected_this, unconnected_this] = self.faces_nonsingularly_connected(
                seed)
            indices[end() + 1] = connected_this
            connected = [connected, connected_this]
            unconnected_this = setdiff(unconnected_this, connected)

        return indices


def separate_into_nonsingular_pieces(data=None, directory='Dir_Name'):
    surface = Surface(data)
    surface.separate_into_nonsingular_pieces()
