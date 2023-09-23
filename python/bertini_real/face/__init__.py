import bertini_real.cell


class Face(Cell):
	"""
	the one-dimensional cell type
	"""
	def __init__(self):
		super(Face, self).__init__(dimension=2)