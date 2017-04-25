#what does it mean to dehomogenize points
def dehomogenize(points):
    """ Dehomogenizes points
        Still need to finish this function, but this works when the dehomogenizing index is 0
        and the dimension is 1
    """
    index = 0
    new_points = []
    for i in xrange(len(points)):
        if i == index:
            continue
        else:
            new_points.append(points[i] / points[index])
    return new_points
