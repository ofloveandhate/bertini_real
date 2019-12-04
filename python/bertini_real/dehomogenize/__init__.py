def dehomogenize(points, index = 0):
    """ Dehomogenizes points
        Still need to finish this function, but this works when the
        dehomogenizing index is 0
        and the dimension is 1
    """
    
    new_points = []

    for i in range(1,len(points)):
        new_points.append(points[i] / points[index])
    return new_points
