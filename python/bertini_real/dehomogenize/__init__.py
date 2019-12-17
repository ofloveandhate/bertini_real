"""

.. module:: demohogenize
    :platform: Unix, Windows
    :synopsis: The dehomogenize module makes points non-homogenized.

"""
def dehomogenize(points, index = 0):
    """ Dehomogenizes points, Still need to finish this function, but this works when the
        dehomogenizing index is 0 and the dimension is 1
        
        :param points: list of points
        :param index: zero index

        :rtype: Alist of new points
    """
    
    new_points = []

    for i in range(1,len(points)):
        new_points.append(points[i] / points[index])
    return new_points
