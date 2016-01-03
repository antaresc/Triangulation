"""A utility package for random algorithms.

@author Antares Chen
@since  2016-1-2
"""

from random import randint, uniform
import math


def quick_select(points, compare):
    """Returns the median pair in points, the list of (x,y) pairs, based on the
    comparator function, fn, via the quickselect algorithm.

    Compare function takes in two arguments p1, p2 and returns a positive int
    if p1 > p2, 0 if p1 = p2, and a negative int if p1 < p2.
    """
    def select(points, position):
        pivot = points[len(points) // 2]

        less = [p for p in points if compare(p, pivot) < 0]
        equal = [p for p in points if compare(p, pivot) == 0]
        more = [p for p in points if compare(p, pivot) > 0]

        if position < len(less):
            return select(less, position)
        elif (position > len(less) + len(equal)):
            return select(more, position - len(less) - len(equal))
        else:
            return pivot
    return select(points, len(points) // 2)

def poisson_sample(size, r, length, width, resample = 30):
    """Returns a list of length size containing a set of randomly sampled points
    in 2D euclidean space.
    """
    min_dist = (int) (r / sqrt(2))
    grid = [[-1] * length for i in range(width)]

    pt = random_point(0, length, 0, width)
    grid[pt[0]][pt[1]] = 0
    index = 1
    active = [pt]
    result = [pt]

    while len(active) != 0:
        pt = active[randint(0, len(active) - 1)]




def random_point(x0, x1, y0, y1):
    """Generates a random point (x, y) where x and y is in [x0, x1] and [y0, y1]
    respectively.

    Assumes that x1 > x0 and y1 > y0
    """
    return randint(x0, x1), randint(y0, y1)

def random_point_around(x, y, r):
    """Returns a random point with at least distance r from (x, y)."""
    r1 = uniform(0, 1)
    r2 = uniform(0, 1)

    radius = (r1 + 1) * r
    angle = r2 * 2 * math.pi
    
    new_x = (int)(x + radius * math.cos(angle))
    new_y = (int)(y + radius * math.sin(angle))
    return new_x, new_y

def dist(x0, y0, x1, y1):
    """Returns the euclidean distance between (x0, y0) and (x1, y1)"""
    return math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
