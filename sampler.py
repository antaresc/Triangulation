"""A set of functions that implement Bridson's fast poisson disk sampling in
two dimensions.

Antares Chen
2016-1-2
"""

from random import randint, uniform
import math


def poisson_sample(r, length, width, resample = 30):
    """Returns a list containing a set of randomly sampled points in 2D
    euclidean space bounded by length and width. Points are at least r distance
    apart
    """
    min_dist = (int) (r / math.sqrt(2))
    grid = [[-1] * (length + 1) for i in range(width + 1)]

    pt = random_point(0, length, 0, width)
    grid[pt[0]][pt[1]] = 0
    index = 1

    active = [pt]
    result = [pt]
    while len(active) != 0:
        pt = active.pop(randint(0, len(active) - 1))
        for i in range(resample):
            new_pt = random_point_around(pt[0], pt[1], r)
            if in_range(new_pt, length, width) and not in_neighborhood(new_pt, grid, r):
                active.append(new_pt)
                result.append(new_pt)
                grid[new_pt[0]][new_pt[1]] = index
                index += 1
    return result


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


def dist(p0, p1):
    """Returns the euclidean distance between (x0, y0) and (x1, y1)"""
    return math.sqrt((p1[0] - p0[0]) ** 2 + (p1[1] - p0[1]) ** 2)


def in_range(point, length, width):
    """Returns if the point is within range of length and width."""
    return point[0] >= 0 and point[0] < length and point[1] >= 0 and point[1] < width


def in_neighborhood(point, grid, min_dist):
    """Returns if any points are within min_dist of the given point in grid."""
    length = len(grid[0])
    width = len(grid)

    for r in range(point[0] - min_dist, point[0] + min_dist + 1):
        for c in range(point[1] - min_dist, point[1] + min_dist + 1):
            pt = (r, c)
            if in_range(pt, length, width) and grid[r][c] != -1:
                 if dist(point, pt) <= min_dist:
                     return True
    return False
