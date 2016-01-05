"""A utility script for random algorithms.

Antares Chen
2016-1-2
"""

from random import randint, uniform
from sampler import poisson_sample
from PIL import Image, ImageDraw
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

def draw_points(points, length, width):
    """Plots the list of points on a canvas sized length x width."""
    image = Image.new('RGB', (length, width))
    draw = ImageDraw.Draw(image)
    draw.rectangle([(0, 0), (length, width)], fill = 'white')
    draw.point(points, fill = 'black')
    image.show()
