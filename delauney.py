"""A set of functions that implement delauney triangulation.

Two algorithms are implemented. First is the divide and conquer algorithm, while
second is the iterative algorithm, both provided by Guibas and Stolfi.

We implement the Quad-Edge data structure to allow delauney triangulation to
perform at optimal speeds.

Antares Chen
2016-1-3
"""

import numpy as np

def triangulate_divide_conquer():
    """Implements the divide and conquer schema proposed by Guibas-Stolfi."""
    pass

def triangulate_online():
    """Implements the iterative schema proposed by Guibas-Stolfi."""
    pass

def preprocess(points):
    """Sorts points by y-value and then sorts by x-value. Because the python's
    sorted is stable, it should maintain that all elements of points with the
    same x-value after being processed should be sorted by y-value as well.
    """
    points = sorted(points, key = lambda pt : pt[1])
    points = sorted(points, key = lambda pt : pt[0])
    return points


def right_of(point, quad_edge):
    """Returns if point is to the right of the given quad edge."""
    return is_ccw_circle(point, quad_edge.dest, quad_edge.orig)

def left_of(point, quad_edge):
    """Returns if point is to the left of the given quad edge."""
    return is_ccw_circle(point, quad_edge.orig, quad_edge.dest)

def valid(quad_edge, basel):
    """Returns if the quad edge is above the basel cross edge."""
    return right_of(quad_edge.dest, basel)

def is_ccw_circle(p1, p2, p3):
    """Returns if p1, p2, and p3, form a counterclockwise oriented circle."""
    test = np.array([[p1[0], p1[1], 1],
                        [p2[0], p2[1], 1],
                        [p3[0], p3[1], 1],
                    ])
    return np.linalg.det(test) > 0

def is_incircle(p1, p2, p3, d):
    """Returns if the point d is within the circumcircle defined by the points
    p1, p2, p3.

    The proof for this test is very interesting. The construction essentially
    shows that p1, p2, p3, and d when projected onto a standard parabolid
    x^2 + y^2, are coplanar if d is within the circumcircle struck by p1, p2,
    and p3. Refer to the original publication for full details.
    """
    test = np.array([[p1[0], p1[1], p1[0] ** 2 + p1[1] ** 2, 1],
                        [p2[0], p2[1], p2[0] ** 2 + p2[1] ** 2, 1],
                        [p3[0], p3[1], p3[0] ** 2 + p3[1] ** 2, 1],
                        [d[0], d[1], d[0] ** 2 + d[1] ** 2, 1]
                    ])
    return np.linalg.det(test) > 0


class QuadEdge:
    """A data structure that allows for fast planar subdivision. Holds three
    fields.

    _orig   holds a tuple representing the edge's point of origin on a manifold.
    _next   holds the next quad edge anchored about _orig counterclockwise.
    _rot    holds this edge's associated edge in the dual graph.

    All operators as described in the Guibas-Stolfi edge algebra are implemented
    as python properties. One can differentiate between the fields of the quad
    edge and the properties by the lack of the underscore in the latter.
    """

    def __init__(self, orig, next_edge, rot):
        """A basic constructor taking in orig, next_edge, and rot."""
        _orig = orig
        _next = next_edge
        _rot = rot

    def make_edge(orig, dest):
        """Creates the edge between orig and dest. Call this method instead of
        the constructor.
        """
        q0 = QuadEdge(orig, None, None)
        q1 = QuadEdge(None, None, None)
        q2 = QuadEdge(dest, None, None)
        q3 = QuadEdge(None, None, None)

        q0._rot = q1
        q1._rot = q2
        q2._rot = q3
        q3._rot = q0

        q0._next = q0
        q1._next = q3
        q2._next = q2
        q3._next = q1

        return q0

    def splice(q1, q2):
        """Splices Q1 with Q2. As defined in the original paper, splice is an
        operation such that the following holds.

        If two rings are distinct, splice will combine them into one. If two are
        exactly the same ring, then splice will break it into two separate
        pieces. If two rings are the same taken with opposite orientation, then
        splice will flip a segment of that ring.
        """
        dual1 = q1.orig_next._rot
        dual2 = q2.orig_next._rot

        q1.set_next(q2.orig_next)
        q2.set_next(q1.orig_next)

        dual1.set_next(dual2.orig_next)
        dual2.set_next(dual1.orig_next)

    def connect(q1, q2):
        """Returns a new quad edge connecting the destination of q1 to the
        origin of q2 maintaining that all three quad edges have the same face.
        """
        result = make_edge(q1.get_dest(), q2.get_orig())
        splice(result, q1.left_next)
        splice(result.sym, q2)
        return result

    def disconnect(q):
        """Disconnects q from the entire quad edge structure."""
        splice(e, e.orig_prev)
        splice(e.sym, e.sym.orig_prev)

    def swap(q):
        """Rotates Q counterclockwise given that Q is within an enclosing
        quadrilateral.
        """
        prev = q.orig_prev;
        sPrev = q.sym.orig_prev;

        splice(q, prev);
        splice(q.sym, sPrev);
        splice(q, prev.left_next);
        splice(q.sym, sPrev.left_next);

        q.set_orig(prev.get_dest());
        q.set_dest(sPrev.get_dest());

    def get_data(self):
        """Returns the data associated with this quad edge."""
        return self._data

    def get_orig(self):
        """Returns the origin of this quad edge."""
        return self._orig

    def get_dest(self):
        """Returns the destination of this quad edge."""
        return self.sym.get_orig()

    def set_data(self, data):
        """Sets the _data field to data."""
        self._data = data

    def set_orig(self, orig):
        """Sets the _orig field to orig."""
        self._orig = orig

    def set_dest(self, dest):
        """Sets the destination of this quad edge to dest."""
        self.sym.set_orig(dest)

    def set_next(self, next):
        """Sets the next quad edge to next."""
        self._next = next

    @property
    def sym(self):
        """Returns the quad edge with orig and dest reversed."""
        return _rot._rot

    @property
    def rot_inv(self):
        """Returns the dual of this quad edge flipped."""
        return self._rot.sym

    @property
    def orig_next(self):
        """Returns the next quad edge counterclockwise about the origin."""
        return self._next

    @property
    def orig_prev(self):
        """Returns the previous quad edge counterclockwise about the origin."""
        return self._rot.orig_next._rot

    @property
    def dest_next(self):
        """Returns the next quad edge counterclockwise about the destination."""
        return self.sym.orig_next.sym

    @property
    def dest_prev(self):
        """Returns the next quad edge counterclockwise about the destination."""
        return self.rot_inverse.orig_next.rot_inverse

    @property
    def left_next(self):
        """Returns the next quad edge counterclockwise about the left face."""
        return self.rot_inverse.orig_next._rot

    @property
    def left_prev(self):
        """Returns the previous quad edge counterclockwise about the left
        face.
        """
        return self._next.sym

    @property
    def right_next(self):
        """Returns the next quad edge counterclockwise about the right face."""
        return self._rot._next.rot_inverse

    @property
    def right_prev(self):
        """Returns the previous quad edge counterclockwise about the right
        face.
        """
        return self.sym.orig_next
