"""Basic 2D geometry primitives.
Including a BezierCurve class that implements a curve approximation
algorithm using biarcs.

Copyright (C) 2012 Claude Zervas, claude@utlco.com

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""
import logging
import math
logger = logging.getLogger(__name__)

DEBUG = True
# This can be set to a svg.SuperEffect during testing/debugging if using
# Inkscape as a visual/interactive test platform.
DEBUG_EFFECT = None
DEBUG_LAYER = None

TransformIdentity2D = [[1.0, 0.0, 0.0],[0.0, 1.0, 0.0]]
"""2D Tranform Identity Matrix"""

EPSILON = 0.00001
EPSILON2 = EPSILON**2
PHASH_RND = 5 # Rounding precision for point hash - should correspond to EPSILON
PHASH_STR = '%%.%df,%%.%df' % (PHASH_RND, PHASH_RND)

def set_epsilon(epsilon):
    """Set the global absolute error value and rounding limit for approximate
    floating point comparison operations. This value is accessible via the
    :attr:`planar.EPSILON` global variable.

    The default value of ``0.00001`` is suitable for values
    that are in the "countable range". You may need a larger
    epsilon when using large absolute values, and a smaller value
    for very small values close to zero. Otherwise approximate
    comparison operations will not behave as expected.
    """
    global EPSILON, EPSILON2, PHASH_RND, PHASH_STR
    EPSILON = float(epsilon)
    EPSILON2 = EPSILON**2
    precision = 0
    while epsilon < 1.0:
        precision += 1
        epsilon *= 10
    PHASH_RND = precision
    PHASH_STR = '%%.%df,%%.%df' % (PHASH_RND, PHASH_RND)

def float_eq(a, b):
    """Compare two floats for equality.
    The two float are considered equal if the difference between them is
    less than EPSILON.
    """
#    if a < EPSILON or b < EPSILON:
#        # for tiny values use a smaller epsilon
#        return abs(a - b) < EPSILON2
    return abs(a - b) < EPSILON


class P(tuple):
    """Two dimensional immutable point (vector).
    
    :param x: x coordinate.
    :type x: float
    :param y: y coordinate.
    :type y: float
    """
    def __new__(self, x, y=None):
        if y is None:
            return tuple.__new__(P, ((float(x[0]), float(x[1]))))
        else:
            return tuple.__new__(P, ((float(x), float(y))))

    @property
    def x(self):
        """The horizontal coordinate."""
        return self[0]

    @property
    def y(self):
        """The vertical coordinate."""
        return self[1]
    
    @staticmethod
    def from_polar(r, angle):
        """Create a Cartesian point from polar coordinates.
        See http://en.wikipedia.org/wiki/Polar_coordinate_system
        """
        x = r * math.cos(angle)
        y = r * math.sin(angle)
        return P(x,y)
    
    def to_polar(self):
        """Return the polar coordinates of this point as a tuple
        containing the length and angle respectively (r, a).
        """
        return (self.length(), self.angle())
        
    def is_zero(self):
        """ Return True if x and y are zero"""
        return self[0] == 0.0 and self[1] == 0.0
    
    def is_almost_zero(self):
        """Return True if x and y are within EPSILON distance to zero"""
        return self[0] + self[1] < EPSILON

    def almost_equal(self, other):
        """Compare vectors for approximate equality.
        :param other: Vector being compared.
        :type other: P
        :return: True if distance between the vectors < ``EPSILON``.
        """
        try:
            dx = self[0] - other[0]
            dy = self[1] - other[1]
            return (dx*dx + dy*dy) < EPSILON2
        except:
            return False

    def length(self):
        """The length or scalar magnitude of the vector."""
        return math.hypot(self[0], self[1])

    def length2(self):
        """The square of the length of the vector."""
        x, y = self
        return x*x + y*y

    def unit(self):
        """Return the vector scaled to unit length. If the vector
        is null, the null vector is returned.
        :rtype: P
        """
        L2 = self.length2()
        if L2 > EPSILON2:
            L = math.sqrt(L2)
            return P(self[0] / L, self[1] / L)
        else:
            return P(0.0, 0.0)
        
    def normal(self):
        """Return a vector perpendicular to this one."""
        return P(-self[1], self[0])

    def dot(self, other):
        """Compute the dot product with another vector.
        See http://en.wikipedia.org/wiki/Dot_product
        :param other: The vector with which to compute the dot product.
        :type other: P
        :rtype: float
        """
        x2, y2 = other
        return self[0] * x2 + self[1] * y2
    
    def cross(self, other):
        """Compute the cross product with another vector.
        Also called the perp-dot product for 2D vectors.
        Also called determinant for 2D matrix.
        See http://mathworld.wolfram.com/PerpDotProduct.html
        See http://www.gamedev.net/topic/441590-2d-cross-product/     
        :param other: The vector with which to compute the cross product.
        :type other: P
        :rtype: float
        """
        x2, y2 = other
        return self[0] * y2 - x2 * self[1]
    
    def angle(self):
        """Return the angle of this vector to x axis in radians."""
        return math.atan2(self[1], self[0])
    
    def angle2(self, p1, p2):
        """Return the angle between the two given vectors using this
        point as the origin.
        """
        v1 = (p1 - self)#.unit()
        v2 = (p2 - self)#.unit()
        #a = math.acos(v1.dot(v2))
        # Apparently this is more accurate for angles near 0 or PI:
        # see http://www.mathworks.com/matlabcentral/newsreader/view_thread/151925
        a = math.atan2(v1.cross(v2), v1.dot(v2))
        return a
            
    def distance(self, other):
        """Euclidean distance to other point."""
        x2, y2 = other
        return math.hypot(self[0] - x2, self[1] - y2)
        
    def distance2(self, other):
        """Euclidean distance squared to other point."""
        a = self[0] - other[0]
        b = self[1] - other[1]
        return (a * a) + (b * b)
        
    def distance_to_line(self, p1, p2):
        """Euclidean distance from this point to it's projection on a line
        that intersects the given points.
        See http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
        http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
        """
        v1 = p2 - p1 # Normalize the line segment
        seglen = v1.length()   # Segment length
        if seglen < EPSILON:       # Degenerate line segment...?
            return self.distance(p1) # TBD: This should probably be undefined...
        
        v2 = p1 - self
        return v1.cross(v2) / seglen
    
    def distance_to_segment(self, p1, p2):
        """Euclidean distance from this point to a line segment defined by
        the given start and end points.
        If this point's projection doesn't fall on the line segment then
        this returns the distance from this point to the closest segment
        endpoint.
        This method is useful for determining bezier curve flatness
        """
        seg = Line(p1, p2)
        return seg.distance_to_point(self, segment=True)
    
    def inside_triangle2D(self, A, B, C):
        """Return True if this point lies inside the triangle defined
        by points A, B, and C.
        Where ABC is clockwise or counter-clockwise.
        See http://www.sunshine2k.de/stuff/Java/PointInTriangle/PointInTriangle.html
        """
        # Using barycentric coordinates
        v1 = B - A
        v2 = C - A
        v3 = self - A
        det = v1.cross(v2)
        s = v1.cross(v3) / det
        t = v2.cross(v3) / det
        return s >= 0.0 and t >= 0.0 and (s+t) <= 1

    def __eq__(self, other):
        """Compare for equality."""
        # Use EPSILON to compare point values so that spatial hash tables
        # and other geometric comparisons work as expected.
        # There may be cases where an exact compare is necessary but for
        # most purposes (like collision detection) this works better.
        return self.almost_equal(other)
#        try:
#            x2, y2 = other
#            return self[0] == x2 and self[1] == y2
#        except:
#            return False

    def __nonzero__(self):
        """Return True if this is not a null vector (x!=0 and y!=0)."""
        return self[0] > EPSILON or self[1] > EPSILON

    def __neg__(self):
        """Return the unary negation of the vector."""
        return P(-self[0], -self[1])
    
    def __add__(self, other):
        """Add the vector to a scalar or another vector.

        :param other: The vector or scalar to add.
        :type other: P
        """
        try:
            n = float(other)
            return P(self[0] + n, self[1] + n)
        except TypeError:
            try:
                x2, y2 = other
                return P(self[0] + x2, self[1] + y2)
            except Exception:
                return NotImplemented

    __iadd__ = __add__

    def __sub__(self, other):
        """Subtract a scalar or another vector from this vector.

        :param other: The vector or scalar to substract.
        :type other: P
        """
        try:
            n = float(other)
            return P(self[0] - n, self[1] - n)
        except TypeError:
            try:
                x2, y2 = other
                return P(self[0] - x2, self[1] - y2)
            except Exception:
                return NotImplemented

    __isub__ = __sub__

    def __mul__(self, other):
        """Multiply the vector by a scalar. This operation is undefined
        for any other type since it doesn't make geometric sense. Use dot()
        or cross() instead.

        :param other: The scalar to multiply by.
        :type other: float
        """
        try:
            n = float(other)
        except TypeError:
            return NotImplemented
        else:
            return P(self[0] * n, self[1] * n)
    
    __rmul__ = __imul__ = __mul__
        
    def __truediv__(self, other):
        """Divide the vector by a scalar.

        :param other: The value to divide by.
        :type other: float
        """
        try:
            other = float(other)
            return tuple.__new__(P, (self[0] / other, self[1] / other))
        except TypeError:
            return NotImplemented

    __idiv__ = __div__ = __itruediv__ = __truediv__

    def __floordiv__(self, other):
        """Divide the vector by a scalar or componentwise by
        another vector, rounding down.

        :param other: The value to divide by.
        :type other: Vec2 or float
        """
        try:
            other = float(other)
        except TypeError:
            return NotImplemented
        else:
            return tuple.__new__(P, (self[0] // other, self[1] // other))

    __ifloordiv__ = __floordiv__

    def __pos__(self):
        return self

    def __abs__(self):
        """Compute the absolute magnitude of the vector."""
        return self.length()

    def __str__(self):
        """Concise string representation."""
        return "P(%.4f, %.4f)" % self

    def __repr__(self):
        """Precise string representation."""
        return "P(%r, %r)" % self
    
    def __hash__(self):
#        hashval = (PHASH_STR % (self[0], self[1])).__hash__()
        # This is from http://www.beosil.com/download/CollisionDetectionHashing_VMV03.pdf
        # TODO: consider http://www.boost.org/doc/libs/1_46_1/boost/functional/hash/detail/hash_float_generic.hpp
        # This seems to work pretty well and is very fast.
        a = int(self[0] * 73856093)
        b = int(self[1] * 83492791)
        hashval = a ^ b
        return hashval
    
#    __hash__ = tuple.__hash__ # hash is not inherited in Py 3
    
    def SVG_plot(self, color='#000000'):
        """Draw a dot. Useful for debugging and testing."""
        if DEBUG_EFFECT is not None and DEBUG_LAYER is not None:
            style = 'fill:%s;stroke:none' % (color,)
            DEBUG_EFFECT.create_circle(self.x, self.y, '.75pt', style, DEBUG_LAYER)
        

# Make some method aliases to be compatible with various Point implementations
P.mag = P.length
P.normalized = P.unit
P.perpendicular = P.normal

class Shape(tuple):
    """Base class of the geometric shapes.
    This just defines some virtual methods that subclasses
    should implement."""
    def start_tangent(self):
        """Return the angle in radians of a line tangent to this shape
        beginning at the first point."""
        return 0.0

class Line(Shape):
    """Two dimensional immutable line segment defined by two points.
    
    :param p1: start point.
    :type p1: P
    :param p2: end point.
    :type p2: P
    """
    def __new__(self, p1, p2):
        return tuple.__new__(Line, (P(p1), P(p2)))

    @property
    def p1(self):
        """The start point of line."""
        return self[0]

    @property
    def p2(self):
        """The end point of line."""
        return self[1]
    
    def length(self):
        """Return the length of this line segment."""
        return (self.p2 - self.p1).length()
    
    def slope(self):
        """Return the slope of this line."""
        return (self.p2.x - self.p1.x) / (self.p2.y - self.p1.y)
    
    def angle(self):
        """Return the angle of this line segment in radians."""
        return (self.p2 - self.p1).angle()
    
    def start_tangent(self):
        """Return the direction of this line segment in radians.
        This is the same as the angle"""
        return self.angle()
    
    def midpoint(self):
        """Return the midpoint of the line segment defined by p1 and p2"""
        #return P((self.p1.x + self.p2.x) / 2, (self.p1.y + self.p2.y) / 2)
        return (self.p1 + self.p2) * 0.5
    
    def bisector(self):
        """Return a line that is perpendicular to and passes through
        the midpoint of this line. Also called the perpendicular bisector.
        Essentially this line segment is rotated 90deg about its midpoint.
        """
        midp = self.midpoint()
        p1 = self.p1 - midp
        p2 = self.p2 - midp
        bp1 = midp + P(p1.y, -p1.x)
        bp2 = midp + P(p2.y, -p2.x)
        return Line(bp1, bp2)

    def distance_to_point(self, p, segment=False):
        """Return the Euclidian distance from the spcified point and
        its projection on to this line.
        If <segment> is True then if the point projection does not lie
        within the two points that define this line segment return the
        shortest distance to either of the two endpoints.
        See http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
        http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
        """
        v1 = self.p2 - self.p1 # Normalize the line segment
        v2 = p - self.p1 # Normalize the point vector
        L2 = v1.length2()
        if L2 < EPSILON2: # Segment points are coincident?
            return self.p1.distance(p)
        u = v2.dot(v1) / v1.length2() # Projection (0->1) on to segment
        if u < 0: # Projection not on segment but nearer to p1?
            return self.p1.distance(p)
        elif u > 1.0: # Projection not on segment but nearer to p2?
            return self.p2.distance(p)
        p_proj = self.p1 + v1 * u # Point of projection on line segment
        d = p.distance(p_proj) # distance between point and projection
        return d
        
    def intersection(self, other, segment=False):
        """Return the intersection point (if any) of this line and another line.
        If <segment> is True then the intersection point must lie on both
        segments.
        Returns a point if it intersects otherwise None.
        See http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/
        and http://mathworld.wolfram.com/Line-LineIntersection.html
        """
        x1, y1 = self[0]
        x2, y2 = self[1]
        x3, y3 = other[0]
        x4, y4 = other[1]
        
        a = (x4-x3) * (y1-y3) - (y4-y3) * (x1-x3)
        b = (x2-x1) * (y1-y3) - (y2-y1) * (x1-x3)
        denom = (y4-y3) * (x2-x1) - (x4-x3) * (y2-y1)
        
        if abs(denom) < EPSILON: # Lines are parallel ?
            if abs(a) < EPSILON and abs(b) < EPSILON: # Lines are coincident ?
                return self.midpoint()
            else:
                return None
        
        mua = a / denom
        mub = b / denom
        if segment and (0. > mua or mua > 1. or 0. > mub or mub > 1.):
            # The intersection lies outside the line segments
            return None
        x = x1 + mua * (x2 - x1)
        y = y1 + mua * (y2 - y1)
        
        return P(x, y)

    def same_side(self, pt1, pt2):
        """Return True if the given points lie on the same side of this line.
        """
        # Normalize the points first
        v1 = self.p2 - self.p1
        v2 = pt1 - self.p1
        v3 = pt2 - self.p1
        # The sign of the perp-dot product determines which side the point lies.
        c1 = v1.cross(v2)
        c2 = v1.cross(v3)
        return (c1 >=0 and c1 >= 0) or (c1 < 0 and c2 < 0)
    
    def reverse(self):
        """Return a Line segment with start and end points reversed."""
        return Line(self.p2, self.p1)
        
    def __eq__(self, other):
        """Compare for segment equality in a geometric sense.
        Useful for culling overlapping segments..."""
        # Compare both directions
        same = (self.p1.almost_equal(other[0]) and self.p2.almost_equal(other[1])) \
            or (self.p1.almost_equal(other[1]) and self.p2.almost_equal(other[0]))
        logger.debug('segment: %s == %s' % (str(self), str(other)))
        if same:
            logger.debug('segment same')
        return same

    def __hash__(self):
        """Create a hash value for this line segment.
        The hash value will be the same if p1 and p2 are reversed.
        """
        return hash(self.p1) ^ hash(self.p2)

    def __str__(self):
        """Concise string representation."""
        return "Line(P(%.4f, %.4f), P(%.4f, %.4f))" % (self.p1[0], self.p1[1], self.p2[0], self.p2[1])

    def __repr__(self):
        """Precise string representation."""
        return "Line(P(%r, %r), P(%r, %r))" % (self.p1[0], self.p1[1], self.p2[0], self.p2[1])

    def to_SVG_path(self):
        """Return a string with the SVG path 'd' attribute
        that corresponds with this line.
        """
        return 'M %f %f L %f %f' % (self.p1.x, self.p1.y, self.p2.x, self.p2.y)
    
    def SVG_plot(self, color='#cccc99'):
        """Draw an SVG line for debugging/testing"""
        if DEBUG_EFFECT is not None and DEBUG_LAYER is not None:
            style = 'fill:none;stroke:%s;stroke-width:0.25pt;stroke-opacity:1' % (color,)
            attrs = {'d': self.to_SVG_path(), 'style': style}
            DEBUG_EFFECT.create_path(attrs, layer=DEBUG_LAYER)
        

class Rectangle(Shape):
    """Two dimensional immutable rectangle defined by two points
    and an optional angle of rotation.
    The first point is the lower left corner,
    the second is the upper right.
    """
    def __new__(self, p1, p2, angle=0.0):
        # Canonicalize the point order so that p1 is
        # always lower left.
        if p1[0] <= p2[0]:
            x1 = p1[0]
            x2 = p2[0]
        else:
            x1 = p2[0]
            x2 = p1[0]
        if p1[1] <= p2[1]:
            y1 = p1[1]
            y2 = p2[1]
        else:
            y1 = p2[1]
            y2 = p1[1]
        return tuple.__new__(Rectangle, (P(x1, y1), P(x2, y2), angle))
    
    @property
    def p1(self):
        """The first corner of rectangle."""
        return self[0]

    @property
    def p2(self):
        """The second corner of rectangle."""
        return self[1]
    
    def height(self):
        """Height of rectangle. (along Y axis)"""
        return self[1][1] - self[0][1]
    
    def width(self):
        """Width of rectangle. (along X axis)"""
        return self[1][0] - self[0][0]
    
    def point_inside(self, p):
        """Return True if the point is inside this rectangle."""
        return p[0] > self[0][0] and p[0] < self[1][0] \
                and p[1] > self[0][1] and p[1] < self[1][1]

    def line_inside(self, line):
        """Return True if the line segment is inside this rectangle."""
        return self.point_inside(line.p1) and self.point_inside(line.p2)
    
    def all_points_inside(self, points):
        """Return True if the given set of points lie inside this rectangle."""
        for p in points:
            if not self.point_inside(p):
                return False
        return True

    def start_tangent(self):
        """Return the angle in radians of a line tangent to this shape
        beginning at the first point.
        The corner point order for rectangles is clockwise from lower left.
        """
        return self[2] + math.pi/2
    

class Arc(Shape):
    """Two dimensional immutable circular arc segment.
    
    :param p1: start point.
    :type p1: P
    :param p2: end point.
    :type p2: P
    """
    def __new__(self, p1, p2, radius, angle, center=None):
        if center is None:
            center = Arc._center(p1, p2, radius, angle)
        return tuple.__new__(Arc, (P(p1), P(p2), radius, angle, P(center)))
    
    @staticmethod
    def from_two_points_and_tangent(p1, ptan, p2, reverse=False):
        """Create an Arc given two points and a tangent vector from p1->ptan.
        Reverse the resulting arc direction if <reverse> is True.
        """
        # The arc angle is 2 * the angle defined by the tangent and the secant.
        # See http://en.wikipedia.org/wiki/Tangent_lines_to_circles
        angle = 2 * p1.angle2(ptan, p2)
        chord_len = p1.distance(p2)
        if chord_len < EPSILON or p1.distance(ptan) < EPSILON:
            logging.debug('degenerate arc: angle=%.5f, p1=(%.5f, %.5f), t=(%.4f,%.4f), p2=(%.4f,%.4f), reverse=%s' % \
                          (angle, p1.x, p1.y, ptan.x, ptan.y, p2.x, p2.y, str(reverse)))
            return None
        radius = abs(chord_len / (2 * math.sin(angle / 2)))
        if reverse:
            return Arc(p2, p1, radius, -angle)
        else:
            return Arc(p1, p2, radius, angle)
         
    @property
    def p1(self):
        """The start point of the arc."""
        return self[0]

    @property
    def p2(self):
        """The end point of the arc."""
        return self[1]
    
    @property
    def radius(self):
        """The radius of the arc."""
        return self[2]

    @property
    def angle(self):
        """Return the angle of this arc. The sign of the angle
        determines its direction.
        """
        return self[3]
    
    @property
    def center(self):
        """Return the center point of this arc."""
        return self[4]
    
    @staticmethod
    def _center(p1, p2, radius, angle):
        """Calculate the center point of an arc given two endpoints,
        the radius, and an angle.
        This method is static so that it can be used by __new__.
        http://math.stackexchange.com/questions/27535/how-to-find-center-of-an-arc-given-start-point-end-point-radius-and-arc-direc
        """
        if p1.almost_equal(p2): # Points coinciedent?
            return p1
        chord = Line(p1, p2)
        # distance between start and endpoint
        d = chord.length()
        # determine mid-point
        mp = chord.midpoint()
        # Unit normal of midpoint
        #n = P(-((p2.y - p1.y) / d), (p2.x - p1.x) / d)
        # distance from center to midpoint
        h = math.sqrt(radius**2 - (d**2 / 4))
        # Determine which side the arc center is
        sign = 1 if angle > 0 else -1
        # calculate the center point
        cx = mp.x - (sign * h * ((p2.y - p1.y) / d))
        cy = mp.y + (sign * h * ((p2.x - p1.x) / d))
        return P(cx, cy)
    
    def length(self):
        """Return the length of this arc segment (radius * angle)"""
        return abs(self.radius * self.angle)
    
    def area(self):
        """Return the area between the arc and its circular center"""
        r2 = self.radius * self.radius
        return r2 * abs(self.angle) / 2
    
    def segment_area(self):
        """Return the area of the shape limited by the arc and a straight line
        between the two end points"""
        r2 = self.radius * self.radius
        return r2 * abs(self.angle - math.sin(self.angle)) / 2
    
    def is_clockwise(self):
        """Return True if this arc is clockwise."""
        return self.angle < 0
    
    def start_tangent(self):
        """Return the direction of this arc segment in radians.
        This is the angle of a tangent vector at the arc segment's
        first point."""
        v = self.center - self.p1 if self.is_clockwise() else self.p1 - self.center
        return v.angle() + math.pi/2
        
    def distance_to_point(self, p, inside_only=True):
        """Return the Hausdorff distance from this arc segment
        to the specified point.
        :inside_only: return -1 if the point->circle projection
        is outside the arc segment
        """
        # Check for degenerate arc case
        if self.radius < EPSILON or abs(self.angle) < EPSILON:
            return self.p1.distance(p)
        # If the point->circle projection is outside the arc segment
        # then return the distance closest to either endpoint.
        # Note: see http://www.blackpawn.com/texts/pointinpoly/default.html
        # http://www.sunshine2k.de/stuff/Java/PointInTriangle/PointInTriangle.html
        # http://blogs.msdn.com/b/rezanour/archive/2011/08/07/barycentric-coordinates-and-point-in-triangle-tests.aspx
        # Using barycentric coordinates
        v1 = self.p1 - self.center
        v2 = self.p2 - self.center
        v3 = p - self.center
        det = v1.cross(v2)
        s = v1.cross(v3) / det
        t = v3.cross(v2) / det
        is_inside_arc = (s >= 0.0 and t >= 0.0)
        if is_inside_arc:
            #Line(self.center, p).SVG_plot('#00cccc')
            # Distance from arc center to point.
            dp = self.center.distance(p)
            # Distance from point to edge of arc.
            d = abs(dp - self.radius)
#            Line(p, (p - self.center) * (d / dp) + p).SVG_plot('#0000cc')
        elif inside_only:
            return -1
        else:
#            Line(self.center, p).SVG_plot('#c0cc00')
            # Otherwise distance to closest arc segment endpoint.
            d = min(self.p1.distance(p), self.p2.distance(p))
        return d
    
    def point_at(self, angle):
        """Return the point on this arc given the specified angle from
        the start point of the arc segment."""
        if angle < EPSILON or angle > abs(self.angle):
            return None
        # TODO: there is surely a faster way to do this....
        v1 = self.p1 - self.center
        v2 = self.p2 - self.center
        a1 = math.fmod(math.atan2(v1.y, v1.x) + math.pi, 2 * math.pi)
        a2 = math.fmod(math.atan2(v2.y, v2.x) + math.pi, 2 * math.pi)
        if a1 > a2:
            angle = a1 - angle
        else:
            angle = a1 + angle
        x = self.center.x + self.radius * math.cos(angle)
        y = self.center.y + self.radius * math.sin(angle)
        return P(x, y)
            
    def subdivide(self, angle):
        """Split this arc into two arcs at the point on this arc given
        by the specified positive arc angle (0-2pi) from the start point.
        Returns a tuple containing one or two Arc objects.
        """
        if angle < EPSILON:
            return (self,)
        angle2 = abs(self.angle) - angle
        p = self.point_at(angle)
        if self.angle < 0:
            angle = -angle
            angle2 = -angle2
        arc1 = Arc(self.p1, p, self.radius, angle, self.center)
        arc2 = Arc(p, self.p2, self.radius, angle2, self.center)
        return (arc1, arc2)

    def to_SVG_path(self):
        """Return a string with the SVG path 'd' attribute
        that corresponds to this arc.
        """
        sweep_flag = 0 if self.angle < 0 else 1
        return 'M %f %f A %f %f 0.0 0 %d %f %f' % \
                (self.p1.x, self.p1.y, self.radius,
                 self.radius, sweep_flag, self.p2.x, self.p2.y)

    def SVG_plot(self, color='#cccc99'):
        """Draw an SVG arc for debugging/testing"""
        if DEBUG_EFFECT is not None and DEBUG_LAYER is not None:
            style = 'fill:none;stroke:%s;stroke-width:0.25pt;stroke-opacity:1' % (color,)
            attrs = {'d': self.to_SVG_path(), 'style': style}
            DEBUG_EFFECT.create_path(attrs, layer=DEBUG_LAYER)
            # Draw the center-arc wedge
            seg1 = Line(self.center, self.p1)
            seg2 = Line(self.center, self.p2)
            self.center.SVG_plot(color=color)
            seg1.SVG_plot(color=color)
            seg2.SVG_plot(color=color)
                

class Circle(Shape):
    """Two dimensional immutable circle.
    """
    def __new__(self, center, radius):
        return tuple.__new__(Circle, (P(center), radius))
    
    @staticmethod
    def from_three_points(p1, p2, p3):
        """Create a Circle given three points on its circumference."""
        raise Exception('Unimplemented')

    @property
    def center(self):
        """Return the center point of this circle."""
        return self[0]

    @property
    def radius(self):
        """The radius of this circle."""
        return self[1]

    def point_inside(self, p):
        """Return True if the given point is inside this circle."""
        return self.center.distance2(p) < self.radius**2

    def all_points_inside(self, points):
        """Return True if all the given points are inside this circle."""
        for p in points:
            if not self.point_inside(p):
                return False
        return True

class CubicBezier(Shape):
    """Two dimensional immutable cubic bezier curve.
    
    :param p1: start point.
    :type p1: P
    :param c1: start control point.
    :type c1: P
    :param c2: end control point.
    :type c2: P
    :param p2: end point.
    :type p2: P
    """
    def __new__(self, p1, c1, c2, p2):
        return tuple.__new__(CubicBezier, (P(p1), P(c1), P(c2), P(p2)))

    @property
    def p1(self):
        """The start point of curve."""
        return self[0]

    @property
    def c1(self):
        """The start control point of curve."""
        return self[1]

    @property
    def c2(self):
        """The end control point of curve."""
        return self[2]
    
    @property
    def p2(self):
        """The end point of curve."""
        return self[3]
    
    def start_tangent(self):
        """Return the tangent direction of this curve in radians.
        This is usually the angle of the first control point vector."""
        return self.tangent_at(0.0).angle()

    def point_at(self, t):
        """A point on the curve corresponding to <t>."""
        return (self.p1 * (1 - t)**3 +
                self.c1 * t * 3 * (1 - t)**2 +
                self.c2 * t**2 * 3 * (1-t) +
                self.p2 * t**3)
        
    def tangent_at(self, t):
        """The tangent vector at the point on the curve corresponding to <t>."""
        return ((self.c1 - self.p1) * 3 * (1 - t)**2 +
                (self.c2 - self.c1) * 6 * t * (1 - t) +
                (self.p2 - self.c2) * 3 * t**2)
        
    def is_straight_line(self, flatness=EPSILON):
        """Return True if curve is essentially a straight line
        :param flatness: The required flatness to be considered a line
            (default is EPSILON)
        :type flatness: float
        """
        return self.flatness() < flatness
    
    def flatness(self):
        """Return the flatness of this curve.
        The maximum distance between the control points and the line segment
        defined by the start and end points of the curve.
        """
        chord = Line(self.p1, self.p2)
        d1 = chord.distance_to_point(self.c1, segment=True)
        d2 = chord.distance_to_point(self.c2, segment=True)
        flatness = max(d1, d2)
        return flatness
    
    def subdivide(self, t):
        """Subdivide this curve at the point corresponding to <t>
        into two cubic bezier curves, where 0<=t<=1.
        Uses De Casteljaus's algorithm.
        Returns a tuple of one or two CubicBezier objects.
        """
        if t <= EPSILON or t >= 1.0:
            return (self,)
        cp0, cp1, p, cp2, cp3 = self.controlpoints_at(t)
#        p.SVG_plot('#00ffff')
        curve1 = CubicBezier(self.p1, cp0, cp1, p)
        curve2 = CubicBezier(p, cp2, cp3, self.p2)
#        curve1.SVG_plot()
#        curve2.SVG_plot()
        return (curve1, curve2)
    
    def subdivide_inflections(self):
        """Subdivide this curve at the inflection points, if any.
        Returns a tuple containing one to three curves depending on whether
        there are no inflections, one inflection, or two inflections.
        """
        t1, t2 = self.find_inflections()
        if t1 > 0.0 and t2 == 0.0:
            return self.subdivide(t1)
        elif t1 == 0.0 and t2 > 0.0:
            return self.subdivide(t2)
        elif t1 > 0.0 and t2 > 0.0:
            curves = []
            curve1, curve2 = self.subdivide(t1)
            # Two inflection points
            t1, t2 = curve1.find_inflections()
            t = max(t1, t2)
            if t > 0.0:
                curves.extend(curve1.subdivide(t))
                curves.append(curve2)
            else:
                t1, t2 = curve1.find_inflections()
                t = max(t1, t2)
                if t > 0.0:
                    curves.append(curve1)
                    curves.extend(curve2.subdivide(t))
                else:
                    curves.extend((curve1, curve2))
            #curves[0].SVG_plot()
            #curves[1].SVG_plot()
            #curves[2].SVG_plot()
            return curves
        else:
            return (self,)
    
    def controlpoints_at(self, t):
        """Get the point on this curve corresponding to <t> plus control points.
        Returns a tuple of the form (C0, C1, P, C2, C3).
        Relevant points found using De Casteljaus's algorithm. Useful for
        subdividing the curve at <t>.
        """
        # First intermediate points
        d01 = (1. - t) * self.p1 + t * self.c1
        d12 = (1. - t) * self.c1 + t * self.c2
        d23 = (1. - t) * self.c2 + t * self.p2
        # Second intermediate points
        d012 = (1. - t) * d01 + t * d12
        d123 = (1. - t) * d12 + t * d23
        # Finally, the split point
        d0123 = (1. - t) * d012 + t * d123
        
        return (d01, d012, d0123, d123, d23)
    
    def derivative1(self, t):
        """Get the 1st derivative of this curve at <t>.
        Returns a vector (point).
        """
        t2 = t**2
        return (self.p1 * ((2 * t - t2 - 1) * 3) +
                self.c1 * ((3 * t2 - 4 * t + 1) * 3) +
                self.c2 * (t * (2 - 3 * t) * 3) +
                self.p2 * (t2 * 3) )    
    
    def derivative2(self, t):
        """Get the 2nd derivative of this curve at <t>.
        Returns a vector (point).
        """
        return 6 * ( (1 - t) * self.p1 + (3*t - 2) * self.c1 +
                     (1 - 3*t) * self.c2 + t * self.p2 )
    
    def find_inflections(self):
        """Find <t1, t2> where the curve has C1 discontinuities
        (ie from convex to concave or vice versa, or a cusp, or a loop).
        There may be none, one, or two inflections on the curve.
        See http://www.caffeineowl.com/graphics/2d/vectorial/cubic-inflexion.html
        Returns a tuple containing the two inflection locations.
        The location values will be 0 if no inflection.
        """
        # Basically the equation to be solved is
        # P' * P'' = 0
        # Where P' and P'' are the first and second derivatives respectively
        
        # Temporary vectors to simplify the math
        v1 = self.c1 - self.p1
        v2 = self.c2 - self.c1 - v1
        v3 = self.p2 - self.c2 - v1 - 2 * v2
        
        # Calculate quadratic coefficients
        # of the form a*x**2 + b*x + c = 0
        a = v2.x * v3.y - v2.y * v3.x
        b = v1.x * v3.y - v1.y * v3.x
        c = v1.x * v2.y - v1.y * v2.x
        
        # Get the two roots of the quadratic - if any
        # These will be the inflection locations if 0 >= t <= 1
        t1 = 0.0
        t2 = 0.0
        D = b * b - 4 * a * c # the discriminant
        aa = 2 * a
        if abs(aa) > 0.0: # Avoid div by zero
            d = math.sqrt(abs(D))
            t1 = (-b - d) / aa
            t2 = (-b + d) / aa
    
            if t1 < EPSILON or t1 >= (1.0 - EPSILON):
                t1 = 0.0
            if t2 < EPSILON or t2 >= (1.0 - EPSILON):
                t2 = 0.0
            
        return (t1, t2)
    
    def biarc_joint_arc(self):
        """Calculate the arc that intersects the two endpoints of this curve
        and the set of possible biarc joints.
        """
        # The center <s> of the circle is the intersection of the bisectors
        # of line segments P1->P2 and (P1+unit(C1))->(P2+unit(C2))
        chord = Line(self.p1, self.p2)
        u1 = (self.c1 - self.p1).unit()
        u2 = (self.c2 - self.p2).unit()
        useg = Line(self.p1 + u1, self.p2 + P(-u2.x, -u2.y))
        center = chord.bisector().intersection(useg.bisector())
        radius = center.distance(self.p1)
        angle = center.angle2(self.p1, self.p2)
        return Arc(self.p1, self.p2, radius, angle, center)
    
    def biarc_approximation(self, tolerance=0.01, max_depth=4, line_flatness=0.01, _recurs_depth=0):
        """Approximate this curve using biarcs.
        This will recursively subdivide the curve into a series of
        G1 connected arcs or lines until the Hausdorff distance between the
        approximation and this bezier curve is within the specified tolerance.
        Returns a list of Arc and/or Line objects.
        """
#        self.SVG_plot() #DEBUG
        # Check for degenerate cases:
        # Bail if the curve endpoints are coincident.
        if self.p1.almost_equal(self.p2):
            return ()                
        # Or if the curve is basically a straight line then return a Line.
        if self.flatness() < line_flatness:
            return (Line(self.p1, self.p2),)
        
        if _recurs_depth == 0:
            # Subdivide this curve at any inflection points to make sure
            # the curve has monotone curvature with no discontinuities.
            # This should only be required once.
            curves = self.subdivide_inflections()
            if len(curves) > 1:
                biarcs = []        
                for curve in curves:
                    biarcs.extend(curve.biarc_approximation(tolerance=tolerance, max_depth=max_depth, line_flatness=line_flatness, _recurs_depth=_recurs_depth+1))
                return biarcs

        # Calculate the arc that intersects the two endpoints of this curve
        # and the set of possible biarc joints.
        j_arc = self.biarc_joint_arc()
#        j_arc.SVG_plot(color='#c00000') #DEBUG
        # Another degenerate case
        if j_arc.radius < EPSILON:
            return ()
        
        # TODO: Use a better method of finding J.
        # A possibly more accurate method (see A. Riskus, 2006)
        # is to find the intersection of the joint arc and this curve
        # and use that point for J. The curve subdivision would occur at J
        # as well.
        # To make this simple for now:
        # The biarc joint J will be the intersection of the line
        # whose endpoints are the center of the joint arc and the
        # 'middle' (t=0.5) of the bezier curve and the joint arc.
        p = self.point_at(0.5)
        v = p - j_arc.center
        pjoint = v * (j_arc.radius / v.length()) + j_arc.center
#        pjoint.SVG_plot(color='#ffff00') #DEBUG
                
        # Create the two arcs that define the biarc. These will be None
        # if the control points are degenerate (ie collinear).
        arc1 = Arc.from_two_points_and_tangent(self.p1, self.c1, pjoint)
        arc2 = Arc.from_two_points_and_tangent(self.p2, self.c2, pjoint, reverse=True)
        if arc1 is None or arc2 is None:
            # In case of degenerate control points just create line segements.
            # In general this should be caught by the test for curve flatness.
            return (arc1 if arc1 is not None else Line(self.p1, pjoint),
                    arc2 if arc2 is not None else Line(pjoint, self.p2))
        
        if _recurs_depth < max_depth:
            # Calculate Hausdorff distances from arcs to this curve
            d1 = self.distance_to_arc(arc1)
            d2 = self.distance_to_arc(arc2)
            # If Hausdorff distance is above tolerance then split this curve and
            # recursively approximate the two sub-curves.
            if d1 > tolerance or d2 > tolerance:
                curve1, curve2 = self.subdivide(0.5)
                biarcs = []
                biarcs.extend(curve1.biarc_approximation(tolerance, max_depth=max_depth, line_flatness=line_flatness, _recurs_depth=_recurs_depth+1))
                biarcs.extend(curve2.biarc_approximation(tolerance, max_depth=max_depth, line_flatness=line_flatness, _recurs_depth=_recurs_depth+1))
                return biarcs
#        arc1.SVG_plot(color='#00c000') #DEBUG
#        arc2.SVG_plot(color='#00c000') #DEBUG
        return (arc1, arc2)

    def distance_to_arc(self, arc, ndiv=9):
        """Calculate an approximate Hausdorff distance between this curve and
        the specified circular arc.
        The approximation accuracy depends on the number of curve subdivisions
        specified by <ndiv>.
        TODO: improve this method
        """
        dmax = 0.0
        for i in range(ndiv+1):
            t = float(i) / ndiv
            p = self.point_at(t)
#            p.SVG_plot('#00cccc') #DEBUG
            d = arc.distance_to_point(p)
            if d > 0.0: # only if arc intersects segment(center,p)
                dmax = max(dmax, d)
        return dmax

    def __str__(self):
        """Concise string representation."""
        return "CubicBezier((%.4f, %.4f), (%.4f, %.4f), (%.4f, %.4f), (%.4f, %.4f))" % \
                (self.p1.x, self.p1.y, self.c1.x, self.c1.y,
                 self.c2.x, self.c2.y, self.p2.x, self.p2.y)

    def __repr__(self):
        """Precise string representation."""
        return "CubicBezier((%r, %r), (%r, %r), (%r, %r), (%r, %r))" % \
                (self.p1.x, self.p1.y, self.c1.x, self.c1.y,
                 self.c2.x, self.c2.y, self.p2.x, self.p2.y)

    def to_SVG_path(self):
        """Return a string with the SVG path 'd' attribute
        that corresponds with this curve.
        """
        return 'M %f %f C %f %f %f %f %f %f' % \
                (self.p1.x, self.p1.y, self.c1.x, self.c1.y,
                 self.c2.x, self.c2.y, self.p2.x, self.p2.y)

    def SVG_plot(self, color='#cccc99'):
        """Draw an SVG version of this curve for debugging/testing.
        Include control points, inflection points, and tangent lines.
        """
        if DEBUG_EFFECT is not None and DEBUG_LAYER is not None:
            style = 'fill:none;stroke:%s;stroke-width:0.25pt;stroke-opacity:1' % (color,)
            attrs = {'d': self.to_SVG_path(), 'style': style}
            DEBUG_EFFECT.create_path(attrs, layer=DEBUG_LAYER)
            # Draw control points and tangents
            self.c1.SVG_plot(color='#0000c0')
            self.c2.SVG_plot(color='#0000c0')
            tseg1 = Line(self.p1, self.c1)
            tseg1.SVG_plot()
            tseg2 = Line(self.p2, self.c2)
            tseg2.SVG_plot()
            # Draw inflection points if any
            t1, t2 = self.find_inflections()
            if t1 > 0.0:
                ip1 = self.point_at(t1)
                ip1.SVG_plot(color='#c00000')
            if t2 > 0.0:
                ip2 = self.point_at(t2)
                ip2.SVG_plot(color='#c00000')
            # Draw midpoint
            mp = self.point_at(0.5)
            mp.SVG_plot(color='#666666')

