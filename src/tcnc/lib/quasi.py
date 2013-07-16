"""
Python port of quasi.c which was originally written by Eric Weeks
weeks@physics.emory.edu

See: http://www.physics.emory.edu/~weeks/software/quasic.html

Mostly unchanged except to make it a little more pythonic and:
    - Removed Postscript output and main()
    - Fixed divide by zero exception for even symmetries
    - Added segment connection to vertices options
    - Added S.J. King's coloring method
    
====
"""

import math
import logging
logger = logging.getLogger(__name__)

class QuasiPlotter(object):
    """Quasi plotter base class.
    Subclass this to produce output.
    
    Does nothing by default.
    """

    def plot_polygon(self, vertices):
        """Draw a polygon.
        
        :param vertices: a list of tuples containing the (x,y) coordinates of the
            polygon vertices.
        :return: True if the polygon is clipped and segments don't need to
            be drawn. False is default.
        """
        return False

    def plot_segment(self, x1, y1, x2, y2):
        """Draw a line segment.
        
        :param x1: x coordinate of start point
        :param y1: y coordinate of start point
        :param x2: x coordinate of end point
        :param y2: y coordinate of end point
        """
        pass
            
    def set_fill_color(self, color):
        """Set the current polygon fill color.
        
        :param color: a grayscale value between 0.0 and 1.0
        :type color: float
        """
        pass
        

class Quasi(object):
    # Segment connection types
    SEG_NONE = 0
    """No segment connection."""
    SEG_MIDP_ACUTE = 1
    """Connect midpoints of polygon edges that meet at an acute angle."""
    SEG_MIDP_OBTUSE = 2
    """Connect midpoints of polygon edges that meet at an obtuse angle."""
    SEG_MIDP_CROSS = 3
    """Connect midpoints of polygon edges to form a cross."""
    SEG_MIDP_RECT = 4
    """Connect midpoints of polygon edges to form a rectangle."""
    SEG_VERT_ACUTE = 5
    """Connect polygon vertices whose edges form an acute angle."""
    SEG_VERT_OBTUSE = 6
    """Connect polygon vertices whose edges form an obtuse angle."""
    SEG_VERT_CROSS = 7
    """Connect polygon vertices to form a cross."""
    
    def _segtype_midp(self, segtype):
        return 0 < segtype < 5

    # Eris Week's original values.
    #offset_salt_y = 0.1132
    #offset_salt_x = 0.2137
    offset_salt_y = 0.1618
    """Random-ish offset to avoid more than one line intersecting."""
    offset_salt_x = 0.314159
    """Random-ish offset to avoid more than one line intersecting."""
    segment_ratio = 0.5
    """Ratio that determines edge midpoint."""
    skinnyfat_ratio = 0.5
    """Ration that determines whether a rhombus is fat or skinny"""
    symmetry = 5
    """Degrees of symmetry."""
    numlines = 30
    """Number of lines. A larger number enables more tiles to be generated."""
    color_fill = False
    """Color fill polygons"""
    color_by_polytype = False
    """Fill color by rhombus type (fat/skinny)."""
    plotter = None
    """Plotter to draw output."""
    
    _rad_incr = 3 #0.4
    
    # note: implemented as a dictionary independent of segtype order
    _SWAPPED_SEGTYPES = {SEG_NONE: SEG_NONE,
                         SEG_MIDP_ACUTE: SEG_MIDP_OBTUSE,
                         SEG_MIDP_OBTUSE: SEG_MIDP_ACUTE,
                         SEG_MIDP_CROSS: SEG_MIDP_CROSS,
                         SEG_MIDP_RECT: SEG_MIDP_RECT,
                         SEG_VERT_ACUTE: SEG_VERT_OBTUSE,
                         SEG_VERT_OBTUSE: SEG_VERT_ACUTE,
                         SEG_VERT_CROSS: SEG_VERT_CROSS}
    
    def __init__(self, symmetry=5, segtype_skinny=SEG_NONE,
                 segtype_fat=SEG_NONE, plotter=None):
        """
        :param symmetry: degrees of symmetry. Must be at least two.
        :param segtype_skinny: segment connection type for skinny rhombuses.
        :param segtype_fat: segment connection type for fat rhombuses.
        :param plotter: plotter to draw output.
        :type plotter: QuasiPlotter
        """
        self.plotter = plotter
        if plotter is None:
            self.plotter = QuasiPlotter()
        self.symmetry = symmetry
        self.segtype_skinny = segtype_skinny
        self.segtype_fat = segtype_fat
        self._current_fill_color = 0.2

    def quasi(self):
        """Draw tiling.
        """    
        # t is 1st direction, r is 2nd.  look for intersection between pairs
        # of lines in these two directions. (will be x0,y0) 
        index = [0,] * self.symmetry
        maxline = self.numlines - 1
        minline = maxline / 2
        minrad = 0.0
        vx, vy, mm, vb = self._init_vectors()
        
        while minrad <= float(maxline):
            rad1 = minrad * minrad
            minrad += self._rad_incr
            rad2 = minrad * minrad
            for n in range(1, maxline):
                n2 = (n - minline) * (n - minline)
                for m in range(1, maxline):
                    rad = float(n2 + (m - minline) * (m - minline))
                    if rad1 <= rad < rad2:
                        for t in range(self.symmetry - 1):
                            for r in range(t + 1, self.symmetry):
                                x0 = (vb[t][n] - vb[r][m]) / (mm[r] - mm[t])
                                y0 = mm[t]*x0 + vb[t][n]
                                do_plot = True
                                for i in range(self.symmetry):
                                    if i != t and i != r:
                                        dx = -x0 * vy[i] + (y0 - vb[i][0]) * vx[i]
                                        index[i] = int(-dx)
                                        if index[i] > self.numlines - 3 or index[i] < 1:
                                            do_plot = False
                                            break
                                if do_plot:
                                    index[t] = n - 1
                                    index[r] = m - 1
                                    self._plot(vx, vy, index, t, r)

    def _init_vectors(self):
        """Initialize vectors."""
        vx = []
        vy = []
        mm = []
        vb = []
        phi = 2*math.pi / self.symmetry
        # This avoids a div by 0 exception for even symmetries
        if self.symmetry % 2 == 0:
            phi += 0.000001
        for t in range(self.symmetry):
            angle = phi * t
            x = math.cos(angle)
            y = math.sin(angle)
            m = y / x
            vx.append(x)
            vy.append(y)
            mm.append(m)
            vr = []
            for r in range(self.numlines):
                y1 = (y * (t * self.offset_salt_y)) - (x * (r - self.numlines/2))
                x1 = (x * (t * self.offset_salt_x)) + (y * (r - self.numlines/2))
                intercept = y1 - m*x1
                vr.append(intercept)
            vb.append(vr)
        return (vx, vy, mm, vb)


    def _plot(self, vx, vy, index, t, r):
        x0 = 0.0
        y0 = 0.0
        for i in range(self.symmetry):
            x0 += vx[i] * index[i]
            y0 += vy[i] * index[i]
        vertices = (
            (x0, y0),
            (x0 + vx[t], y0 + vy[t]),
            (x0 + vx[t] + vx[r], y0 + vy[t] + vy[r]),
            (x0 + vx[r], y0 + vy[r]))
        if self.color_fill:
            self.plotter.set_fill_color(self._get_fill_color(vx, vy,
                                                             t, r, index[i]))
        polygon_clipped = self.plotter.plot_polygon(vertices)
        if (self.segtype_skinny > 0 or self.segtype_fat > 0) \
        and not polygon_clipped:
            # Determine if the polygon is fat or skinny
            d1 = self._distance(vertices[0], vertices[2])
            d2 = self._distance(vertices[1], vertices[3])
            is_skinny = (min(d1, d2) / max(d1, d2)) < self.skinnyfat_ratio
            segtype = self.segtype_skinny if is_skinny else self.segtype_fat
            if d1 > d2:
                segtype = Quasi._SWAPPED_SEGTYPES[segtype]
            if self._segtype_midp(segtype):
                # Calculate segment endpoints
                midx1 = x0 + vx[t] * self.segment_ratio
                midy1 = y0 + vy[t] * self.segment_ratio
                midx2 = x0 + vx[t] + vx[r] * self.segment_ratio
                midy2 = y0 + vy[t] + vy[r] * self.segment_ratio
                midx3 = x0 + vx[r] + vx[t] * self.segment_ratio
                midy3 = y0 + vy[r] + vy[t] * self.segment_ratio
                midx4 = x0 + vx[r] * self.segment_ratio
                midy4 = y0 + vy[r] * self.segment_ratio
                if segtype == Quasi.SEG_MIDP_ACUTE:
                    self.plotter.plot_segment(midx1, midy1, midx2, midy2)
                    self.plotter.plot_segment(midx3, midy3, midx4, midy4)
                elif segtype == Quasi.SEG_MIDP_OBTUSE:
                    self.plotter.plot_segment(midx1, midy1, midx4, midy4)
                    self.plotter.plot_segment(midx2, midy2, midx3, midy3)
                elif segtype == Quasi.SEG_MIDP_CROSS:
                    self.plotter.plot_segment(midx1, midy1, midx3, midy3)
                    self.plotter.plot_segment(midx2, midy2, midx4, midy4)
                elif segtype == Quasi.SEG_MIDP_RECT:
                    self.plotter.plot_segment(midx1, midy1, midx2, midy2)
                    self.plotter.plot_segment(midx3, midy3, midx4, midy4)
                    self.plotter.plot_segment(midx1, midy1, midx4, midy4)
                    self.plotter.plot_segment(midx2, midy2, midx3, midy3)
            else:
                if segtype == Quasi.SEG_VERT_ACUTE:
                    self.plotter.plot_segment(vertices[1][0], vertices[1][1],
                                              vertices[3][0], vertices[3][1])
                elif segtype == Quasi.SEG_VERT_OBTUSE:
                    self.plotter.plot_segment(vertices[0][0], vertices[0][1],
                                              vertices[2][0], vertices[2][1])
                elif segtype == Quasi.SEG_VERT_CROSS:
                    self.plotter.plot_segment(vertices[0][0], vertices[0][1],
                                              vertices[2][0], vertices[2][1])
                    self.plotter.plot_segment(vertices[1][0], vertices[1][1],
                                              vertices[3][0], vertices[3][1])

    def _distance(self, p1, p2):
        """Distance between two points."""
        return math.hypot(p2[0] - p1[0], p2[1] - p1[1])  
     
    def _get_fill_color(self, vx, vy, t, r, i):
        self._current_fill_color += .05
        if self._current_fill_color > 1.0:
            self._current_fill_color = 0.2
        if self.color_by_polytype:
            self._current_fill_color = 0.0
            for i in range(self.symmetry):
                self._current_fill_color += i
            while self._current_fill_color > (self.symmetry - 1) / 2.0:
                self._current_fill_color -= (self.symmetry - 1) / 2.0
            self._current_fill_color = self._current_fill_color / ((self.symmetry - 1) / 2.0) * 0.8 + 0.1
            self._current_fill_color += math.fabs(vx[t] * vx[r] + vy[t] * vy[r]) # dot product
            if self._current_fill_color > 1.0:
                self._current_fill_color -= 1.0
        return self._current_fill_color
    
    @staticmethod
    def get_fill_color_rgb(fill_color):
        """Return an RGB fill color given a grayscale color value.
        
        Based on quasi_colour.c by S.J. King.
        
        :param fill_color: a grayscale color value in the range 0.0 to 1.0
        :return: a 3-tuple (R, G, B)
        """
        phi = fill_color * 2.0 * math.pi
        theta = (math.pi / 2.0) * math.sin(3.0 * phi) + math.pi / 2.0
        r = 0.5 * math.sin(theta) * math.cos(phi) + 0.5
        g = 0.5 * math.sin(theta) * math.sin(phi) + 0.5
        b = 0.5 * math.cos(theta) + 0.5
        return (r, g, b)
