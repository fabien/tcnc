'''
Python port of quasi.c which was originally written by Eric Weeks
link: http://www.physics.emory.edu/~weeks/software/quasi.html
email: weeks@physics.emory.edu

Mostly unchanged except to make it a little more pythonic...

- Added segment connection to vertices options
'''

import math
from optparse import OptionParser
import logging
logger = logging.getLogger(__name__)

from lib import geom

class QuasiPlotter(object):
    ''''''
    def __init__(self, transform=geom.IdentityTransform2D, clip_region=None):
        self.polygons = []
        self.segments = []
        self.transform = transform
        self.clip_region = clip_region
        
    def plot_polygon(self, vertices):
        '''Draw a polygon. Return True if the polygon is clipped.'''
        clipped = self._is_clipped(vertices)
        if not clipped:
            vertices = map(lambda p: geom.matrix_apply_to_point(self.transform, p), vertices)
            self.polygons.append(self._transform(vertices))
        return clipped

    def plot_segment(self, x1, y1, x2, y2):
        '''Draw a segment.'''
        endpoints = geom.Line(geom.P(x1, y1), geom.P(x2, y2))
        if not self._is_clipped(endpoints):
            self.segments.append(endpoints.transform(self.transform))
        
    def _is_clipped(self, points):
        return self.clip_region is not None and \
            not self.clip_region.all_points_inside(points)
        

class Quasi(object):
    # Segment connection types
    SEG_NONE = 0
    SEG_MIDP_ACUTE = 1
    SEG_MIDP_OBTUSE = 2
    SEG_MIDP_CROSS = 3
    SEG_MIDP_RECT = 4
    SEG_VERT_ACUTE = 5
    SEG_VERT_OBTUSE = 6
    SEG_VERT_CROSS = 7
    
    def _segtype_midp(self, segtype):
        return 0 < segtype < 5

    # Eris Week's original values.
    #offset_salt_y = 0.1132
    #offset_salt_x = 0.2137
    offset_salt_y = 0.1618
    offset_salt_x = 0.314159
    do_polygons = True
    do_segments = True
    segment_ratio = 0.2
    skinnyfat_ratio = 0.2
    symmetry = 5
    numlines = 30
    color_by_polytype = False
    plotter = QuasiPlotter()
    
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
        """"""
        if plotter is not None:
            self.plotter = plotter
        self.symmetry = symmetry
        self.segtype_skinny = segtype_skinny
        self.segtype_fat = segtype_fat
        self._current_fill_color = 0.2

    def quasi(self):
        """"""    
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
                                    self._plot_segments(vx, vy, index, t, r)

    def _init_vectors(self):
        """Initialize vectors."""
        vx = []
        vy = []
        mm = []
        vb = []
        T0 = 2.0
        if self.symmetry % 2 == 0:
            # This avoids a div by 0 exception for even symmetries
            T0 = 2.00001
        for t in range(self.symmetry):
            phi = ((t * T0) / self.symmetry) * math.pi
            x = math.cos(phi)
            y = math.sin(phi)
            m = y / x
            vx.append(x)
            vy.append(y)
            mm.append(m)
#            logger.debug('x,y,m: %.5f, %.5f, %.5f' % (x, y, m))
            vr = []
            for r in range(self.numlines):
                y1 = y * (t * self.offset_salt_y) - x * (r - self.numlines/2) 
                x1 = x * (t * self.offset_salt_x) + y * (r - self.numlines/2)
                intercept = y1 - m*x1
                vr.append(intercept)
            vb.append(vr)
        return (vx, vy, mm, vb)

    def _plot_segments(self, vx, vy, index, t, r):
        x0 = 0.0
        y0 = 0.0
        for i in range(self.symmetry):
            x0 += vx[i] * index[i]
            y0 += vy[i] * index[i]
        vertices = (
            geom.P(x0, y0),
            geom.P(x0 + vx[t], y0 + vy[t]),
            geom.P(x0 + vx[t] + vx[r], y0 + vy[t] + vy[r]),
            geom.P(x0 + vx[r], y0 + vy[r]))
        polygon_clipped = False
        if self.do_polygons:
            if self.plotter.fill_polygons:
                self.plotter.fill_color = self.get_fill_color(vx, vy,
                                                              t, r, index[i])
            polygon_clipped = self.plotter.plot_polygon(vertices)
        if self.do_segments \
        and (self.segtype_skinny > 0 or self.segtype_fat > 0) \
        and not polygon_clipped:
            # Determine if the polygon is fat or skinny
            d1 = vertices[0].distance2(vertices[2])
            d2 = vertices[1].distance2(vertices[3])
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
    
    def get_fill_color(self, vx, vy, t, r, i):
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
        # Based on quasi_colour.c by S.J. King zh84@lycos.com
        phi = fill_color * 2.0 * math.pi
        theta = (math.pi / 2.0) * math.sin(3.0 * phi) + math.pi / 2.0
        r = 0.5 * math.sin(theta) * math.cos(phi) + 0.5
        g = 0.5 * math.sin(theta) * math.sin(phi) + 0.5
        b = 0.5 * math.cos(theta) + 0.5
        return (r, g, b)


class PSQuasiPlotter(QuasiPlotter):
    '''Postscript plotter'''
    clipped = False
    scale = 1.0
    offsetx = 0.0
    offsety = 0.0
    xcenter = 0.0
    ycenter = 0.0
    rotate = False
    window = 1.0
    stroke_width = 0.015
    
    def __init__(self, *args, **kwargs):
        super(PSQuasiPlotter, self).__init__(*args, **kwargs)        
    
    def plot_polygon(self, vertices):
        ''''''
        self._psplot(vertices[0][0], vertices[0][1], 0)
        for p in vertices[1:-1]:
            self._psplot(p[0], p[1], 1)
        self._psplot(vertices[-1][0], vertices[-1][1], 2)
        return False

    def plot_segment(self, x1, y1, x2, y2):
        ''''''
        self._psplot(x1, y1, 0)
        self._psplot(x2, y2, 2)

    def ps_header(self):
        '''Output header'''
        # PostScript Header (taken from CGLE output)
        print("%!PS-Adobe-1.0 ")
        print("%%BoundingBox: -1 -1 766.354 567.929 ")
        print("%%EndComments ")
        print("%%EndProlog ")
        print("gsave ")
        print(" ")
        print("/f {findfont exch scalefont setfont} bind def ")
        print("/s {show} bind def ")
        print("/ps {true charpath} bind def ")
        print("/l {lineto} bind def ")
        print("/m {newpath moveto} bind def ")
        print("/sg {setgray} bind def")
        print("/a {stroke} bind def")
        print("/cp {closepath} bind def")
        print("/g {gsave} bind def")
        print("/h {grestore} bind def")
        print("matrix currentmatrix /originmat exch def ")
        print("/umatrix {originmat matrix concatmatrix setmatrix} def ")
        print(" ")
        print("% Flipping coord system ")
        print("[8.35928e-09 28.3465 -28.3465 8.35928e-09 609.449 28.6299] umatrix ")
        print("[] 0 setdash ")
        print("0 0 0 setrgbcolor ")
        print("0 0 m ")
        print("%.6f setlinewidth " % self.stroke_width)
    
    def ps_footer(self):
        '''Output footer'''
        print("showpage grestore ")
        print("%%Trailer")
    
    def getdx(self, x, center):
        dx = (x - center) / self.window
        dx = 0.5 * (dx + 1.0)
        return dx

    def _psplot(self, x, y, plotflag):
        swap = 0.0
        cmx = 0.0
        cmy = 0.0        # x,y in centimeters 
        # plotflag:  0 = start line; 1 = lineto; 2 = endpoint 
    
        dx = self.getdx(x, self.xcenter)
        dy = self.getdx(y, self.ycenter)
        if self.rotate:
            swap = dx
            dx = dy
            dy = swap
    
        if ((dx < 1.3) and (dy < 1.0) and (dx > 0) and (dy > 0)):        # in window 
            cmx = (dx * self.scale) + self.offsetx
            cmy = (dy * self.scale) + self.offsety
            if plotflag < 1:
                print("%.2f %.2f m" % (cmx, cmy))
            else:
                if self.clipped:
                    print("%.2f %.2f m" % (cmx, cmy))
                print("%.2f %.2f l" % (cmx, cmy))
                if plotflag == 2:
                    print("cp")
                    if self.fillon:
                        print("g")
                        if self.use_color:
                            print("%.2f %.2f %.2f sc" % self.get_fill_color_rgb())
                        else:
                            print("%.1f sg" % self.get_fill_color())
                        print("fill")
                        print("h")
                    else:
                        print("%.1f sg" % self.get_stroke_color())
                    print("a")
            self.clipped = False
        else:
            self.clipped = True


def main():
    parser = OptionParser()

    parser.add_option('-s', type='float', default='15.0', dest='scale', help='size of output (width) in cm [15]')
    parser.add_option('-f', type='int', default='0', dest='fillon', help='fill polygons')
    parser.add_option('-z', type='int', default='0', dest='zfill', help='fill color is according to polygon type')
    parser.add_option('-F', type='int', default='0', dest='rotate', help='flip picture 90 degrees')
    parser.add_option('-m', type='float', default='1.0', dest='magnifya', help='magnification factor')
    parser.add_option('-M', type='float', default='0.0', dest='midon0', help='midpoint type for skinny diamonds')
    parser.add_option('-N', type='float', default='0.0', dest='midon1', help='midpoint type for fat diamonds')
    parser.add_option('-S', type='int', default='5', dest='symmetry', help='degrees of symmetry [5]')
    parser.add_option('-n', type='int', default='30', dest='maxmax', help='number of lines to use [30]')
    parser.add_option('-w', type='float', default='20.0', dest='window', help='clipping window')

    (options, args) = parser.parse_args()

    fillon = bool(options.fillon)
    zfill = bool(options.zfill)
    midon = [int(options.midon0), int(options.midon1)]
    symmetry = int(options.symmetry)
    maxmax = int(options.maxmax)

    plotter = PSQuasiPlotter()
    plotter.offsetx = 2.0                    # lower left corner of picture */
    plotter.offsety = 2.0
    plotter.scale = float(options.scale)
    plotter.rotate = bool(options.rotate)
    plotter.window = float(options.window) / float(options.magnifya)
    plotter.fillon = fillon
    plotter.ps_header()
    quasi = Quasi()
    quasi.quasi(zfill, midon, symmetry, maxmax, plotter)
    plotter.ps_footer()
    
    exit(0)


# Uncomment main() to use as command line tool compatible with the original
# quasi.c which produces postscript to stdout.
#main()
