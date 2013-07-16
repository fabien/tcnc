#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

"""
Created on Apr 20, 2013

@author: claude

====
"""
import geom
    
#==============================================================================
# Chan's Convex Hull O(n log h) - Tom Switzer <thomas.switzer@gmail.com>
# See http://tomswitzer.net/2010/12/2d-convex-hulls-chans-algorithm/
#==============================================================================
TURN_LEFT, TURN_RIGHT, TURN_NONE = (1, -1, 0)
 
def turn(p, q, r):
    """Returns -1, 0, 1 if p,q,r forms a right, straight, or left turn.
    """
    return cmp((q[0] - p[0])*(r[1] - p[1]) - (r[0] - p[0])*(q[1] - p[1]), 0)
 
def _keep_left(hull, r):
    while len(hull) > 1 and turn(hull[-2], hull[-1], r) != TURN_LEFT:
        hull.pop()
    if not len(hull) or hull[-1] != r:
        hull.append(r)
    return hull
 
def convex_hull_graham_scan(points):
    """Returns points on convex hull of an array of points in CCW order.
    
    Uses the Graham Scan algorithm.
    """
    points = sorted(points)
    lh = reduce(_keep_left, points, [])
    uh = reduce(_keep_left, reversed(points), [])
#    lh.extend(uh[i] for i in xrange(1, len(uh) - 1))
    lh.extend(uh[1:-1])
    return lh
 
def convex_hull_chan(pts):
    """Returns the points on the convex hull of pts in CCW order.
    
    Uses Chan's algorithm.
    """
#    for m in (1 << (1 << t) for t in xrange(len(pts))):
    for m in ((1 << t) for t in xrange(len(pts))):
        hulls = [convex_hull_graham_scan(pts[i:i + m])
                 for i in xrange(0, len(pts), m)]
        hull = [_min_hull_pt_pair(hulls)]
        for unused in xrange(m):
            p = _next_hull_pt_pair(hulls, hull[-1])
            if p == hull[0]:
                return [hulls[h][i] for h, i in hull]
            hull.append(p)
    return hull

def _rtangent(hull, p):
    """Return the index of the point in hull that the right tangent line from p
    to hull touches.
    """
    l, r = 0, len(hull)
    l_prev = turn(p, hull[0], hull[-1])
    l_next = turn(p, hull[0], hull[(l + 1) % r])
    while l < r:
        c = (l + r) / 2
        c_prev = turn(p, hull[c], hull[(c - 1) % len(hull)])
        c_next = turn(p, hull[c], hull[(c + 1) % len(hull)])
        c_side = turn(p, hull[l], hull[c])
        if c_prev != TURN_RIGHT and c_next != TURN_RIGHT:
            return c
        elif c_side == TURN_LEFT and (l_next == TURN_RIGHT or
                                      l_prev == l_next) or \
                c_side == TURN_RIGHT and c_prev == TURN_RIGHT:
            r = c               # Tangent touches left chain
        else:
            l = c + 1           # Tangent touches right chain
            l_prev = -c_next    # Switch sides
            l_next = turn(p, hull[l], hull[(l + 1) % len(hull)])
    return l
 
def _min_hull_pt_pair(hulls):
    """Returns the hull, point index pair that is minimal."""
    h, p = 0, 0
    for i in xrange(len(hulls)):
        j = min(xrange(len(hulls[i])), key=lambda j: hulls[i][j])
        if hulls[i][j] < hulls[h][p]:
            h, p = i, j
    return (h, p)
 
def _dist(p1, p2):
    """Euclidean distance squared between two points."""
    a = p1[0] - p2[0]
    b = p1[1] - p2[1]
    return (a * a) + (b * b)

def _next_hull_pt_pair(hulls, pair):
    """
    Returns the (hull, point) index pair of the next point in the convex
    hull.
    """
    p = hulls[pair[0]][pair[1]]
    nextpair = (pair[0], (pair[1] + 1) % len(hulls[pair[0]]))
    for h in (i for i in xrange(len(hulls)) if i != pair[0]):
        s = _rtangent(hulls[h], p)
        q, r = hulls[nextpair[0]][nextpair[1]], hulls[h][s]
        t = turn(p, q, r)
        if t == TURN_RIGHT or t == TURN_NONE and _dist(p,r) > _dist(p,q):
            nextpair = (h, s)
    return nextpair
 

#==============================================================================
# Area and centroid calculations for non self-intersecting closed polygons.
# See http://paulbourke.net/geometry/polygonmesh/
#==============================================================================

def area(points):
    """Return the area of a simple polygon.
    
    :param points: the polygon vertices.
    """
    area=0
    nmax = len(points) - 1
    for n in xrange(0 if points[0] == points[nmax] else -1, nmax):
        area += points[n].x*points[n+1].y - points[n+1].x*points[n].y
    return abs(area) / 2

def centroid(points):
    """Return the centroid of a simple polygon.
    
    :param points: the polygon vertices.
    """
    nmax = len(points) - 1
    x = 0
    y = 0
    area = 0
    for n in xrange(0 if points[0] == points[nmax] else -1, nmax):
        cross = points[n].x*points[n+1].y - points[n+1].x*points[n].y
        area += cross
        x += (points[n].x + points[n+1].x) * cross
        y += (points[n].y + points[n+1].y) * cross
    t = abs(area) * 6
    return geom.P(x/t, y/t)


#==============================================================================
#Portions of this code (point in polygon test) are derived from:
#http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
#and uses the following license:
#Copyright (c) 1970-2003, Wm. Randolph Franklin
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal in
#the Software without restriction, including without limitation the rights to
#use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
#of the Software, and to permit persons to whom the Software is furnished to do
#so, subject to the following conditions:
#Redistributions of source code must retain the above copyright notice,
#this list of conditions and the following disclaimers.
#Redistributions in binary form must reproduce the above copyright notice in the
#documentation and/or other materials provided with the distribution.
#The name of W. Randolph Franklin may not be used to endorse or promote products
#derived from this Software without specific prior written permission.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#==============================================================================

def point_inside(points, p):
    """Return True if point `p` is inside the polygon defined by `points`.
    
    See http://paulbourke.net/geometry/polygonmesh/
    
    Also: http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
    
    :param points:
    :param p:
    """
    is_inside = False
    j = -1
    for i in range(len(points)):
        p1 = points[i]
        p2 = points[j]
        # This is a tricky conditional - see W. R. Franklin's comments
        if ((p1.y > p.y) != (p2.y > p.y)) and \
           p.x < (((p2.x - p1.x) * (p.y - p1.y) / (p2.y - p1.y)) + p1.x):
            is_inside = not is_inside
        j = i
    return is_inside
        
        
