"""Voronoi diagram calculator/ Delaunay triangulator

Based on Steve Fortune's original code:
    http://netlib.bell-labs.com/cm/cs/who/sjf/index.html
    
Translated to Python by Bill Simons September, 2005

Simplified a bit by Claude Zervas, 2013

Calculate Delaunay triangulation or the Voronoi polygons for a set of 
2D input points.

Derived from code bearing the following notice::

    The author of this software is Steven Fortune. Copyright (c) 1994 by AT&T
    Bell Laboratories.
    
    Permission to use, copy, modify, and distribute this software for any
    purpose without fee is hereby granted, provided that this entire notice
    is included in all copies of any software which is or includes a copy
    or modification of this software and in all copies of the supporting
    documentation for such software.
    THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
    WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
    REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
    OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
    
Comments were incorporated from Shane O'Sullivan's translation of the 
original code into C++: http://mapviewer.skynet.ie/voronoi.html

====
"""
import math
import sys
import logging
logger = logging.getLogger(__name__)

# Tolerance for floating point comparisons
_EPSILON = 1e-9

def voronoi_diagram(points, triangulate=False):
    """Takes a list of point objects (2-tuples with x and y coordinates).
    
    Returns a 4-tuple of:

    1.  A list of 2-tuples, which are the x,y coordinates of the 
        Voronoi diagram vertices
    2.  A list of 3-tuples (a,b,c) which are the equations of the
        lines in the Voronoi diagram: a*x + b*y = c
    3.  A list of 3-tuples, (i, v1, v2) representing edges of the 
        Voronoi diagram.  i is the index of the line, v1 and v2 are
        the indices of the vertices at the end of the edge. If 
        v1 or v2 is -1, the line extends to infinity.
    4.  An optional list of edges describing the Delauney Triangulation
        of the voronoi diagram. Each edge is a 4-tuple containing the
        endpoints of the edge line segment.
    """
    diagram = voronoi(points, triangulate)
    return (diagram.vertices, diagram.lines, diagram.edges,
            diagram.delauney_edges)


class DiagramContext(object):
    """The Voronoi diagram and optional Delauney triangulation.
    """
    def __init__(self):
        """"""
        self.vertices = []
        """list of vertex 2-tuples: (x,y)"""
        self.lines = []
        """List of lines. An edge is a 3-tuple (a, b, c) containing
        line equation of the form a*x+b*y = c"""
        self.edges = []
        """List of edges. An edge is a 3-tuple:
        (line index, vertex 1 index, vertex 2 index)
        if either vertex index is -1, the edge extends to infinity"""
        self.triangles = []
        """List of 3-tuples of vertex indices"""
        self.delauney_edges = []
        """List of 4-tuples of Delauney edges"""

    def add_vertex(self, s):
        self.vertices.append((s.x, s.y))
        logger.debug("vertex(%d) at %f %f" % (s.sitenum, s.x, s.y))

    def add_triple(self, s1, s2, s3):
        self.triangles.append((s1.sitenum, s2.sitenum, s3.sitenum))
        logger.debug("circle through left=%d right=%d bottom=%d" % \
                     (s1.sitenum, s2.sitenum, s3.sitenum))

    def add_bisector(self, edge, triangulate):
        self.lines.append((edge.a, edge.b, edge.c))
        if triangulate:
            self.delauney_edges.append((edge.reg[0].x, edge.reg[0].y,
                                        edge.reg[1].x, edge.reg[1].y))
        logger.debug("line(%d) %gx+%gy=%g, bisecting %d %d" % \
                     (edge.edgenum, edge.a, edge.b, edge.c,
                      edge.reg[0].sitenum, edge.reg[1].sitenum))

    def add_edge(self, edge):
        sitenumL = -1
        if edge.ep[_Edge.LE] is not None:
            sitenumL = edge.ep[_Edge.LE].sitenum
        sitenumR = -1
        if edge.ep[_Edge.RE] is not None:
            sitenumR = edge.ep[_Edge.RE].sitenum
        self.edges.append((edge.edgenum, sitenumL, sitenumR))


def voronoi(points, triangulate=False):
    """Create a Voronoi diagram.
    
    :param points: a list of point 2-tuples (x, y)
    :param triangulate: True if Delauney triangulation is to be performed.
        Default is False.
        
    :return: a :class:`DiagramContext` containing the Voronoi diagram output.
    """
    diagram = DiagramContext()
    site_list = _SiteList(points)
    nsites = len(site_list)
    if nsites < 2:
        return
    edgeList  = _EdgeList(site_list.xmin, site_list.xmax, nsites)
    priorityQ = _PriorityQueue(site_list.ymin, site_list.ymax, nsites)
    siteIter = iter(site_list)
    bottomsite = siteIter.next()
    newsite = siteIter.next()
    minpt = _Site(sys.float_info.min, sys.float_info.min)
    site_count = 0
    while True:
        if not priorityQ.isEmpty():
            minpt = priorityQ.getMinPt()

        if newsite and (priorityQ.isEmpty() or cmp(newsite, minpt) < 0):
            # get first Halfedge to the LEFT and RIGHT of the new site 
            lbnd = edgeList.pop_leftbnd(newsite) 
            rbnd = lbnd.right                    
            
            # if this halfedge has no edge, bot = bottom site (whatever that is)
            # create a new edge that bisects
            bot  = lbnd.rightreg(bottomsite)     
            edge = _Edge.bisect(bot, newsite)      
            diagram.add_bisector(edge, triangulate)
            
            # create a new Halfedge, setting its pm field to 0 and insert 
            # this new bisector edge between the left and right vectors in
            # a linked list
            bisector = _Halfedge(edge, _Edge.LE)    
            edgeList.insert(lbnd, bisector)       

            # if the new bisector intersects with the left edge, remove 
            # the left edge's vertex, and put in the new one
            p = lbnd.intersect(bisector)
            if p is not None:
                priorityQ.delete(lbnd)
                priorityQ.insert(lbnd, p, newsite.distance(p))

            # create a new _Halfedge, setting its pm field to 1
            # insert the new _Halfedge to the right of the original bisector
            lbnd = bisector
            bisector = _Halfedge(edge, _Edge.RE)     
            edgeList.insert(lbnd, bisector)        

            # if this new bisector intersects with the right Halfedge
            p = bisector.intersect(rbnd)
            if p is not None:
                # push the Halfedge into the ordered linked list of vertices
                priorityQ.insert(bisector, p, newsite.distance(p))
            
            try:
                newsite = siteIter.next()
            except StopIteration:
                newsite = None

        elif not priorityQ.isEmpty():
            # intersection is smallest - this is a vector (circle) event 

            # pop the Halfedge with the lowest vector off the ordered list of 
            # vectors.  Get the Halfedge to the left and right of the above HE
            # and also the Halfedge to the right of the right HE
            lbnd  = priorityQ.popMinHalfedge()      
            llbnd = lbnd.left               
            rbnd  = lbnd.right              
            rrbnd = rbnd.right              
            
            # get the _Site to the left of the left HE and to the right of
            # the right HE which it bisects
            bot = lbnd.leftreg(bottomsite)  
            top = rbnd.rightreg(bottomsite) 
            
            # output the triple of sites, stating that a circle goes through them
            mid = lbnd.rightreg(bottomsite)
            diagram.add_triple(bot, top, mid)          

            # get the vertex that caused this event and set the vertex number
            # couldn't do this earlier since we didn't know when it would be processed
            v = lbnd.vertex
            v.sitenum = site_count
            site_count += 1
            diagram.add_vertex(v)
            
            # set the endpoint of the left and right Halfedge to be this vector
            if lbnd.edge.setEndpoint(lbnd.pm, v):
                diagram.add_edge(lbnd.edge)
            
            if rbnd.edge.setEndpoint(rbnd.pm, v):
                diagram.add_edge(rbnd.edge)

            
            # delete the lowest HE, remove all vertex events to do with the 
            # right HE and delete the right HE
            edgeList.delete(lbnd)           
            priorityQ.delete(rbnd)
            edgeList.delete(rbnd)
            
            
            # if the site to the left of the event is higher than the Site
            # to the right of it, then swap them and set 'pm' to RIGHT
            pm = _Edge.LE
            if bot.y > top.y:
                bot,top = top,bot
                pm = _Edge.RE

            # Create an _Edge (or line) that is between the two Sites.  This 
            # creates the formula of the line, and assigns a line number to it
            edge = _Edge.bisect(bot, top)     
            diagram.add_bisector(edge, triangulate)

            # create a HE from the edge 
            bisector = _Halfedge(edge, pm)    
            
            # insert the new bisector to the right of the left HE
            # set one endpoint to the new edge to be the vector point 'v'
            # If the site to the left of this bisector is higher than the right
            # Site, then this endpoint is put in position 0; otherwise in pos 1
            edgeList.insert(llbnd, bisector) 
            if edge.setEndpoint(_Edge.RE - pm, v):
                diagram.add_edge(edge)
            
            # if left HE and the new bisector don't intersect, then delete 
            # the left HE, and reinsert it 
            p = llbnd.intersect(bisector)
            if p is not None:
                priorityQ.delete(llbnd);
                priorityQ.insert(llbnd, p, bot.distance(p))

            # if right HE and the new bisector don't intersect, then reinsert it 
            p = bisector.intersect(rrbnd)
            if p is not None:
                priorityQ.insert(bisector, p, bot.distance(p))
        else:
            break

    he = edgeList.leftend.right
    while he is not edgeList.rightend:
        diagram.add_edge(he.edge)
        he = he.right
    
    return diagram

#------------------------------------------------------------------
class _Site(object):
    def __init__(self, x=0.0, y=0.0, sitenum=0):
        self.x = x
        self.y = y
        self.sitenum = sitenum

    def __str__(self):
        return "_Site(x=%g, y=%g, sitenum=%d)" % (self.x, self.y, self.sitenum)

    def __cmp__(self, other):
        if self.y < other.y:
            return -1
        elif self.y > other.y:
            return 1
        elif self.x < other.x:
            return -1
        elif self.x > other.x:
            return 1
        else:
            return 0

    def distance(self,other):
        dx = self.x - other.x
        dy = self.y - other.y
        return math.sqrt(dx*dx + dy*dy)

#------------------------------------------------------------------
class _Edge(object):
    LE = 0
    RE = 1
    DELETED = {}   # marker value
    _edgeNum = 0
 
    def __init__(self):
        self.a = 0.0
        self.b = 0.0
        self.c = 0.0
        self.ep  = [None, None]
        self.reg = [None, None]
        self.edgenum = 0

    def __str__(self):
        return "_Edge(edgenum=%d, a=%g, b=%g, c=%g, ep=%s, reg=%s)" % \
                (self.edgenum, self.a, self.b, self.c,
                 str(self.ep), str(self.reg))

    def setEndpoint(self, lrFlag, site):
        self.ep[lrFlag] = site
        if self.ep[_Edge.RE - lrFlag] is None:
            return False
        return True

    @staticmethod
    def bisect(s1, s2):
        newedge = _Edge()
        newedge.reg[0] = s1 # store the sites that this edge is bisecting
        newedge.reg[1] = s2

        # to begin with, there are no endpoints on the bisector
        # - it goes to infinity
        # ep[0] and ep[1] are None

        # get the difference in x dist between the sites
        dx = s2.x - s1.x
        dy = s2.y - s1.y
        adx = abs(dx)  # make sure that the difference in positive
        ady = abs(dy)
        
        # get the slope of the line
        newedge.c = s1.x*dx + s1.y*dy + (dx*dx + dy*dy)*0.5
        if adx > ady :
            # set formula of line, with x fixed to 1
            newedge.a = 1.0
            newedge.b = dy / dx
            newedge.c /= dx
        else:
            # set formula of line, with y fixed to 1
            newedge.b = 1.0
            newedge.a = dx / dy
            newedge.c /= dy

        newedge.edgenum = _Edge._edgeNum
        _Edge._edgeNum += 1
        return newedge


#------------------------------------------------------------------
class _Halfedge(object):
    def __init__(self, edge=None, pm=_Edge.LE):
        self.left = None   # left _Halfedge in the edge list
        self.right = None   # right _Halfedge in the edge list
        self.qnext = None   # priority queue linked list pointer
        self.edge = edge   # edge list _Edge
        self.pm = pm
        self.vertex = None  # _Site()
        self.ystar = sys.float_info.max

    def __str__(self):
        return "_Halfedge(left=%s, right=%s, edge=%s, pm=%s, vertex=%s, ystar=%s)" % \
            (str(self.left), str(self.right), str(self.edge), str(self.pm),
             str(self.vertex), str(self.ystar))

    def __cmp__(self, other):
        if self.ystar > other.ystar:
            return 1
        elif self.ystar < other.ystar:
            return -1
        elif self.vertex.x > other.vertex.x:
            return 1
        elif self.vertex.x < other.vertex.x:
            return -1
        else:
            return 0

    def leftreg(self, default):
        if not self.edge: 
            return default
        elif self.pm == _Edge.LE:
            return self.edge.reg[_Edge.LE]
        else:
            return self.edge.reg[_Edge.RE]

    def rightreg(self, default):
        if not self.edge: 
            return default
        elif self.pm == _Edge.LE:
            return self.edge.reg[_Edge.RE]
        else:
            return self.edge.reg[_Edge.LE]

    # returns True if p is to right of halfedge self
    def isPointRightOf(self, pt):
        e = self.edge
        topsite = e.reg[1]
        right_of_site = pt.x > topsite.x
        
        if right_of_site and self.pm == _Edge.LE: 
            return True
        
        if not right_of_site and self.pm == _Edge.RE:
            return False
        
        if _float_eq(e.a, 1.0):
            dyp = pt.y - topsite.y
            dxp = pt.x - topsite.x
            fast = False
            if (not right_of_site and e.b < 0.0) or (right_of_site and e.b >= 0.0):
                above = dyp >= e.b * dxp
                fast = above
            else:
                above = pt.x + pt.y * e.b > e.c
                if(e.b < 0.0):
                    above = not above
                if not above:
                    fast = True
            if not fast:
                dxs = topsite.x - (e.reg[0]).x
                above = e.b * (dxp*dxp - dyp*dyp) < \
                        dxs*dyp*(1.0 + 2.0*dxp/dxs + e.b*e.b)
                if e.b < 0.0:
                    above = not above
        else:  # e.b == 1.0 
            yl = e.c - e.a * pt.x
            t1 = pt.y - yl
            t2 = pt.x - topsite.x
            t3 = yl - topsite.y
            above = t1*t1 > (t2*t2 + t3*t3)
        
        if self.pm == _Edge.LE:
            return above
        else:
            return not above

    #--------------------------
    # create a new site where the _Halfedges el1 and el2 intersect
    def intersect(self, other):
        e1 = self.edge
        e2 = other.edge
        if (e1 is None) or (e2 is None):
            return None

        # if the two edges bisect the same parent return None
        if e1.reg[1] is e2.reg[1]:
            return None

        d = e1.a * e2.b - e1.b * e2.a
        if _float_eq(d, 0.0):
            return None

        xint = (e1.c*e2.b - e2.c*e1.b) / d
        yint = (e2.c*e1.a - e1.c*e2.a) / d
        if cmp(e1.reg[1], e2.reg[1]) < 0:
            he = self
            e = e1
        else:
            he = other
            e = e2

        rightOfSite = xint >= e.reg[1].x
        if ((rightOfSite and he.pm == _Edge.LE) or
            (not rightOfSite and he.pm == _Edge.RE)):
            return None

        # create a new site at the point of intersection - this is a new 
        # vector event waiting to happen
        return _Site(xint, yint)

#------------------------------------------------------------------
class _EdgeList(object):
    def __init__(self, xmin, xmax, nsites):
        if xmin > xmax:
            xmin, xmax = xmax, xmin
        self.hashsize = int(2 * math.sqrt(nsites + 4))
        self.xmin = xmin
        self.dx = xmax - xmin
        self.hash = [None] * self.hashsize
        self.leftend  = _Halfedge()
        self.rightend = _Halfedge()
        self.leftend.right = self.rightend
        self.rightend.left = self.leftend
        self.hash[0]  = self.leftend
        self.hash[-1] = self.rightend

    def insert(self,left,he):
        he.left  = left
        he.right = left.right
        left.right.left = he
        left.right = he

    def delete(self,he):
        he.left.right = he.right
        he.right.left = he.left
        he.edge = _Edge.DELETED

    def pop_leftbnd(self, pt):
        # Use hash table to get close to desired halfedge 
        bucket = int(((pt.x - self.xmin) / self.dx*self.hashsize))
        
        if bucket < 0: 
            bucket = 0;
        elif bucket >= self.hashsize: 
            bucket = self.hashsize - 1

        he = self._gethash(bucket)
        if he is None:
            i = 1
            while True:
                he = self._gethash(bucket-i)
                if he is not None:
                    break
                he = self._gethash(bucket+i)
                if he is not None:
                    break
                i += 1
    
        # Now search linear list of halfedges for the correct one
        if (he is self.leftend) or (he is not self.rightend and he.isPointRightOf(pt)):
            he = he.right
            while he is not self.rightend and he.isPointRightOf(pt):
                he = he.right
            he = he.left;
        else:
            he = he.left
            while he is not self.leftend and not he.isPointRightOf(pt):
                he = he.left

        if 0 < bucket < self.hashsize-1:
            self.hash[bucket] = he
            
        return he

    # Get entry from hash table, pruning any deleted nodes 
    def _gethash(self, b):
        if b < 0 or b >= self.hashsize:
            return None
        he = self.hash[b]
        if he is None or he.edge is not _Edge.DELETED:
            return he
        #  Hash table points to deleted half edge.  Patch as necessary.
        self.hash[b] = None
        return None



#------------------------------------------------------------------
class _PriorityQueue(object):
    def __init__(self, ymin, ymax, nsites):
        self.ymin = ymin
        self.deltay = ymax - ymin
        self.hashsize = int(4 * math.sqrt(nsites))
        self.count = 0
        self.minidx = 0
        self.hash = []
        for unused in range(self.hashsize):
            self.hash.append(_Halfedge())

    def __len__(self):
        return self.count

    def isEmpty(self):
        return self.count == 0

    def insert(self, he, site, offset):
        he.vertex = site
        he.ystar  = site.y + offset
        last = self.hash[self._get_bucket(he)]
        qnext = last.qnext
        while (qnext is not None) and cmp(he, qnext) > 0:
            last = qnext
            qnext = last.qnext
        he.qnext = last.qnext
        last.qnext = he
        self.count += 1

    def delete(self, he):
        if he.vertex is not None:
            last = self.hash[self._get_bucket(he)]
            while last.qnext is not he:
                last = last.qnext
            last.qnext = he.qnext
            he.vertex = None
            self.count -= 1

    def getMinPt(self):
        while self.hash[self.minidx].qnext is None:
            self.minidx += 1
        he = self.hash[self.minidx].qnext
        return _Site(he.vertex.x, he.ystar)

    def popMinHalfedge(self):
        curr = self.hash[self.minidx].qnext
        self.hash[self.minidx].qnext = curr.qnext
        self.count -= 1
        return curr

    def _get_bucket(self, he):
        bucket = int(((he.ystar - self.ymin) / self.deltay) * self.hashsize)
        if bucket < 0:
            bucket = 0
        if bucket >= self.hashsize:
            bucket = self.hashsize-1
        if bucket < self.minidx:
            self.minidx = bucket
        return bucket


#------------------------------------------------------------------
class _SiteList(list):
    """A sorted list of sites with min/max point values."""
    def __init__(self, points):
        """Points should be 2-tuples with x and y value"""
        self.xmin = sys.float_info.max
        self.ymin = sys.float_info.max
        self.xmax = sys.float_info.min
        self.ymax = sys.float_info.min
        for i, p in enumerate(points):
            site = _Site(p[0], p[1], i)
            self.append(site)
            self.xmin = min(site.x, self.xmin)
            self.ymin = min(site.y, self.ymin)
            self.xmax = max(site.x, self.xmax)
            self.ymax = max(site.y, self.ymax)
        self.sort()

def _float_eq(a, b):
    """Compare two floats for relative equality.
    
    See: http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
    for a discussion of floating point comparisons.
    """
    norm = max(abs(a), abs(b))
    return (norm < _EPSILON) or (abs(a - b) < (_EPSILON * norm))


