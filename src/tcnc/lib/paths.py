'''
Created on Apr 5, 2013

@author: claude
'''

import sys

import simplepath

import geom

class Path(list):
    """A list of Line, Arc, and CubicBeziers.
    """
    name_prefix = 'pathid_'
    name_autoindex = 1000
    transform = geom.IdentityTransform2D
    
    def __init__(self, name=None, transform=None):
        """:path_id: path identifier or label, such as an SVG path id.
        """
        self.name = name
        if self.name is None:
            self.name = '%s%d' % (Path.name_prefix, Path.name_autoindex)
            Path.name_autoindex += 1
        if transform is not None:
            self.transform = transform
    
    @property
    def startp(self):
        """The start point of this path."""
        return self[0].p1
    
    @property
    def endp(self):
        """The end point of this path."""
        return self[-1].p2
            
    def parse_svg(self, d):
        """Extract path segments from an SVG path specifier
        (the 'd' attribute value) and append them to this path."""
        
    def to_svg_segments(self):
        """Return path segments as SVG path segments (simplepath)."""
        sp = [('M', (self.startp.x, self.startp.y)),]
        for segment in self:
            sp.append(segment.to_svg_segment())
        return sp
    
    def to_svg_path(self):
        """Return path segments as SVG path specifier ('d' attribute value.)"""
        return simplepath.formatPath(self.to_svg_segments())
    
    def is_closed(self):
        """Return True if this path forms a closed polygon."""
        return len(self) > 2 and self.startp == self.endp
    
    def bounding_box(self, recalc=False):
        """Return the bounding rectangle aligned to the XY axes.
        The value is cached unless :recalc: is True.
        """
        if recalc or not hasattr(self, 'bbox'):
            min_x = sys.float_info.max
            max_x = 0.0
            min_y = sys.float_info.max
            max_y = 0.0
            for segment in self:
                min_x = min(min(segment.p1.x, segment.p2.x), min_x)
                min_y = min(min(segment.p1.y, segment.p2.y), min_y)
                max_x = max(max(segment.p1.x, segment.p2.x), max_x)
                max_y = max(max(segment.p1.y, segment.p2.y), max_y)
            self.bbox = geom.Box(geom.P(min_x, min_y), geom.P(max_x, max_y))
        return self.bbox
    
    
#    def convert_beziers_to_biarcs(self, tolerance=0.01, max_depth=4, line_flatness=0.01):
#        """Convert any CubicBeziers to Arcs using biarc approximation.
#        Returns a new Path.
#        """
#        newpath = Path(self.name, self.transform)
#        for shape in self:
#            if isinstance(shape, geom.CubicBezier):
#                biarcs = shape.biarc_approximation(tolerance=tolerance,
#                                                   max_depth=max_depth,
#                                                   line_flatness=line_flatness)
#                newpath.extend(biarcs)
#            else:
#                newpath.append(shape)
#        return newpath

