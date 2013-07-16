#-----------------------------------------------------------------------------#
#    Copyright (C) 2012,2013 Claude Zervas
#    email: claude@utlco.com
#    
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#    
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#-----------------------------------------------------------------------------#

"""A simple library for manipulating SVG.

Useful for creating Inkscape plugins

====
"""
from lxml import etree

def svgns(tag):
    return '{http://www.w3.org/2000/svg}%s' % tag

class SVGContext(object):
    """"""
    # Document root
    docroot = None
    # Current node under which new SVG elements will be added
    current_parent = None

    def __init__(self, docroot, parent=None, *args, **kwargs):
        self.docroot = docroot
        self.current_parent = docroot if parent is None else parent
        
    def create_path(self, attrs, style=None, parent=None):
        """Create an SVG path element."""
        if parent is None:
            parent = self.current_parent
        if parent is not None:
            if style:
                attrs['style'] = style
            return etree.SubElement(parent, svgns('path'), attrs)

    def create_circle(self, cx, cy, radius, style=None, parent=None):
        """Create an SVG circle element."""
        if parent is None:
            parent = self.current_parent
        if parent is not None:
            attrs = {'r': str(radius), 'cx': str(cx), 'cy': str(cy)}
            if style:
                attrs['style'] = style
            return etree.SubElement(parent, svgns('circle'), attrs)
    
    def create_line(self, x1, y1, x2, y2, style=None, parent=None):
        """Shortcut to create an SVG path consisting of one line segment."""
        attrs = {'d': 'M %5f %5f L %5f %5f' % (x1, y1, x2, y2)}
        return self.create_path(attrs, style, parent)

    def create_arc(self, x1, y1, x2, y2, radius, angle, style=None, parent=None):
        """Create an SVG arc."""
        sweep_flag = 0 if angle < 0 else 1
        attrs = { 'd': 'M %5f %5f A %5f %5f 0 0 %d %5f %5f' % \
                  (x1, y1, radius, radius,
                   sweep_flag, x2, y2),}
        return self.create_path(attrs, style, parent)

    def create_polygon(self, vertices, style=None, parent=None):
        """Create an SVG path describing a polygon."""
        d = 'M %f %f L' % (vertices[0].x, vertices[0].y)
        for p in vertices[1:]:
            d = d + ' %f,%f' % (p.x, p.y)
        if vertices[0] != vertices[-1]:
            d = d + ' %f,%f' % (vertices[0].x, vertices[0].y)
        attrs = {'d': d}
        return self.create_path(attrs, style, parent)

    def create_simple_marker(self, marker_id, d, style, transform):
        """Create an Inkscape line end marker glyph."""
        defs = self.docroot.find(svgns('defs'))
        if defs is None:
            defs = etree.SubElement(self.docroot, svgns('defs'))
        marker = etree.SubElement(defs, svgns('marker'),
                        {'id': marker_id, 'orient': 'auto', 'refX':  '0.0',
                         'refY': '0.0', 'style': 'overflow:visible',})
        etree.SubElement(marker, svgns('path'), 
                        { 'd': d, 'style': style, 'transform': transform, })
        return marker
    
