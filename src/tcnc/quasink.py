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

"""An Inkscape extension to create Penrose tilings.

====
"""
import os
import sys
import logging
import math
import gettext
_ = gettext.gettext

import inkex
import simplestyle

from lib import geom
from lib import supereffect
from lib import quasi
from lib import polygon
from lib import voronoi

VERSION = "0.1"

DEFAULT_LOG_LEVEL = logging.DEBUG
logger = logging.getLogger(__name__)

class IdentityProjector(object):
    """Identity projection. No distortion."""
    def project(self, p):
        return p
    
class SphericalProjector(IdentityProjector):
    """Project a point on to a sphere."""
    def __init__(self, center, radius, invert=False):
        self.center = center
        self.radius = radius
        self.invert = invert
    
    def project(self, p):
        v = p - self.center
        h = v.length()
        if h > self.radius:
            return p
        scale = math.sin((h * (math.pi / 2)) / self.radius)
        if not self.invert:
            scale = (self.radius * scale) / h
        return (v * scale) + self.center


class Quasink(supereffect.SuperEffect):
    """Inkscape plugin that creates quasi-crystal-like patterns.
    Based on quasi.c by Eric Weeks.
    """
    styles = {
        'polygon':
            'fill:none;stroke:#333333;stroke-width:0.5pt;',
        'polyseg':
            'fill:none;stroke-opacity:0.8;stroke-width:%s;stroke:%s;',
        'polysegpath':
            'fill:none;stroke-opacity:0.8;stroke-linejoin:round;'
            'stroke-width:%s;stroke:%s;',
        'segment':
            'fill:none;stroke-opacity:0.8;stroke-width:%s;stroke:%s;',
        'segpath':
            'fill:none;stroke-opacity:0.8;stroke-linejoin:round;'
            'stroke-width:%s;stroke:%s;',
        'voronoi':
            'fill:none;stroke-opacity:0.8;stroke-width:%s;stroke:%s;',
    }
    
    # Scale multiplier. This should be about right to get the whole
    # thing on a A4 size sheet of paper.
    _SCALE_SCALE = 8.0
    
    projector = None
        
    def effect(self):
        """Main entry point for Inkscape plugins.
        """
        if self.options.log_create_log and self.options.log_filename:
            log_path = os.path.expanduser(self.options.log_filename)
            log_path = os.path.abspath(log_path)
            log_level = getattr(logging, self.options.log_level,
                                DEFAULT_LOG_LEVEL)
            logging.basicConfig(filename=log_path, filemode='w', level=log_level)
            
#        logger.debug('scale: %f' % self.options.scale)
#        logger.debug('fillon: %s' % self.options.fillon)
#        logger.debug('zfill: %s' % self.options.zfill)
#        logger.debug('flip: %s' % self.options.flip)
#        logger.debug('symmetry: %d' % self.options.symmetry)
#        logger.debug('maxmax: %d' % self.options.numlines)
        logger.debug('midon0: %d' % int(self.options.mid_skinny))
        logger.debug('midon1: %d' % int(self.options.mid_fat))
#        logger.debug('polygon-stroke: %s' % self.options.polygon_stroke)
#        logger.debug('segment-stroke: %s' % self.options.segment_stroke)

        units = self.get_document_units()
        unit_scale = inkex.uuconv[units]
        # Perform any necessary unit conversion on plugin options
        self.convert_option_units(default_unit_scale=unit_scale)
        
        geom.set_epsilon(self.options.epsilon)
#        logger.debug('EPSILON= %f' % geom.EPSILON)

        self.scale = self.options.scale * self._SCALE_SCALE
        self.offset = geom.P(self.view_center) + geom.P(self.options.offset_x,
                                                        self.options.offset_y)
        mrot = geom.matrix_rotate_2d(self.options.rotate)
        self.transform = geom.matrix_scale_translate_2d(self.scale, self.scale,
                                                self.offset.x, self.offset.y)
        self.transform = geom.compose_transform(self.transform, mrot)       
        doc_size = geom.P(self.get_document_size())
        self.doc_center = doc_size / 2
                
        plotter = SVGQuasiPlotter()
#        if DEFAULT_LOG_LEVEL == logging.DEBUG:
#            geom.DEBUG_EFFECT = self
#            geom.DEBUG_LAYER = self.create_layer('debug')
#            plotter.debug_layer = geom.DEBUG_LAYER
        
        # Set up document and margin clipping regions
        mBL = geom.P(self.options.margin_left, self.options.margin_bottom)
        mTR = doc_size - geom.P(self.options.margin_right,
                                self.options.margin_top)
        margin_clip_rect = geom.Box(mBL, mTR)
        if self.options.clip_to_circle:
            r = min(doc_size.x, doc_size.y) / 2.0
            plotter.doc_clip_region = geom.Circle(self.doc_center, r)
            r = min(margin_clip_rect.width(), margin_clip_rect.height()) / 2.0
            plotter.margin_clip_region = geom.Circle(self.doc_center, r)
        else:
            if self.options.clip_to_doc:
                plotter.doc_clip_region = geom.Box(geom.P(0,0), doc_size)
            if self.options.clip_to_margins:
                plotter.margin_clip_region = margin_clip_rect
            
        plotter.doc_center = self.doc_center
        plotter.blowup_scale = self.options.blowup_scale
#        plotter.clip_to_doc = self.options.clip_to_doc
#        plotter.clip_polygons_to_margins = self.options.clip_poly
#        plotter.clip_segments_to_margins = self.options.clip_segments
            
        if self.options.project_sphere:
            self.projector = SphericalProjector(self.doc_center,
                                                self.options.project_radius,
                                                invert=self.options.project_invert)
        if self.projector is None:
            self.projector = IdentityProjector()
        plotter.projector = self.projector
        if self.options.polygon_mult == 1:
            plotter.polygon_layers.append(self.create_layer('quasink polygons'))
        else:
            for n in range(self.options.polygon_mult):
                layer = self.create_layer('quasink polygons %02d' % (n+1))
                plotter.polygon_layers.append(layer)
        if self.options.segment_draw:
            plotter.segment_layer = self.create_layer('quasink segments')
        plotter.scale = self.scale
        plotter.offset = self.offset
        plotter.transform = self.transform
#        plotter.offsetx = self.offset[0]
#        plotter.offsety = self.offset[1]
        plotter.fill_polygons = bool(self.options.polygon_fill)
        plotter.use_color = bool(self.options.polygon_colorfill)
        plotter.polygon_stroke = self.options.polygon_stroke
        plotter.segment_stroke = self.options.segment_stroke
        plotter.polygon_stroke_width = self.options.polygon_stroke_width
        plotter.segment_stroke_width = self.options.segment_stroke_width
        
        plotter.polygon_mult = self.options.polygon_mult
        plotter.polygon_mult_spacing = self.options.polygon_mult_spacing
        
        q = quasi.Quasi()
        q.offset_salt_x = self.options.salt_x
        q.offset_salt_y = self.options.salt_y
        q.skinnyfat_ratio = self.options.skinnyfat_ratio
        q.segment_ratio = self.options.segment_ratio
        q.segtype_skinny = self.options.mid_skinny
        q.segtype_fat = self.options.mid_fat
        q.symmetry = self.options.symmetry
        q.numlines = self.options.numlines
        q.plotter = plotter
        q.quasi()
        
        def segsort_key(seg):
#            mp = seg.midpoint() - self.doc_center
#            return mp.angle()
            return seg.midpoint().distance(self.doc_center)
                    
        if self.options.polyseg_draw:
            polyseg_layer = self.create_layer('quasink_polygon_segments')
            color = self.parse_color(self.options.polyseg_stroke)
            style = self.styles['polyseg'] % (self.options.polyseg_stroke_width,
                                              color,)
            for segment in plotter.poly_segments:
                self.draw_segment(segment, polyseg_layer, style)
                
        if self.options.polysegpath_draw:
            polysegpath_layer = self.create_layer('quasink_polygon_segment_paths')
            segtable = SegmentChainer2()
#            seglist = sorted(plotter.poly_segments, key=segsort_key)
            segpath_list = segtable.create_chains(plotter.poly_segments)
            # Sort segment paths so that the largest are at the bottom of the Z-order
#            segpath_list.sort(key=SegmentPath.bounding_box_area, reverse=True)
            segpath_list.sort(key=len, reverse=True)
            style = self.styles['polysegpath'] % (self.options.polysegpath_stroke_width,
                                                  self.parse_color(self.options.polysegpath_stroke),)
            for segpath in segpath_list:
                plotter.plot_segpath(polysegpath_layer, segpath, style)
                
#        for p in plotter.poly_segments.nodes.keys():
#            plotter.debug_plot_point(p)

        if self.options.convex_hull_draw:
            hull_layer = self.create_layer('quasink_convex_hull')
            plotter.polygon_layers[0] = hull_layer
            
            hull = plotter.poly_segments.convex_hull()
            plotter.polygon_stroke = "#f08080"
            plotter.plot_polygon(hull, nofill=True, draw_only=True)
        
        if self.options.hull_draw or self.options.voronoi_draw:
            polygon_hulls = plotter.poly_segments.boundary_polygons(
                                peel_hulls=self.options.hull_inner_draw)

        if self.options.hull_draw:
            hull_layer = self.create_layer('quasink_polygon_hulls')
            plotter.polygon_layers[0] = hull_layer
            for hull in polygon_hulls:
                plotter.polygon_stroke = self.options.hull_stroke
                plotter.polygon_stroke_width = self.options.hull_stroke_width
                plotter.plot_polygon(hull, nofill=True, draw_only=True)
                                
        if self.options.segpath_draw:
            segpath_layer = self.create_layer('quasink_segment_paths')
            segtable = SegmentChainer2()
            segpath_list = segtable.create_chains(plotter.path_segments)
            style = self.styles['segpath'] % (self.options.segpath_stroke_width,
                                              self.parse_color(self.options.segpath_stroke),)
            for segpath in segpath_list:
                plotter.plot_segpath(segpath_layer, segpath, style)
                
        if self.options.voronoi_draw:
            voronoi_layer = self.create_layer('quasink_voronoi')
#            voronoi_path_layer = self.create_layer('quasink_voronoi_paths')
            style = self.styles['voronoi'] % (self.options.voronoi_stroke_width,
                                              self.parse_color(self.options.voronoi_stroke),)
#            points = []
#            for p in plotter.poly_segments.nodes.keys():
#                pt = self.projector.project(p * self.scale + self.offset)
#                points.append(pt)
            points = plotter.poly_segments.nodes.keys()
            voronoi_segments = self.voronoi(points)
            phull = []
            for p in polygon_hulls[0]:
                phull.append(self.projector.project(p * self.scale + self.offset))
            for segment in voronoi_segments:
                p1 = self.projector.project(segment.p1 * self.scale + self.offset)
                p2 = self.projector.project(segment.p2 * self.scale + self.offset)
                line = geom.Line(p1, p2)
                line = margin_clip_rect.clip_line(line)
                if line and polygon.point_inside(phull, line.p1) and \
                polygon.point_inside(phull, line.p2):
                    attrs = {'d': 'M %5f %5f L %5f %5f' % \
                             (line.p1.x, line.p1.y, line.p2.x, line.p2.y),
                             'style': style}
                    inkex.etree.SubElement(voronoi_layer, inkex.addNS('path', 'svg'), attrs)
                    
    def draw_convex_hull(self, plotter):
        hull_layer = self.create_layer('quasink_convex_hull')
        plotter.polygon_layers[0] = hull_layer
        
        hull = plotter.poly_segments.convex_hull()
        plotter.polygon_stroke = "#f08080"
        plotter.plot_polygon(hull, nofill=True, draw_only=True)
                
                        
    def parse_color(self, color_text):
        """Parse the color text input from the extension dialog."""
        try:
            if color_text and color_text.lower() != 'none':
                return simplestyle.formatColoria(simplestyle.parseColor(color_text))
        except:
            pass
        return 'none'

    def draw_segment(self, segment, layer, style):
        p1 = self.projector.project(segment.p1 * self.scale + self.offset)
        p2 = self.projector.project(segment.p2 * self.scale + self.offset)
        attrs = {'d': 'M %5f %5f L %5f %5f' % (p1.x, p1.y, p2.x, p2.y)}
        attrs['style'] = style
        inkex.etree.SubElement(layer, inkex.addNS('path', 'svg'), attrs)
        
    def voronoi(self, points):
        """Create a Voronoi diagram of the given points.
        Returns a list of segments.
        """
        [vertices, unused, edges, unused] = voronoi.voronoi_diagram(points)
#        [vertices, unused, edges] = voronoi_inkscape.computeVoronoiDiagram(points)
        voronoi_segments = []
        for edge in edges:
            if edge[1] >= 0 and edge[2] >= 0:
                v1 = vertices[edge[1]]
                v2 = vertices[edge[2]]
                line = geom.Line(geom.P(v1[0], v1[1]), geom.P(v2[0], v2[1]))
                voronoi_segments.append(line)
        return voronoi_segments

class SegmentPath(object):
    """"""
    def __init__(self, segments):
        super(SegmentPath, self).__init__()
        self.segment_list = segments
        
    def len(self):
        return len(self.segment_list)
    
    def bounding_box_area(self):
        """Return the area of the bounding box"""
        bbox = self.bounding_box()
        return bbox.width() * bbox.height()
        
    def bounding_box(self):
        """Return the bounding rectangle aligned to the XY axes."""
        if getattr(self, 'bbox', None) is None:
            min_x = sys.float_info.max
            max_x = 0.0
            min_y = sys.float_info.max
            max_y = 0.0
            for segment in self.segment_list:
                min_x = min(min(segment.p1.x, segment.p2.x), min_x)
                min_y = min(min(segment.p1.y, segment.p2.y), min_y)
                max_x = max(max(segment.p1.x, segment.p2.x), max_x)
                max_y = max(max(segment.p1.y, segment.p2.y), max_y)
            self.bbox = geom.Box(geom.P(min_x, min_y), geom.P(max_x, max_y))
        return self.bbox
    
    
class SegmentSet(set):
    """A set of line segments.
    """
#    def __init__(self):
#        super(SegmentSet, self).__init__()
#        self.insertion_order = []
        
    def add(self, segment):
        if segment not in self:
            super(SegmentSet, self).add(segment)
            self._added_unique(segment)
            
    def add_polygon(self, vertices):
        """Add a polygon to the list.
        Constructs segments and closes polygon if necessary."""
        p1 = vertices[0]
        for p2 in vertices[1:]:
            self.add(geom.Line(geom.P(p1), geom.P(p2)))
            p1 = p2
        if vertices[0] != vertices[-1]:
            self.add(geom.Line(geom.P(vertices[-1]), geom.P(vertices[0])))
            
    def _added_unique(self, segment):
        """This is called when a new unique segment is added to this set.
        """
#        self.insertion_order.append(segment)
                
class SegmentGraph(SegmentSet):
    """A segment graph.
    """
    def __init__(self):
        super(SegmentGraph, self).__init__()
        self.nodes = {} # Hashtable of graph vertices
        self.sum_x = 0
        self.sum_y = 0
        # point with min y
        self.p_min_y = geom.P.Max()
        
    def order(self):
        """Number of nodes (vertices.)"""
        return len(self.nodes)

    def size(self):
        """Number of edges (segments.)"""
        return len(self)
    
    def centroid(self):
        k = len(self) * 2
        return geom.P(self.sum_x / k, self.sum_y / k)
        
    def _added_unique(self, segment):
        """Overrides SegmentSet..."""
        super(SegmentGraph, self)._added_unique(segment)
        # Build the node graph
        if segment.p1 in self.nodes:
            self.nodes[segment.p1].append(segment.p2)
        else:
            # New node from segment start point
            self.nodes[segment.p1] = [segment.p2,]
            self.sum_x += segment.p1.x
            self.sum_y += segment.p1.y
        if segment.p2 in self.nodes:
            self.nodes[segment.p2].append(segment.p1)
        else:
            # New node from segment end point
            self.nodes[segment.p2] = [segment.p1,]
            self.sum_x += segment.p2.x
            self.sum_y += segment.p2.y
        # Update min y node
        if segment.p1.y < self.p_min_y.y:
            self.p_min_y = segment.p1
        if segment.p2.y < self.p_min_y.y:
            self.p_min_y = segment.p2
        
    def convex_hull(self):
        """Return the convex hull of the endpoints of this segment graph.
        Returns a list of points.
        """
        return map(geom.P, polygon.convex_hull_graham_scan(self.nodes.keys()))
        
    def boundary_polygon(self):
        """Return a polygon defining the outer edges of this segment graph.
        """
        return self._calc_boundary_polygon(self.p_min_y, self.nodes)
        
    def boundary_polygons(self, peel_hulls=True):
        """Similar to convex hull peeling except using possibly
        non-convex boundary polygons."""
        nodes = dict(self.nodes)
        poly_list = []
        p_min_y = self.p_min_y
        while len(nodes) > 3:
            poly = self._calc_boundary_polygon(p_min_y, nodes)
            poly_list.append(poly)
            self._prune_nodes(nodes, poly)
            p_min_y = self._find_p_min_y(nodes)
            if len(poly) < 4 or not peel_hulls:
                break
        return poly_list
    
    def _prune_nodes(self, nodes, hullpoints):
        """Prune the nodes corresponding to the list of points
        in the outer polygon."""
        # Delete the outer nodes (polygon)
        for p in hullpoints:
            if p in nodes:# and len(nodes[p]) < 4:
                del nodes[p]
        # Remove any edges that have only one endpoint
        for node in nodes.values():
            pts = list(node)
            for p in pts:
                if p not in nodes:
                    node.remove(p)
        # Remove any spike nodes
        while True:
            nspikes = 0
            for p, node in nodes.items():
                if len(node) < 2:
                    nspikes += 1
                    if len(node) == 1:
                        nodes[node[0]].remove(p)
                    del nodes[p]
#            logger.debug('spikes: %d' % nspikes)
            if nspikes == 0:
                break
    
    def _find_p_min_y(self, nodes):
        """Given a collection of graph nodes,
        return the point that has the minimum Y value."""
        p_min_y = geom.P.Max()
        for p in nodes.keys():
            if p.y < p_min_y.y:
                p_min_y = p
        return p_min_y
        
    def _calc_boundary_polygon(self, start_node, nodes):
        """Return a polygon defining the outer edges of this segment graph.
        """
        node_p = start_node
        ref_p = geom.P(node_p.x + 1, node_p.y)
        poly = [node_p,]
        for unused in range(len(nodes)):
            next_node_p = self._get_ccw_edge(node_p, nodes[node_p], ref_p)
            ref_p = node_p
            node_p = next_node_p
            poly.append(node_p)
            if node_p == start_node:
                break
        return poly
        
    def _get_ccw_edge(self, node_p, edge_endpoints, ref_p):
        """Return the endpoint of the edge whose angle with repsect to ref_p
        is the most counterclockwise of all edges emanating from this node."""
        def polarkey(p):
            return node_p.ccw_angle2(ref_p, p)
        ns = sorted(edge_endpoints, key=polarkey)
        return ns[-1]
                         
class SegmentChain(list):
    """A connected chain of line segments."""
    min_corner_angle = math.pi * .3
    
#    def reverse(self):
#        super(SegmentChain, self).reverse()
#        # Reverse every segment in the list
#        for n in range(len(self)):
#            self[n] = self[n].reversed()
    
    @property        
    def startp(self):
        return self[0].p1

    @property
    def endp(self):
        return self[-1].p2
    
    def add_segment(self, segment):
        """Try to add the segment to the chain. Return True if successful other False."""
        if len(self) == 0 or segment.p1 == self.endp:
            return self._append_segment(segment)
#            self.append(segment)
        elif segment.p2 == self.endp:
            return self._append_segment(segment.reversed())
#            self.append(segment.reversed())
        elif segment.p2 == self.startp:
            return self._prepend_segment(segment)
#            self.insert(0, segment)
        elif segment.p1 == self.startp:
            return self._prepend_segment(segment.reversed())
#            self.insert(0, segment.reversed())
        else:
            return False
#        return True
    
    def _prepend_segment(self, segment):
        if len(self) == 0 or self._angle_is_ok(segment, self[0]):
            self.insert(0, segment)
            return True
        return False
    
    def _append_segment(self, segment):
        if len(self) == 0 or self._angle_is_ok(self[-1], segment):
            self.append(segment)
            return True
        return False
    
    def _angle_is_ok(self, seg1, seg2):
        if isinstance(seg1, geom.Line) and isinstance(seg2, geom.Line):
            a = abs(seg1.p2.angle2(seg1.p1, seg2.p2))
            return a > self.min_corner_angle
        return True
    

class SegmentChainer2(object):
    """A list of line segments chains."""    
    def __init__(self):
        self.chain_list = []
        
    def create_chains(self, segments):
        while segments:
            chain = SegmentChain()
            n = 1
            while n > 0:
                unchained_segments = []
                for segment in segments:
                    if not chain.add_segment(segment):
                        unchained_segments.append(segment)
                n = len(segments) - len(unchained_segments)
                segments = unchained_segments
            self.chain_list.append(chain)
        return self.chain_list
    

class Document(object):
    def __init__(self, transform=geom.IdentityTransform2D, clip_region=None,
                 projector=IdentityProjector()):
        self.polygons = []
        self.segments = []
        self.transform = transform
        self.clip_region = clip_region
        self.projector = projector
        self._style = None
        self._style_changed = True
        
    def plot_polygon(self, vertices):
        xvertices = self.transform_points(vertices)
        clipped = self._is_clipped(xvertices)
        if not clipped:
            self.polygons.append(xvertices)
        return clipped

    def plot_segment(self, x1, y1, x2, y2):
        endpoints = geom.Line(self.transform_point(geom.P(x1, y1)), 
                              self.transform_point(geom.P(x2, y2)))
        self.segments.append(endpoints)
    
    def transform_points(self, points):
        return map(lambda p: self.transform_point(p), points)
    
    def transform_point(self, p):
        return self.projector.project(p.transform(self.transform))
                    
    def _is_clipped(self, points):
        return self.clip_region is not None and \
            not self.clip_region.all_points_inside(points)    
    
        
class SimpleQuasiPlotter(quasi.QuasiPlotter):
    """Just accumulates the original quasi geometry with no clipping.
    """
    def __init__(self):
        self.polygons = []
        self.segments = []
        
    def plot_polygon(self, vertices):
        self.polygons.append(map(geom.P, vertices))
        return False

    def plot_segment(self, x1, y1, x2, y2):
        self.segments.append(geom.Line(geom.P(x1, y1), geom.P(x2, y2)))
    

class SVGPlotter(object):
    """Plot segments and polygons.
    """
    STYLE = {'fill': 'none',
             'stroke-opacity': '0.8',
             'stroke-linejoin': 'round',
             'stroke-width': '1pt',
             'stroke': '#000000'}
    
    layer = None

    @property
    def stroke_width(self):
        return self.STYLE['stroke-width']
    
    @property
    def stroke_color(self):
        return self.STYLE['stroke']
    
    @property
    def stroke_opacity(self):
        return self.STYLE['stroke-opacity']
    
    @property
    def stroke_linejoin(self):
        return self.STYLE['stroke-linejoin']
    
    def __init__(self, effect):
        self.effect = effect
        
    def draw_segment(self, p1, p2):
        self.effect.create_line(p1.x, p1.y, p2.x, p2.y, self._get_style(), self.layer)
        
    def draw_polygon(self, vertices):
        self.effect.create_polygon(vertices, self._get_style(), self.layer)
        
    def _get_style(self):
        return simplestyle.formatStyle(self.STYLE)

class SVGQuasiPlotter(quasi.QuasiPlotter):
    """SVG output for quasi"""
    scale = 1.0
    offset = geom.P(0, 0)
    transform = geom.IdentityTransform2D
    blowup_scale = 1.0
    doc_center = geom.P(0, 0)
    debug_layer = None
    polygon_layers = []
    segment_layer = None
    polyseg_layer = None
    doc_clip_region = None
    margin_clip_region = None
#    clip_polygons_to_margins = False
#    clip_segments_to_margins = False
#    clip_to_doc = True
    fill_polygons = False
    polygon_stroke = '#333333'
    segment_stroke = '#666666'
    polygon_stroke_width = '.5pt'
    segment_stroke_width = '.5pt'
    projector = IdentityProjector()
    path_segments = SegmentSet()
    poly_segments = SegmentGraph()
    
    polygon_mult = 1
    polygon_mult_spacing = 0.0

    def plot_segpath(self, layer, segpath, style=None):
        p = self.projector.project(self._transform_point(segpath[0].p1))
        d = 'M %f %f L' % (p.x, p.y)
        for segment in segpath:
            p = self.projector.project(self._transform_point(segment.p2))
            d = d + ' %f,%f' % (p.x, p.y)
        attrs = {'d': d, 'style': style}
        inkex.etree.SubElement(layer, inkex.addNS('path', 'svg'), attrs)
    
    def plot_polygon(self, vertices, nofill=False, draw_only=False):
        """Polygons are always four sided..."""
        if not self.polygon_layers:
            return
        # Test for degenerate rhombus
        if vertices[0] == vertices[2] or vertices[1] == vertices[3]:
            return True
        xvertices = [self._transform_point(geom.P(p)) for p in vertices]
        polygon_clipped = self._is_clipped(xvertices)
        if not polygon_clipped:
            xvertices = [self._transform_point(geom.P(p)) for p in vertices]
            xvertices = self._blowout_polygon(xvertices)
            xvertices = [self.projector.project(p) for p in xvertices]
            polygon_clipped = self._is_clipped(xvertices)
            if not polygon_clipped:
                fill = not nofill and self.fill_polygons
                self._draw_polygon(xvertices, fill=fill)
                # Draw concentric polygons if any
                for n in range(1, self.polygon_mult):
                    xvertices = self._inset_polygon(xvertices, self.polygon_mult_spacing)
                    if xvertices is None:
                        break
                    self._draw_polygon(xvertices, fill=fill, nlayer=n)
                if not draw_only:
                    self.poly_segments.add_polygon(vertices)
        return polygon_clipped

    def plot_segment(self, x1, y1, x2, y2):
        if self.segment_layer is None:
            return
        segment = geom.Line(geom.P(x1, y1), geom.P(x2, y2))
        p1 = self._transform_point(segment.p1)
        p2 = self._transform_point(segment.p2)
        if not self._is_clipped((p1, p2)):
            p1 = self.projector.project(p1)
            p2 = self.projector.project(p2)
            attrs = {'d': 'M %5f %5f L %5f %5f' % (p1.x, p1.y, p2.x, p2.y)}
            attrs['style'] = 'fill:none;stroke:%s;stroke-width:%s;' % (self.segment_stroke, self.segment_stroke_width)
            inkex.etree.SubElement(self.segment_layer, inkex.addNS('path', 'svg'), attrs)
            self.path_segments.add(segment)

    def debug_plot_point(self, p, color="#000000"):
        if self.debug_layer is not None:
            p = self._transform_point(p)
            p = self.projector.project(p)
            attrs = {'r': '5pt', 'cx': str(p.x), 'cy': str(p.y),
                     'style': 'fill:%s;stroke:none' % (color,),}
            inkex.etree.SubElement(self.debug_layer, inkex.addNS('circle', 'svg'), attrs)
        
    def rgb2css(self, rgb):
        """convert the rgb float tuple to a CSS compatible color string"""
        return '#%02x%02x%02x' % (int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))

    def _is_clipped(self, points):
        """Return True if any of the points are clipped by the clipping region."""
        return (self.doc_clip_region is not None and not (self.doc_clip_region.all_points_inside(points))) \
        or (self.margin_clip_region is not None and not (self.margin_clip_region.all_points_inside(points)))
        
    def _transform_point(self, p):
        return (p * self.scale) + self.offset
#        return p.transform(self.transform)
    
    def _blowout_polygon(self, vertices):
        if geom.float_eq(self.blowup_scale, 1.0):
            return vertices
        centroid = self._polygon_centroid(vertices) - self.doc_center
        offset = ((centroid * self.blowup_scale) - centroid)
        return [geom.P(v) + offset for v in vertices]
    
    def _draw_polygon(self, vertices, fill=False, nlayer=0):
        d = 'M %f %f L' % (vertices[0].x, vertices[0].y)
        for p in vertices[1:]:
            d = d + ' %f,%f' % (p.x, p.y)
        if vertices[0] != vertices[-1]:
            d = d + ' %f,%f' % (vertices[0].x, vertices[0].y)
        attrs = {'d': d}
        fillcolor = 'none' if not fill else self.rgb2css(self.get_fill_color_rgb())
        attrs['style'] = 'stroke-opacity:.8;fill:%s;stroke:%s;stroke-width:%s;stroke-linejoin:round;' % \
                          (fillcolor, self.polygon_stroke, self.polygon_stroke_width)
        inkex.etree.SubElement(self.polygon_layers[nlayer], inkex.addNS('path', 'svg'), attrs)
    
    def _inset_polygon(self, vertices, offset):
        """Inset the polygon by the amount :offset:"""
        L1 = geom.Line(vertices[0], vertices[1])
        L2 = geom.Line(vertices[1], vertices[2])
        L3 = geom.Line(vertices[2], vertices[3])
        L4 = geom.Line(vertices[3], vertices[0])
        d1 = L1.distance_to_point(L3.p1)
        d2 = L2.distance_to_point(L1.p1)
        if d1 >= offset*3 and d2 >= offset*3:
            offset *= L1.which_side(L2.p2)
            L1_o = L1.offset(offset)
            L2_o = L2.offset(offset)
            L3_o = L3.offset(offset)
            L4_o = L4.offset(offset)
            p1 = L4_o.intersection(L1_o)
            p2 = L1_o.intersection(L2_o)
            p3 = L2_o.intersection(L3_o)
            p4 = L3_o.intersection(L4_o)
            return (p1, p2, p3, p4)
        return None
    
    def _polygon_centroid(self, vertices):
        L1 = geom.Line(vertices[0], vertices[2])
        L2 = geom.Line(vertices[1], vertices[3])
        return L1.intersection(L2)
    
    
    _number = 0
    def _draw_number(self, vertices):
        centroid = self._polygon_centroid(vertices)
        attrs = {'x': str(centroid.x), 'y': str(centroid.y),
                 inkex.addNS('space','xml'): 'preserve',
                 'style': 'font-size:24px;font-family:Sans;fill:#000000'}
        #inkex.addNS('space','xml'): 'preserve'
        t = inkex.etree.SubElement(self.polygon_layers[0], inkex.addNS('text', 'svg'), attrs)
        attrs = {'x': str(centroid.x), 'y': str(centroid.y), inkex.addNS('role','sodipodi'):'line'}
        span = inkex.etree.SubElement(t, inkex.addNS('tspan', 'svg'), attrs)
        span.text = str(self._number)
        self._number += 1

option_info = (
    supereffect.optargs('--active-tab',),
    
    supereffect.optargs('--scale', '-s', type='float', default=1.0, help='Output scale.'),
    supereffect.optargs('--rotate', '-r', type='float', convert_to='rad', default=0.0, help='Rotation.'),
#    supereffect.optargs('--flip', '-F', type='inkbool', default=False, help='Flip by 90 deg.'),
    supereffect.optargs('--symmetry', '-S', type='int', default=5, help='Degrees of symmetry.'),
    supereffect.optargs('--numlines', '-n', type='int', default=30, help='Number of lines.'),
    supereffect.optargs('--mid-skinny', '-M', type='int', default=1, help='Midpoint type for skinny diamonds.'),
    supereffect.optargs('--mid-fat', '-N', type='int', default=1, help='Midpoint type for fat diamonds.'),
    supereffect.optargs('--skinnyfat-ratio', type='float', default=0.2, help='Skinny/fat ratio'),
    supereffect.optargs('--segment-ratio', type='float', default=0.5, help='Segment ratio'),
    supereffect.optargs('--offset-x', type='float', convert_to='world', default=0.0, help='X offset'),
    supereffect.optargs('--offset-y', type='float', convert_to='world', default=0.0, help='Y offset'),
    supereffect.optargs('--salt-x', type='float', default=0.1132, help='X offset salt'),
    supereffect.optargs('--salt-y', type='float', default=0.2137, help='Y offset salt'),
    supereffect.optargs('--epsilon', type='float', default=0.0001, help='Epsilon'),

    supereffect.optargs('--polygon-draw', type='inkbool', default=True, help='Draw polygons.'),
    supereffect.optargs('--polygon-mult', type='int', default=1, help='Number of concentric polygons.'),
    supereffect.optargs('--polygon-mult-spacing', type='float', convert_to='world', default=0.0, help='Concentric polygon spacing.'),
    supereffect.optargs('--polygon-fill', '-f', type='inkbool', default=False, help='Fill polygons.'),
    supereffect.optargs('--polygon-colorfill', type='inkbool', default=False, help='Use color fill.'),
    supereffect.optargs('--polygon-zfill', '-z', type='inkbool', default=False, help='Fill color is according to polygon type.'),
    supereffect.optargs('--polygon-stroke', default='#0000ff', help='Polygon CSS stroke color.'),
    supereffect.optargs('--polygon-stroke-width', default='.5pt', help='Polygon CSS stroke width.'),

    supereffect.optargs('--polyseg-draw', type='inkbool', default=True, help='Create segment-only polygon layer.'),
    supereffect.optargs('--polyseg-stroke', default='#0000ff', help='Polygon CSS stroke color.'),
    supereffect.optargs('--polyseg-stroke-width', default='.5pt', help='Polygon CSS stroke width.'),

    supereffect.optargs('--polysegpath-draw', type='inkbool', default=True, help='Create paths from polygon segments.'),
    supereffect.optargs('--polysegpath-stroke', default='#0000ff', help='Polygon CSS stroke color.'),
    supereffect.optargs('--polysegpath-stroke-width', default='.5pt', help='Polygon CSS stroke width.'),

    supereffect.optargs('--segment-draw', type='inkbool', default=True, help='Draw segments.'),
    supereffect.optargs('--segment-stroke', default='#0000ff', help='Segment CSS stroke color.'),
    supereffect.optargs('--segment-stroke-width', default='.5pt', help='Segment CSS stroke width.'),

    supereffect.optargs('--segpath-draw', type='inkbool', default=True, help='Draw segment paths.'),
    supereffect.optargs('--segpath-min-segments', '-m', type='int', default=1, help='Min segments in path.'),
    supereffect.optargs('--segpath-stroke', default='#0000ff', help='Segment CSS stroke color.'),
    supereffect.optargs('--segpath-stroke-width', default='.5pt', help='Segment CSS stroke width.'),

    supereffect.optargs('--convex-hull-draw', type='inkbool', default=True, help='Draw convex hull.'),
    supereffect.optargs('--hull-draw', type='inkbool', default=True, help='Draw polygon hull.'),
    supereffect.optargs('--hull-inner-draw', type='inkbool', default=True, help='Draw inner polygon hulls.'),
    supereffect.optargs('--hull-stroke', default='#0000ff', help='Polygon CSS stroke color.'),
    supereffect.optargs('--hull-stroke-width', default='.5pt', help='Polygon CSS stroke width.'),

    supereffect.optargs('--voronoi-draw', type='inkbool', default=True, help='Draw Voronoi.'),
    supereffect.optargs('--voronoi-clip-to-hull', type='inkbool', default=True, help='Clip to hull ploygon.'),
    supereffect.optargs('--voronoi-make-paths', type='inkbool', default=True, help='Make paths from segments.'),
    supereffect.optargs('--voronoi-stroke', default='#0000ff', help='Voronoi CSS stroke color.'),
    supereffect.optargs('--voronoi-stroke-width', default='.5pt', help='Voronoi CSS stroke width.'),

    supereffect.optargs('--clip-to-doc', type='inkbool', default=False, help='Clip to document.'),
    supereffect.optargs('--clip-to-circle', type='inkbool', default=False, help='Circular clip region.'),
    supereffect.optargs('--clip-to-margins', '-C', type='inkbool', default=True, help='Clip to document margins.'),
#    supereffect.optargs('--clip-segments', type='inkbool', default=True, help='Clip segments to document margins.'),
    
    supereffect.optargs('--margin-left', type='float',   convert_to='world', default=1.0, help='Left margin'),
    supereffect.optargs('--margin-right', type='float',  convert_to='world', default=1.0, help='Right margin'),
    supereffect.optargs('--margin-top', type='float',   convert_to='world', default=1.0, help='Top margin'),
    supereffect.optargs('--margin-bottom', type='float',   convert_to='world', default=1.0, help='Bottom margin'),

    supereffect.optargs('--project-sphere', type='inkbool', default=False, help='Project on to sphere.'),
    supereffect.optargs('--project-invert', type='inkbool', default=False, help='Invert projection.'),
    supereffect.optargs('--project-radius-useclip', type='inkbool', default=True, help='Use clipping circle for radius.'),
    supereffect.optargs('--project-radius', type='float',   convert_to='world', default=0.0, help='Projection radius.'),
    supereffect.optargs('--blowup-scale', type='float', default=1.0, help='Blow up scale.'),

    supereffect.optargs('--log-create-log', '-D', type='inkbool', default=True, help='Create log files'),
    supereffect.optargs('--log-level', '-L', default='DEBUG', help='Log level'),
    supereffect.optargs('--log-filename', '-O', default='quasink.log', help='Full pathname of log file'),
)

quasink = Quasink(option_info=option_info)
quasink.affect()
