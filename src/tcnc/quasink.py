'''An Inkscape extension to create Penrose tiles.
Based on code by Eric R. Weeks
See http://www.physics.emory.edu/~weeks/software/quasic.html

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
'''
import os
import sys
import logging
import math
import gettext
_ = gettext.gettext

import inkex
import simplestyle

from lib import geom
from lib import svg
from lib.geom import P, Line

import quasi

VERSION = "0.1"

DEFAULT_LOG_LEVEL = logging.DEBUG
logger = logging.getLogger(__name__)

class Quasink(svg.SuperEffect):
    '''Inkscape plugin that creates quasi-crystal-like patterns.
    Based on quasi.c by Eric Weeks.
    '''
    styles = {
              'segment': 'fill:none;stroke:#000000;stroke-width:0.35pt;',
              'segpath': 'fill:none;stroke:#ff0000;stroke-width:0.5pt;',
              'segpath1': 'fill:none;stroke:#00ff00;stroke-width:0.5pt;',
              'segpath_open': 'fill:none;stroke:#000ac0;stroke-width:1pt;',
              'segpath_closed': 'fill:none;stroke:#c00a00;stroke-width:1pt;',
              'polygon': 'fill:none;stroke:#333333;stroke-width:0.5pt;',
              }
    
    # Scale multiplier. This should be about right to get the whole
    # thing on a A4 size sheet of paper.
    SCALE_SCALE = 8.0
        
    def effect(self):
        '''Main entry point for Inkscape plugins.
        '''
        if self.options.log_create_log and self.options.log_filename:
            log_path = os.path.abspath(os.path.expanduser(self.options.log_filename))
            log_level = getattr(logging, self.options.log_level, DEFAULT_LOG_LEVEL)
            logging.basicConfig(filename=log_path, filemode='w', level=log_level)
            
#        logger.debug('scale: %f' % self.options.scale)
#        logger.debug('fillon: %s' % self.options.fillon)
#        logger.debug('zfill: %s' % self.options.zfill)
#        logger.debug('flip: %s' % self.options.flip)
#        logger.debug('symmetry: %d' % self.options.symmetry)
#        logger.debug('maxmax: %d' % self.options.numlines)
#        logger.debug('midon0: %d' % int(self.options.mid_skinny))
#        logger.debug('midon1: %d' % int(self.options.mid_fat))
#        logger.debug('polygon-stroke: %s' % self.options.polygon_stroke)
#        logger.debug('segment-stroke: %s' % self.options.segment_stroke)

        units = self.get_document_units()
        unit_scale = inkex.uuconv[units]
        self.options.margin_right *= unit_scale
        self.options.margin_left *= unit_scale
        self.options.margin_top *= unit_scale
        self.options.margin_bottom *= unit_scale

        self.scale = self.options.scale * self.SCALE_SCALE
        self.offset = P(self.view_center)
                
        plotter = SVGQuasiPlotter()
        top_right = P(self.get_document_size())
        plotter.doc_clip_rect = geom.Rectangle(P(0,0), top_right)
        bottom_left = P(self.options.margin_left, self.options.margin_bottom)
        top_right -= P(self.options.margin_right, self.options.margin_top)
        plotter.margin_clip_rect = geom.Rectangle(bottom_left, top_right)
        if self.options.draw_polygons:
            plotter.polygon_layer = self.create_layer('quasink_polygons')
        if self.options.draw_segments:
            plotter.segment_layer = self.create_layer('quasink_segments')
        if self.options.make_paths:
            plotter.segpath_layer1 = self.create_layer('quasink_segpaths_open')
            plotter.segpath_layer2 = self.create_layer('quasink_segpaths_closed')
        if self.options.polygon_segments:
            plotter.polyseg_layer = self.create_layer('quasink_polysegs')
        plotter.scale = self.scale
        plotter.offset = self.offset
#        plotter.offsetx = self.offset[0]
#        plotter.offsety = self.offset[1]
        plotter.clip_polygons_to_margins = bool(self.options.clip_poly)
        plotter.clip_segments_to_margins = bool(self.options.clip_segments)
        plotter.fill_polygons = bool(self.options.fillon)
        plotter.use_color = bool(self.options.colorfill)
        plotter.polygon_stroke = self.options.polygon_stroke
        plotter.segment_stroke = self.options.segment_stroke
        plotter.polygon_stroke_width = self.options.polygon_stroke_width
        plotter.segment_stroke_width = self.options.segment_stroke_width
        quasi.quasi(bool(self.options.zfill),
                     [int(self.options.mid_skinny), int(self.options.mid_fat)],
                     int(self.options.symmetry),
                     int(self.options.numlines),
                     plotter)
        
        if self.options.polygon_segments:
            style_polyseg = simplestyle.formatStyle({ 'fill': 'none',
                                 'stroke': self.parse_color(self.options.polygon_stroke),
                                 'stroke-width': self.options.polygon_stroke_width,
                                 'stroke-linejoin': 'round',
                                 })
            for segchain in plotter.poly_segments.get_segchains():
                for segment in segchain:
                    self.draw_segment(segment, plotter.polyseg_layer, style_polyseg)
#            for segment in plotter.poly_segments.get_segments():
#                self.draw_segment(segment, plotter.polyseg_layer)
                
        if self.options.make_paths:
            style_open_path = simplestyle.formatStyle({ 'fill': 'none',
                                 'stroke': self.parse_color(self.options.open_path_stroke),
                                 'stroke-width': self.options.open_path_stroke_width,
                                 'stroke-linecap': 'butt',
                                 'stroke-linejoin': 'round',
                                 })
            style_closed_path = simplestyle.formatStyle({ 'fill': self.parse_color(self.options.closed_path_fill),
                                 'stroke': self.parse_color(self.options.closed_path_stroke),
                                 'stroke-width': self.options.closed_path_stroke_width,
                                 'stroke-linejoin': 'round',
                                 })
            segment_chains = SegmentTable()
            for segment in plotter.path_segments:
                segment_chains.chain_segment(segment)
            segpath_list = segment_chains.get_segpaths()
            # Sort segment paths so that the largest are at the bottom of the Z-order
            segpath_list.sort(key=SegmentPath.bounding_box_area, reverse=True)
            for segpath in segpath_list:
                segchain = segpath.segment_list
                if len(segchain) >= self.options.minsegments:
                    if segchain[0].p1 == segchain[-1].p2:
                        plotter.plot_segpath(plotter.segpath_layer2, segchain, style_closed_path)
                    elif not self.options.closed_paths:
                        plotter.plot_segpath(plotter.segpath_layer1, segchain, style_open_path)
                        
    def parse_color(self, color_text):
        '''Parse the color text input from the extension dialog.'''
        try:
            if color_text and color_text.lower() != 'none':
                return simplestyle.formatColoria(simplestyle.parseColor(color_text))
        except:
            pass
        return 'none'

    def draw_segment(self, segment, layer, style):
        p1 = segment.p1 * self.scale + self.offset
        p2 = segment.p2 * self.scale + self.offset
        attrs = {'d': 'M %5f %5f L %5f %5f' % (p1.x, p1.y, p2.x, p2.y)}
        attrs['style'] = style
        inkex.etree.SubElement(layer, inkex.addNS('path', 'svg'), attrs)

class SegmentPath(object):
    ''''''
    def __init__(self, segments):
        super(SegmentPath, self).__init__()
        self.segment_list = segments
    
    def bounding_box_area(self):
        '''Return the area of the bounding box'''
        bbox = self.bounding_box()
        return bbox.width() * bbox.height()
        
    def bounding_box(self):
        '''Return the bounding rectangle aligned to the XY axes.'''
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
            self.bbox = geom.Rectangle(P(min_x, min_y), P(max_x, max_y))
        return self.bbox

class SegmentTable(object):
    '''A set of hashtables used to chain line segments
    and remove duplicates.'''    
    def __init__(self):
        super(SegmentTable, self).__init__()
        self.chain_startp = {}   # Hashtable of segment chain start coordinates
        self.chain_endp = {}     # Hashtable of segment chain end coordinates

    def add_polygon(self, vertices):
        '''Add polygon segments to the hashtable and discard duplicates.'''
        p1 = None
        p2 = P(vertices[0])
        for p in vertices[1:]:
            p1 = p2
            p2 = P(p)
            segset1 = self.chain_startp.get(p1)
            segset2 = self.chain_endp.get(p1)
            if not self.match_segment(p2, 1, segset1) and not self.match_segment(p2, 0, segset2):
                segment = Line(p1, p2)
                if segset1 is None:
                    self.chain_startp[p1] = [segment,]
                else:
                    segset1.append(segment)
                segset2 = self.chain_endp.get(p2)
                if segset2 is None:
                    self.chain_endp[p2] = [segment,]
                else:
                    segset2.append(segment)
    
    def match_segment(self, p, i, segments):
        '''Compare point <p> with point [i] in each segment of <segments>'''
        if segments:
            for seg in segments:
                if p == seg[i]:
                    return True
        return False
    
    def chain_segment(self, segment):
        '''Add a segment.
        Tries to add it to an existing segment chain.
        This will build polygons or open paths on the fly.
        '''
        # Does this segment match the start point of an existing chain but reversed?
        if segment.p1 in self.chain_startp:
            seg_first = self.chain_startp[segment.p1][0]
            if segment.p2 == seg_first.p2:
#                logger.debug('duplicate segment: %s' % (str(segment),))
                # Discard duplicate segments
                return
            else:
#                logger.debug('reversing segment')
                segment = segment.reverse()
                # Does thus segment match another chain in reverse?
                if segment.p1 in self.chain_startp:
#                    logger.debug('reversing chain 1')
                    # Pop the chain off the hash tables
                    chain = self.chain_startp.pop(segment.p1)
                    del self.chain_endp[chain[-1].p2]
                    self.reverse_chain(chain)
        # Does this segment match the end point of an existing chain but reversed?
        elif segment.p2 in self.chain_endp:
            seg_last = self.chain_endp[segment.p2][-1]
            if seg_last.p1 == segment.p1:
#                logger.debug('duplicate segment: %s' % (str(segment),))
                return
            else:
#                logger.debug('reversing segment')
                segment = segment.reverse()
                if segment.p2 in self.chain_endp:
#                    logger.debug('reversing chain 2')
                    # Pop the chain off the hash tables
                    chain = self.chain_endp.pop(segment.p2)
                    del self.chain_startp[chain[0].p1]
                    self.reverse_chain(chain)
        # Does this segment match the start of an existing chain?
        if segment.p2 in self.chain_startp:
#            logger.debug('matches start chain')
            chain1 = self.chain_startp.pop(segment.p2)
            chain1.insert(0, segment)
            # Can it link two existing chains?
            if segment.p1 in self.chain_endp:
#                logger.debug('linking chains 1')
                chain2 = self.chain_endp.pop(segment.p1)
                chain2.extend(chain1)
                self.chain_startp[chain2[0].p1] = chain2
                self.chain_endp[chain2[-1].p2] = chain2
            else:
                self.chain_startp[segment.p1] = chain1
        # Does the start point of this segment match the end point of an existing chain?
        elif segment.p1 in self.chain_endp:
#            logger.debug('matches end chain')
            chain1 = self.chain_endp.pop(segment.p1)
            chain1.append(segment)
            # Can it link another chain
            if segment.p2 in self.chain_startp:
#                logger.debug('linking chains 2')
                chain2 = self.chain_startp.pop(segment.p2)
                chain1.extend(chain2)
                self.chain_startp[chain1[0].p1] = chain1
                self.chain_endp[chain1[-1].p2] = chain1
            else:
                self.chain_endp[segment.p2] = chain1
        else:
#            logger.debug('inserting segment')
            chain = [segment,]
            self.chain_startp[segment.p1] = chain
            self.chain_endp[segment.p2] = chain
    
    def reverse_chain(self, chain):
        '''Reverse the chain and re-insert into hashtable'''
        chain_rev = []
        for segment in chain:
            chain_rev.insert(0, segment.reverse())
        # And re-insert into hash table
        self.chain_startp[chain_rev[0].p1] = chain_rev
        self.chain_endp[chain_rev[-1].p2] = chain_rev
        
    def get_segchains(self):
        return self.chain_startp.values()
    
    def get_segpaths(self):
        #return self.paths
        return map(SegmentPath, self.chain_startp.itervalues())
    
    
class BufQuasiPlotter(quasi.QuasiPlotter):
    '''A plotter that just accumulates the quasi geometry.'''
    polygons = []
    segments = []
    
    def plot_polygon(self, vertices):
        self.polygons.append(vertices)

    def plot_segment(self, x1, y1, x2, y2):
        self.segments.append(Line(P(x1, y1), P(x2, y2)))


class SVGQuasiPlotter(quasi.QuasiPlotter):
    '''SVG output for quasi'''
    scale = 1.0
    offset = P(0, 0)
#    offsetx = 0.0
#    offsety = 0.0
    polygon_layer = None
    segment_layer = None
    segpath_layer1 = None
    segpath_layer2 = None
    polyseg_layer = None
    doc_clip_rect = [P(0, 0), P(1.0, 1.0)]
    margin_clip_rect = [P(0, 0), P(1.0, 1.0)]
    clip_polygons_to_margins = False
    clip_segments_to_margins = False
    clip_polygons_to_doc = True
    clip_segments_to_doc = True
    segment_chains = SegmentTable()
    path_segments = set()
    poly_segments = SegmentTable()
    fill_polygons = False
    polygon_stroke = '#333333'
    segment_stroke = '#666666'
    polygon_stroke_width = '.5pt'
    segment_stroke_width = '.5pt'
    
    def plot_segpath(self, layer, segpath, style=None):
        p = (segpath[0].p1 * self.scale) + self.offset
        d = 'M %f %f L' % (p.x, p.y)
        for segment in segpath:
            p = (segment.p2 * self.scale) + self.offset
            d = d + ' %f,%f' % (p.x, p.y)
        attrs = {'d': d, 'style': style}
        inkex.etree.SubElement(layer, inkex.addNS('path', 'svg'), attrs)
    
    def plot_polygon(self, vertices):
        ''''''
        p = (P(vertices[0]) * self.scale) + self.offset
        if (not self.clip_polygons_to_doc or self.doc_clip_rect.point_inside(p)) \
        and (not self.clip_polygons_to_margins or self.margin_clip_rect.point_inside(p)):
            d = 'M %f %f L' % (p.x, p.y)
            for p in vertices[1:]:
                p = (P(p) * self.scale) + self.offset
                if (self.clip_polygons_to_doc and not self.doc_clip_rect.point_inside(p)) \
                or (self.clip_polygons_to_margins and not self.margin_clip_rect.point_inside(p)):
                    return
                d = d + ' %f,%f' % (p.x, p.y)
            attrs = {'d': d}
            if self.polygon_layer is not None:
                if self.fill_polygons:
                    fill = self.rgb2css(self.get_fill_color_rgb())
                else:
                    fill = 'none'
                attrs['style'] = 'fill:%s;stroke:%s;stroke-width:%s;stroke-linejoin:round;' % \
                                  (fill, self.polygon_stroke, self.polygon_stroke_width)
                inkex.etree.SubElement(self.polygon_layer, inkex.addNS('path', 'svg'), attrs)
            if self.polyseg_layer is not None:
                self.poly_segments.add_polygon(vertices)

    def plot_segment(self, x1, y1, x2, y2):
        segment = Line(P(x1, y1), P(x2, y2))
        p1 = (segment.p1 * self.scale) + self.offset
        p2 = (segment.p2 * self.scale) + self.offset
        if (not self.clip_segments_to_doc or (self.doc_clip_rect.point_inside(p1) and self.doc_clip_rect.point_inside(p2))) \
        and  (not self.clip_segments_to_margins or (self.margin_clip_rect.point_inside(p1) and self.margin_clip_rect.point_inside(p2))):
            if self.segment_layer is not None:
                attrs = {'d': 'M %5f %5f L %5f %5f' % (p1.x, p1.y, p2.x, p2.y)}
                attrs['style'] = 'fill:none;stroke:%s;stroke-width:%s;' % (self.segment_stroke, self.segment_stroke_width)
                inkex.etree.SubElement(self.segment_layer, inkex.addNS('path', 'svg'), attrs)
            if self.segpath_layer1 is not None:
                #self.segment_chains.chain_segment(segment)
                self.path_segments.add(segment)

    def rgb2css(self, rgb):
        '''convert the rgb float tuple to a CSS compatible color string'''
        return '#%02x%02x%02x' % (int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))
        

option_info = (
    ('--active-tab', '', 'store', 'string', 'active_tab', '', ''),
    
    ('--scale', '-s', 'store', 'float', 'scale', '1.0', 'Output scale.'),
    ('--fillon', '-f', 'store', 'inkbool', 'fillon', False, 'Fill polygons.'),
    ('--colorfill', '', 'store', 'inkbool', 'colorfill', False, 'Use color fill.'),
    ('--zfill', '-z', 'store', 'inkbool', 'zfill', False, 'Fill color is according to polygon type.'),
    ('--flip', '-F', 'store', 'inkbool', 'flip', False, 'Flip by 90 deg.'),
    ('--symmetry', '-S', 'store', 'int', 'symmetry', '5', 'Degrees of symmetry.'),
    ('--numlines', '-n', 'store', 'int', 'numlines', '30', 'Number of lines.'),
    ('--mid-skinny', '-M', 'store', 'string', 'mid_skinny', '1', 'Midpoint type for skinny diamonds.'),
    ('--mid-fat', '-N', 'store', 'string', 'mid_fat', '1', 'Midpoint type for fat diamonds.'),

    ('--draw-polygons', '', 'store', 'inkbool', 'draw_polygons', True, 'Draw polygons.'),
    ('--polygon-segments', '', 'store', 'inkbool', 'polygon_segments', True, 'Create segment-only polygon layer.'),
    ('--polygon-stroke', '', 'store', 'string', 'polygon_stroke', '#0000ff', 'Polygon CSS stroke color.'),
    ('--polygon-stroke-width', '', 'store', 'string', 'polygon_stroke_width', '.5pt', 'Polygon CSS stroke width.'),

    ('--draw-segments', '', 'store', 'inkbool', 'draw_segments', True, 'Draw segments.'),
    ('--segment-stroke', '', 'store', 'string', 'segment_stroke', '#0000ff', 'Segment CSS stroke color.'),
    ('--segment-stroke-width', '', 'store', 'string', 'segment_stroke_width', '.5pt', 'Segment CSS stroke width.'),

    ('--clip-poly', '-C', 'store', 'inkbool', 'clip_poly', True, 'Clip polygons to document margins.'),
    ('--clip-segments', '', 'store', 'inkbool', 'clip_segments', True, 'Clip segments to document margins.'),
    ('--make-paths', '', 'store', 'inkbool', 'make_paths', True, 'Chain segments into paths.'),
    ('--min-segments', '-m', 'store', 'int', 'minsegments', '1', 'Min segments in path.'),
    ('--closed-paths', '-P', 'store', 'inkbool', 'closed_paths', True, 'Draw only closed paths.'),
    ('--open-path-stroke', '', 'store', 'string', 'open_path_stroke', '#0000ff', 'Open path CSS stroke color.'),
    ('--open-path-stroke-width', '', 'store', 'string', 'open_path_stroke_width', '1pt', 'Open path CSS stroke width.'),
    ('--closed-path-fill', '', 'store', 'string', 'closed_path_fill', '#ffffab', 'Closed path CSS fill color.'),
    ('--closed-path-stroke', '', 'store', 'string', 'closed_path_stroke', '#ff0000', 'Closed path CSS stroke color.'),
    ('--closed-path-stroke-width', '', 'store', 'string', 'closed_path_stroke_width', '1pt', 'Closed path CSS stroke width.'),
    
    ('--margin-left', '', 'store', 'float', 'margin_left', '1', 'Left margin'),
    ('--margin-right', '', 'store', 'float', 'margin_right', '1', 'Right margin'),
    ('--margin-top', '', 'store', 'float', 'margin_top', '1', 'Top margin'),
    ('--margin-bottom', '', 'store', 'float', 'margin_bottom', '1', 'Bottom margin'),

    ('--create-log', '-D', 'store', 'inkbool', 'log_create_log', True, 'Create log files'),
    ('--log-level', '-L', 'store', 'string', 'log_level', 'DEBUG', 'Log level'),
    ('--log-filename', '-O', 'store', 'string', 'log_filename', 'quasink.log', 'Full pathname of log file'),
)

quasink = Quasink(option_info=option_info)
quasink.affect()
