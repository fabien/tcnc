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
from lib.geom import P

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
        logger.debug('polygon-stroke: %s' % self.options.polygon_stroke)
        logger.debug('segment-stroke: %s' % self.options.segment_stroke)

        self.scale = self.options.scale * self.SCALE_SCALE
        self.offset = P(self.view_center)
        
        doc_size = P(self.get_document_size())
        
        plotter = SVGQuasiPlotter()
        if self.options.draw_polygons:
            plotter.polygon_layer = self.create_layer('quasink_polygons')
        if self.options.draw_segments:
            plotter.segment_layer = self.create_layer('quasink_segments')
        if self.options.make_paths:
            plotter.segpath_layer = self.create_layer('quasink_segpaths')
        plotter.scale = self.scale
        plotter.offsetx = self.offset[0]
        plotter.offsety = self.offset[1]
        plotter.clip_rect = geom.Rectangle(P(0.0, 0.0), doc_size)
        plotter.is_clipped = bool(self.options.clip_poly)
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
        
        if self.options.make_paths:
            style_open_path = simplestyle.formatStyle({ 'fill': 'none',
                                 'stroke': simplestyle.formatColoria(simplestyle.parseColor(self.options.open_path_stroke)),
                                 'stroke-width': self.options.open_path_stroke_width,
                                 'stroke-linecap': 'butt',
                                 'stroke-linejoin': 'round'
                                 })
            style_closed_path = simplestyle.formatStyle({ 'fill': simplestyle.formatColoria(simplestyle.parseColor(self.options.closed_path_fill)),
                                 'stroke': simplestyle.formatColoria(simplestyle.parseColor(self.options.closed_path_stroke)),
                                 'stroke-width': self.options.closed_path_stroke_width,
                                 'stroke-linejoin': 'round'
                                 })
            segpath_list = plotter.pathtable.get_paths()
            # Sort segment paths so that the largest are at the bottom of the Z-order
            segpath_list.sort(key=SegmentPath.bounding_box_area, reverse=True)
            for segpath in segpath_list:
                segchain = segpath.segment_list
                if len(segchain) >= self.options.minsegments:
                    if segchain[0].p1 == segchain[-1].p2:
                        plotter.plot_segpath(segchain, style_closed_path)
                    elif not self.options.closed_paths:
                        plotter.plot_segpath(segchain, style_open_path)

class SegmentPath(object):
    def __init__(self, segments):
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
    chain_startp = {}   # Hashtable of segment start coordinates
    chain_endp = {}     # Hashtable of segment end coordinates
    
    def reverse_chain(self, chain):
        '''Reverse the chain and re-insert into hashtable'''
        chain_rev = []
        for segment in chain:
            chain_rev.insert(0, segment.reverse())
        # And re-insert into hash table
        self.chain_startp[chain_rev[0].p1] = chain_rev
        self.chain_endp[chain_rev[-1].p2] = chain_rev
        
    
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
    
    def get_paths(self):
        #return self.paths
        return map(SegmentPath, self.chain_startp.itervalues())
    
    
class BufQuasiPlotter(quasi.QuasiPlotter):
    '''A plotter that just accumulates the quasi geometry.'''
    polygons = []
    segments = []
    
    def plot_polygon(self, vertices):
        self.polygons.append(vertices)

    def plot_segment(self, x1, y1, x2, y2):
        self.segments.append(geom.Line(P(x1, y1), P(x2, y2)))


class SVGQuasiPlotter(quasi.QuasiPlotter):
    '''SVG output for quasi'''
    scale = 1.0
    offsetx = 0.0
    offsety = 0.0
    polygon_layer = None
    segment_layer = None
    segpath_layer = None
    clip_rect = [P(0, 0), P(1.0, 1.0)]
    is_clipped = False
    pathtable = SegmentTable()
    fill_polygons = False
    polygon_stroke = '#333333'
    segment_stroke = '#666666'
    polygon_stroke_width = '.5pt'
    segment_stroke_width = '.5pt'
    
    def plot_segpath(self, segpath, style=None):
        x = (segpath[0].p1.x * self.scale) + self.offsetx
        y = (segpath[0].p1.y * self.scale) + self.offsety
        #if not self.is_clipped or self.clip_rect.point_inside(P(x, y)):
        d = 'M %f %f L' % (x, y)
        for segment in segpath:
            x = (segment.p2.x * self.scale) + self.offsetx
            y = (segment.p2.y * self.scale) + self.offsety
            d = d + ' %f,%f' % (x, y)
        attrs = {'d': d}
        if style:
            attrs['style'] = style
        else:
            attrs['style'] = Quasink.styles['segpath']
        inkex.etree.SubElement(self.segpath_layer, inkex.addNS('path', 'svg'), attrs)
    
    def plot_polygon(self, vertices):
        ''''''
        if self.polygon_layer is None:
            return
        x = (vertices[0][0] * self.scale) + self.offsetx
        y = (vertices[0][1] * self.scale) + self.offsety
        if not self.is_clipped or self.clip_rect.point_inside(P(x, y)):
            d = 'M %f %f L' % (x, y)
            for p in vertices[1:]:
                x = (p[0] * self.scale) + self.offsetx
                y = (p[1] * self.scale) + self.offsety
                if self.is_clipped and not self.clip_rect.point_inside(P(x, y)):
                    return
                d = d + ' %f,%f' % (x, y)
            attrs = {'d': d}
            if self.fill_polygons:
                fill = self.rgb2css(self.get_fill_color_rgb())
            else:
                fill = 'none'
            attrs['style'] = 'fill:%s;stroke:%s;stroke-width:%s;' % (fill, self.polygon_stroke, self.polygon_stroke_width)
            inkex.etree.SubElement(self.polygon_layer, inkex.addNS('path', 'svg'), attrs)

    def plot_segment(self, x1, y1, x2, y2):
        segment = geom.Line(P(x1,y1),P(x2,y2))
        x1 = (x1 * self.scale) + self.offsetx
        y1 = (y1 * self.scale) + self.offsety
        x2 = (x2 * self.scale) + self.offsetx
        y2 = (y2 * self.scale) + self.offsety
        if not self.is_clipped or self.clip_rect.line_inside(geom.Line(P(x1,y1),P(x2,y2))):
            if self.segment_layer is not None:
                attrs = {'d': 'M %5f %5f L %5f %5f' % (x1, y1, x2, y2)}
                attrs['style'] = 'fill:none;stroke:%s;stroke-width:%s;' % (self.segment_stroke, self.segment_stroke_width)
                inkex.etree.SubElement(self.segment_layer, inkex.addNS('path', 'svg'), attrs)
            if self.segpath_layer is not None:
                self.pathtable.chain_segment(segment)

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
    ('--polygon-stroke', '', 'store', 'string', 'polygon_stroke', '#0000ff', 'Polygon CSS stroke color.'),
    ('--polygon-stroke-width', '', 'store', 'string', 'polygon_stroke_width', '.5pt', 'Polygon CSS stroke width.'),

    ('--draw-segments', '', 'store', 'inkbool', 'draw_segments', True, 'Draw segments.'),
    ('--segment-stroke', '', 'store', 'string', 'segment_stroke', '#0000ff', 'Segment CSS stroke color.'),
    ('--segment-stroke-width', '', 'store', 'string', 'segment_stroke_width', '.5pt', 'Segment CSS stroke width.'),

    ('--clip-poly', '-C', 'store', 'inkbool', 'clip_poly', True, 'Clip polygons to document bounds.'),
    ('--make-paths', '', 'store', 'inkbool', 'make_paths', True, 'Chain segments into paths.'),
    ('--min-segments', '-m', 'store', 'int', 'minsegments', '1', 'Min segments in path.'),
    ('--closed-paths', '-P', 'store', 'inkbool', 'closed_paths', True, 'Draw only closed paths.'),
    ('--open-path-stroke', '', 'store', 'string', 'open_path_stroke', '#0000ff', 'Open path CSS stroke color.'),
    ('--open-path-stroke-width', '', 'store', 'string', 'open_path_stroke_width', '1pt', 'Open path CSS stroke width.'),
    ('--closed-path-fill', '', 'store', 'string', 'closed_path_fill', '#ffffab', 'Closed path CSS fill color.'),
    ('--closed-path-stroke', '', 'store', 'string', 'closed_path_stroke', '#ff0000', 'Closed path CSS stroke color.'),
    ('--closed-path-stroke-width', '', 'store', 'string', 'closed_path_stroke_width', '1pt', 'Closed path CSS stroke width.'),
    
    ('--create-log', '-D', 'store', 'inkbool', 'log_create_log', True, 'Create log files'),
    ('--log-level', '-L', 'store', 'string', 'log_level', 'DEBUG', 'Log level'),
    ('--log-filename', '-O', 'store', 'string', 'log_filename', 'quasink.log', 'Full pathname of log file'),
)

quasink = Quasink(option_info=option_info)
quasink.affect()
