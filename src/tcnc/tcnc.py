"""An Inkscape extension that will output G-code from
selected paths. The G-code is suited to a CNC cutting machine
that has a tangent tool (ie a knife or a brush).

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
import os
import logging
import math
import gettext

import inkex
import cubicsuperpath
import simpletransform

from lib import supereffect
from lib import geom
from lib import paths
from lib import gcode
from lib import simplecam
from lib import paintcam
from lib import svg

_ = gettext.gettext

VERSION = "0.1"

DEFAULT_LOG_LEVEL = logging.DEBUG
logger = logging.getLogger(__name__)

PREVIEW_LAYER_NAME = 'tcnc-preview'
DEBUG_LAYER_NAME = 'tcnc-debug'
# Layers to ignore during input processing
IGNORED_LAYERS=(PREVIEW_LAYER_NAME, DEBUG_LAYER_NAME)

class SVGPreviewPlotter(gcode.PreviewPlotter):
    """Provides a graphical preview of the G-code output.
    Outputs SVG."""
    feed_line_width = 3
    feed_line_opacity = 0.5
    toolmark_offset = 10
    toolmark_width = 30
    toolmark_line_interval = 10
    toolmark_angle_interval = math.pi / 30
    
    # Current XYZA location
    _current_xy = geom.P(0.0, 0.0)
    _current_z = 0.0
    _current_a = 0.0
           
    _styles = {
        'cutpath_end_marker': 'fill-rule:evenodd;fill:#c03030;stroke:#903030;stroke-width:1px;marker-start:none',
        'movepath_end_marker': 'fill-rule:evenodd;fill:#30ff30;stroke:#30c030;stroke-width:1px;marker-start:none',
        'feedline': 'fill:none;stroke:#c03030;stroke-width:%dpx;stroke-opacity:%f%s',
        'feedarc': 'fill:none;stroke:#c03030;stroke-width:%dpx;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:4;stroke-opacity:%f%s',
#        'moveline': 'fill:none;stroke:#30c030;stroke-width:3px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;stroke-miterlimit:4;stroke-dasharray:1.2, 5;stroke-dashoffset:0;marker-end:url(#PreviewLineEnd2)',
        'moveline': 'fill:none;stroke:#30ff30;stroke-width:3px;marker-end:url(#PreviewLineEnd2)',
        'toolmark': 'fill:none;stroke:#cc60c0;stroke-width:1px;stroke-opacity:.75',
    }
    _line_end_markers = (
        ('PreviewLineEnd0', 'M 5.77,0.0 L -2.88,5.0 L -2.88,-5.0 L 5.77,0.0 z',
         'cutpath_end_marker', 'scale(0.4) translate(-4.5,0)'),
        ('PreviewLineEnd1', 'M 5.77,0.0 L -2.88,5.0 L -2.88,-5.0 L 5.77,0.0 z',
         'cutpath_end_marker', 'scale(-0.4) translate(-4.5,0)'),
        ('PreviewLineEnd2', 'M 5.77,0.0 L -2.88,5.0 L -2.88,-5.0 L 5.77,0.0 z',
         'movepath_end_marker', 'scale(0.5) translate(-4.5,0)'),
    )

    def __init__(self, docroot, layer=None, *args, **kwargs):
        super(SVGPreviewPlotter, self).__init__(*args, **kwargs)
        self.svg = svg.SVGContext(docroot, layer)
        # Create Inkscape line end marker glyphs and insert them into the document.
        for marker in self._line_end_markers:
            self.svg.create_simple_marker(marker[0], marker[1],
                                          self._styles[marker[2]], marker[3])

    def plot_move(self, endp):
        """Plot G00 - rapid move to :endp:(x,y,z,a)."""
        p2 = geom.P(endp)
        self.svg.create_line(self._current_xy.x, self._current_xy.y,
                             p2.x, p2.y, self._styles['moveline'])
        self._update_location(endp)
    
    def plot_feed(self, endp):
        """Plot G01 - linear feed to :endp:(x,y,z,a)."""
        p2 = geom.P(endp)
        style = self._build_style('feedline', 'PreviewLineEnd0',
                                  self.feed_line_width, depth=endp[2])
        self.svg.create_line(self._current_xy.x, self._current_xy.y,
                             p2.x, p2.y, style)
        self._draw_tool_marks(geom.Line(self._current_xy, p2),
                              start_angle=self._current_a, end_angle=endp[3])
        self._update_location(endp)
    
    def plot_arc(self, endp, clockwise, arc_x, arc_y):
        """Plot G02/G03 - arc feed to :endp:(x,y,z,a)."""
        p2 = geom.P(endp)
        center = geom.P(arc_x, arc_y) + self._current_xy
        angle = center.angle2(self._current_xy, p2)
        radius = center.distance(self._current_xy)
        sweep_flag = 0 if clockwise else 1
        style = self._build_style('feedarc', 'PreviewLineEnd' + str(sweep_flag),
                                  self.feed_line_width, depth=endp[2])
        attrs = { 'd': 'M %5f %5f A %5f %5f 0 0 %d %5f %5f' % \
                  (self._current_xy.x, self._current_xy.y,
                   radius, radius, sweep_flag, p2.x, p2.y),
                 'style': style }
        self.svg.create_path(attrs)
        self._draw_tool_marks(geom.Arc(self._current_xy, p2, radius,
                                       angle, center),
                              start_angle=self._current_a, end_angle=endp[3])
        self._update_location(endp)
        
    def _draw_tool_marks(self, segment, start_angle, end_angle):
        seglen = segment.length()
        rotation = end_angle - start_angle
        if seglen > 0:
            num_markers = max(1, int(seglen / max(self.toolmark_line_interval,
                                                  self.toolmark_offset)))
        else:
            num_markers = max(1, int(rotation / self.toolmark_angle_interval))
        angle_incr = rotation / num_markers
        point_incr = 1.0 / num_markers
        angle = start_angle
        u = 0
        for unused in range(num_markers):
            p = segment.point_at(u)
            px = p + geom.P.from_polar(self.toolmark_offset, angle - math.pi)
            r = self.toolmark_width / 2
            p1 = px + geom.P.from_polar(r, angle + math.pi/2)
            p2 = px + geom.P.from_polar(r, angle - math.pi/2)
            self.svg.create_line(p.x, p.y, px.x, px.y,
                                   self._styles['toolmark'])
            self.svg.create_line(p1.x, p1.y, p2.x, p2.y,
                                   self._styles['toolmark'])
            angle += angle_incr
            u += point_incr

    def _build_style(self, style_id, marker_id, width, depth):
        if width > 5:
            suffix = ''
            opacity = self.feed_line_opacity
        else:
            suffix = ';marker-end:url(#%s)' % marker_id
            opacity = 1.0
            if width < 1.0:
                width = 1
        return self._styles[style_id] % (width, opacity, suffix)
            
    def _update_location(self, endp):
        self._current_xy = geom.P(endp)
        self._current_z = endp[2]
        self._current_a = endp[3]
    

class TCnc(supereffect.SuperEffect):
    """Inkscape plugin that converts selected SVG elements into gcode suitable for a
    four axis (XYZA) CNC machine with a tangential tool (ie a knife or a brush) as the A axis.
    """
        
    class SVGPreviewPlotter(simplecam.SimpleCAM.PreviewPlotter):
        """Plotter used by SimpleCAM to generate G code preview output."""
        _OPACITY = 0.5
        styles = {'cutpath_end_marker': 'fill-rule:evenodd;fill:#FFFFFF;stroke:#ff0000;stroke-width:1.0pt;marker-start:none',
                  'movepath_end_marker': 'fill-rule:evenodd;fill:#00ff00;stroke:#00ff00;stroke-width:1.0pt;marker-start:none',
                  'feedline': 'fill:none;stroke:#6060c0;stroke-width:%dpx;stroke-opacity:%f%s',
                  'brushline': 'fill:none;stroke:#cc60c0;stroke-width:1px;stroke-opacity:.75',
                  'feedrotate': 'fill:#6060c0;stroke:none;stroke-width:%dpx;fill-opacity:%f',
                  'feedarc': 'fill:none;stroke:#6060c0;stroke-width:%dpx;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:4;stroke-opacity:%f%s',
                  'moveline': 'fill:none;stroke:#00ff00;stroke-width:5px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;stroke-miterlimit:4;stroke-dasharray:1.2,5;stroke-dashoffset:0;marker-end:url(#PreviewLineEnd2)',
                  'anglemark': 'fill:none;stroke:#c060c0;stroke-width:5px;stroke-opacity:0.75',
                  'intervalmark': 'fill:red;stroke:red',}
        
        min_brush_marker_interval = 2
        max_brush_marker_interval = 100
        nominal_brush_marker_count = 10
        
        def __init__(self, inkex, feed_line_width, brush_trail):
            self.inkex = inkex
            self.feed_line_width = feed_line_width
            self.brush_trail = brush_trail
            # Create Inkscape line end marker glyphs and insert them into the document.
            self.inkex.create_simple_marker('PreviewLineEnd0',
                            'M 5.77,0.0 L -2.88,5.0 L -2.88,-5.0 L 5.77,0.0 z',
                            self.styles['cutpath_end_marker'], 'scale(0.4) translate(-4.5,0)')
            self.inkex.create_simple_marker('PreviewLineEnd1',
                            'M 5.77,0.0 L -2.88,5.0 L -2.88,-5.0 L 5.77,0.0 z',
                            self.styles['cutpath_end_marker'], 'scale(-0.4) translate(-4.5,0)')
            self.inkex.create_simple_marker('PreviewLineEnd2',
                            'M 5.77,0.0 L -2.88,5.0 L -2.88,-5.0 L 5.77,0.0 z',
                            self.styles['movepath_end_marker'], 'scale(0.5) translate(-4.5,0)')
        
        def _build_style(self, style_id, marker_end_url, width, depth):
            if width > 5:
                suffix = ''
                opacity = self._OPACITY
            else:
                suffix = ';marker-end:url(%s)' % marker_end_url
                opacity = 1.0
            if width < 1.0:
                width = 1.0
            return self.styles[style_id] % (width, opacity, suffix)
            
        def draw_move_line(self, p1, p2):
            """Draw a line from :p1: to :p2: on the preview layer
            that represents a rapid move.
            """
            self.inkex.create_line(p1.x, p1.y, p2.x, p2.y,
                                   self.styles['moveline'])
        
        def draw_feed_line(self, p1, p2, start_angle, end_angle=None, depth=0.0):
            """Draw a line from :p1: to :p2: on the preview layer
            that represents a linear feed.
            :depth: Current tool depth in machine units
            """
            style = self._build_style('feedline', '#PreviewLineEnd0', self.feed_line_width, depth)
            self.inkex.create_line(p1.x, p1.y, p2.x, p2.y, style)
            if end_angle is not None and not geom.float_eq(end_angle, (p2 -p1).angle()):
                pass
            if end_angle is not None:
                line = geom.Line(p1, p2)
                self._draw_brush_markers(line, start_angle, end_angle)
        
        def draw_feed_rotate(self, center, angle1, angle2, depth=0.0):
            """Draw a tool rotation at :center: from :angle1: to :angle2:
            on the preview layer.
            :depth: Current tool depth in machine units
            This draws a more or less hourglass shape at the point of rotation.
            """
            if self.feed_line_width > 3:
                r = self.feed_line_width / 2
                a90 = math.pi / 2
                p1 = center + geom.P.from_polar(r, angle1 + a90)
                p2 = center + geom.P.from_polar(r, angle2 + a90)
                p3 = center + geom.P.from_polar(r, angle2 - a90)
                p4 = center + geom.P.from_polar(r, angle1 - a90)
                arc1 = geom.Arc(p1, p2, r, angle2 - angle1, center)
                arc2 = geom.Arc(p3, p4, r, -(angle2 - angle1), center)
                style = self.styles['feedrotate'] % (self.feed_line_width, self._OPACITY)
                attrs = { 'd': 'M %5f %5f L %5f %5f A %5f %5f 0 0 %d %5f %5f L %5f %5f L %5f %5f A %5f %5f 0 0 %d %5f %5f L %5f %5f' % \
                          (center.x, center.y, arc1.p1.x, arc1.p1.y,
                           arc1.radius, arc1.radius,
                           0 if arc1.angle < 0 else 1,
                           arc1.p2.x, arc1.p2.y, center.x, center.y,
                           arc2.p1.x, arc2.p1.y,
                           arc1.radius, arc1.radius,
                           0 if arc2.angle < 0 else 1,
                           arc2.p2.x, arc2.p2.y,
                           center.x, center.y,
                           ),
                         'style': style }
                self.inkex.create_path(attrs)
            self._draw_brush_markers(geom.Line(center, center), angle1, angle2)
        
        def draw_feed_arc(self, arc, start_angle, end_angle=None, depth=0.0):
            """Draw an arc on the preview layer that represents a circular feed.
            :arc: A geom.Arc object
            :depth: Current tool depth in machine units
            """
            sweep_flag = 0 if arc.angle < 0 else 1
            style = self._build_style('feedarc',
                                      '#PreviewLineEnd' + str(sweep_flag),
                                      self.feed_line_width, depth)
            attrs = { 'd': 'M %5f %5f A %5f %5f 0 0 %d %5f %5f' % \
                      (arc.p1.x, arc.p1.y, arc.radius, arc.radius,
                       sweep_flag, arc.p2.x, arc.p2.y),
                     'style': style }
            self.inkex.create_path(attrs)
            if end_angle is not None:
                self._draw_brush_markers(arc, start_angle, end_angle)
        
        def draw_angle_marker(self, p, angle, depth=0.0):
            """Draw a mark on the preview layer
            that represents a the current tool angle.
            :depth: Current tool depth in machine units
            """
            r = self.feed_line_width / 2
            p1 = p + geom.P.from_polar(r, angle + math.pi/2)
            p2 = p + geom.P.from_polar(r, angle - math.pi/2)
            self.inkex.create_line(p1.x, p1.y, p2.x, p2.y,
                                   self.styles['anglemark'])
        
        def draw_interval_marker(self, p, depth=0.0):
            """Draw a marker glyph (a small filled circle for example)
            at the specified location :p: on the preview layer.
            :depth: Current tool depth in machine units
            """
            self.inkex.create_circle(p.x, p.y, '2pt', self.styles['intervalmark'])
            
        def _draw_brush_markers(self, segment, start_angle, end_angle):
            seglen = segment.length()
#            if seglen <= geom.EPSILON:
#                return
            num_markers = self.nominal_brush_marker_count
            if seglen > 0.0 and seglen/num_markers < self.min_brush_marker_interval:
                num_markers = int(seglen / self.min_brush_marker_interval)
            num_markers = max(num_markers, self.nominal_brush_marker_count)
            angle_incr = (end_angle - start_angle) / num_markers
            point_incr = 1.0 / num_markers
            angle = start_angle
            u = 0
            for unused in range(num_markers + 1):
                p = segment.point_at(u)
                px = p + geom.P.from_polar(self.brush_trail, angle - math.pi)
                r = self.feed_line_width / 2
                p1 = px + geom.P.from_polar(r, angle + math.pi/2)
                p2 = px + geom.P.from_polar(r, angle - math.pi/2)
                self.inkex.create_line(p.x, p.y, px.x, px.y,
                                       self.styles['brushline'])
                self.inkex.create_line(p1.x, p1.y, p2.x, p2.y,
                                       self.styles['brushline'])
                angle += angle_incr
                u += point_incr
    
    def effect(self):
        """Main entry point for Inkscape plugins.
        """
        # Create a log file for debugging
        if self.options.log_create_log and self.options.log_filename:
            log_path = os.path.abspath(os.path.expanduser(self.options.log_filename))
            log_level = getattr(logging, self.options.log_level, DEFAULT_LOG_LEVEL)
            logging.basicConfig(filename=log_path, filemode='w', level=log_level)
        
        self.process_options()
        
        geom.set_epsilon(self.options.epsilon)
        
        # Create a transform to flip the coordinates
        # from cartesian to SVG (flip Y axis from lower left to upper left).
        page_height = self.get_document_size()[1]
        #page_width = float(self.docroot.get('width'))
        flip_transform_attr = 'translate(0, %f) scale(1, -1)' % page_height
        flip_transform = simpletransform.parseTransform(flip_transform_attr)
        
        # Clear out any previous preview and debug layers if any
        self.clear_layer(PREVIEW_LAYER_NAME)
        self.clear_layer(DEBUG_LAYER_NAME)
        
        # Create a new layer that will contain the G code preview
        self.preview_layer = self.create_layer(PREVIEW_LAYER_NAME)
        self.preview_layer.set('transform', flip_transform_attr)
        self.current_layer = self.preview_layer        
        
        if self.options.debug_layer:
            debug_layer = self.create_layer(DEBUG_LAYER_NAME)
            debug_layer.set('transform', flip_transform_attr)
            # setup geom module for debug output
            geom.DEBUG_EFFECT = self
            geom.DEBUG_LAYER = debug_layer        

        # Get selected SVG elements if any
        rootnodes = self.selected.values()
        if rootnodes:
            # Get the parent layer transform for the
            # selected elements in order to get them in the right place
            shapelist = []
            for node in rootnodes:
                parent_layer = self.get_parent_layer(node)
                matrix = geom.IdentityTransform2D
                if parent_layer is not None:
                    parent_transform = parent_layer.get('transform')
                    matrix = simpletransform.parseTransform(parent_transform)
                shapelist.extend(supereffect.flatten_nodetree((node,), mtransform=matrix,
                                                      ignored_layers=IGNORED_LAYERS))
        else:
            # Use entire document if nothing is selected
            rootnodes = self.document.getroot()
            shapelist = supereffect.flatten_nodetree(rootnodes,
                                             ignored_layers=IGNORED_LAYERS)
        
        # Process the SVG shape/path elements
        cutpath_list = self.process_shapes(shapelist, flip_transform)
        
#        self.trystuff(cutpath_list)
#        return
        preview_plotter=SVGPreviewPlotter(self.document.getroot(),
                                          self.preview_layer)
        cam = self.generate_gcode(cutpath_list, preview_plotter)
#        try:
#            # Generate and export G code
#            cam = self.generate_gcode(cutpath_list)
#        except Exception, e:
#            inkex.errormsg(str(e))
#            return
        try:
            cam.export(self.options.filename, self.options.directory,
                      append_suffix=self.options.append_suffix)
        except IOError, e:
            inkex.errormsg(str(e))            
            

#    def trystuff(self, cutpath_list):
#        fillet_radius = 50
#        min_angle = math.pi * .33
#        min_segment_length = .001#fillet_radius * 2
#        brush_trail_down = 20
#        style = 'fill:none;stroke:#ff0000;stroke-width:3px;'
#        style2 = 'fill:none;stroke:#0000ff;stroke-width:1px;'
#        for path_id, cutpath in cutpath_list:
#            line1 = None
#            for line2 in cutpath:
#                if line1 is not None and isinstance(line1, geom.Line):
#                    if isinstance(line2, geom.Line):
#                        lineside = line1.which_side(line2.p2)
#                        offset = fillet_radius * lineside
#                        offset_line1 = line1.offset(offset)
#                        offset_line2 = line2.offset(offset)
#                        center = offset_line1.intersection(offset_line2)
#                        p1 = line1.projection(center, segment=True)
#                        if p1.distance(line1.p2) > line1.length()/2:
#                            p1 = line1.midpoint()
#                        p2 = line2.projection(center, segment=True)
#                        if p2.distance(line1.p2) > line2.length()/2:
#                            p2 = line2.midpoint()
#                        angle = line1.p2.angle2(line1.p1, line2.p2)
#                        logger.debug('angle=%.2f' % math.degrees(angle))
#                        if line1.length() > min_segment_length and line2.length() > min_segment_length:
#                            if abs(angle) > min_angle:
#                                fillet_angle = center.angle2(p1, p2)
#                                logger.debug('fillet_angle=%.2f' % math.degrees(fillet_angle))
#                                self.create_line(center.x, center.y, p1.x, p1.y, style=style2)
#                                self.create_line(center.x, center.y, p2.x, p2.y, style=style2)
#                                self.create_arc(p1.x, p1.y, p2.x, p2.y, fillet_radius, fillet_angle, style=style)
#                            else:
#                                px = line1.offset(brush_trail_down * lineside).p2
#                                arc1 = geom.Arc.from_two_points_and_tangent(p1, line1.p2, px)
#                                arc2 = geom.Arc.from_two_points_and_tangent(p2, line1.p2, px, reverse=True)
#                                self.create_arc(arc1.p1.x, arc1.p1.y, arc1.p2.x, arc1.p2.y, arc1.radius, arc1.angle, style=style)
#                                self.create_arc(arc2.p1.x, arc2.p1.y, arc2.p2.x, arc2.p2.y, arc2.radius, arc2.angle, style=style)
#                line1 = line2
    
    def process_options(self):
        """Convert option units, etc..."""
        if self.options.gcode_units == 'doc':
            units = self.get_document_units()
            # Determine if the units are metric or imperial.
            # Pica and pixel units are considered imperial for now...
            if units in ('in', 'ft', 'yd', 'pc', 'pt', 'px'):
                self.options.gcode_units = 'in'
            elif units in ('mm', 'cm', 'm', 'km'):
                self.options.gcode_units = 'mm'
            else:
                inkex.errormsg(_('Document units must be inches or mm.'))
                raise Exception()
        # Perform any necessary unit conversion on plugin options
        self.convert_option_units(default_unit_scale=self.get_unit_scale())
        
    def process_shapes(self, shapelist, transform):
        """Convert the SVG shape elements to Lines and CubicBeziers,
        and apply object/layer transforms.
        Returns a list of tuples ((path-id, path), ...).
        """
        cutpath_list = []
        for node, layer_transform in shapelist:
            # Convert the shape element to a simplepath
            path = supereffect.convert_element_to_path(node)
            
            # Convert the simplepath to a 'cubicsuperpath' which is
            # just a list of cubic bezier curves.
            # This seems to be how most Inkscape plugins do things...
            # TODO: This is really an unnecessary step since
            # we could deal with SVG shapes directly...
            csp = cubicsuperpath.CubicSuperPath(path)
            
            # Apply the SVG element transform and it's layer transform to the
            # path segments so that we are working in absolute coordinates.
            # Transform SVG coordinates into cartesian (ie G code) coordinates
            # (flip the Y axis from upper left to lower left).
#            node_transform = simpletransform.parseTransform(node.get('transform'))
#            node_transform = simpletransform.composeTransform(node_transform, transform)
#            node_transform = simpletransform.composeTransform(node_transform, layer_transform)
            node_transform = simpletransform.composeTransform(transform, layer_transform)
            simpletransform.applyTransformToPath(node_transform, csp)
                        
            # Convert cubic path segments to curves and lines
            # TODO: get bounding boxes and calculate offset if doc origin not used.
            cutpath = paths.Path(name=node.get('id'))
            for subcsp in csp:
                for i in range(1,len(subcsp)):
                    p1 = geom.P(subcsp[i-1][1])
                    c1 = geom.P(subcsp[i-1][2])
                    p2 = geom.P(subcsp[i][0])
                    c2 = geom.P(subcsp[i][1])
                    if p1 == c1 and p2 == c2:
                        segment = geom.Line(p1, p2)
                    else:
                        segment = geom.CubicBezier(p1, c1, p2, c2)
                    cutpath.append(segment)

            cutpath_list.append(cutpath)
        
        return cutpath_list
            
    def generate_gcode(self, cutpath_list, preview_plotter):
        """Generate G code from cutpaths."""
        preview_plotter.toolmark_width = self.options.brush_size
        preview_plotter.toolmark_offset = self.options.brush_trail_down
        gc = gcode.GCode(zsafe=float(self.options.z_safe),
                         zfeed=float(self.options.z_feed),
                         xyfeed=float(self.options.xy_feed),
                         afeed=float(self.options.a_feed * 60),
                         plotter=preview_plotter)
        gc.set_axis_offsets('XYZA', (self.options.x_offset,
                            self.options.y_offset, self.options.z_offset,
                            self.options.a_offset + self.options.brush_angle))
        gc.set_axis_scales('XYZA', (self.options.x_scale,
                           self.options.y_scale, self.options.z_scale,
                           self.options.a_scale))
        gc.units = self.options.gcode_units
        gc.unit_scale = 1.0 / self.get_unit_scale(units=gc.units)
        gc.trajectory_mode = self.options.traj_mode
        gc.trajectory_tolerance = self.options.traj_tolerance
        
#        cam = paintcam.PaintCAM(gc, preview_plotter=TCnc.SVGPreviewPlotter(self, self.options.brush_size, self.options.brush_trail_down))
        cam = paintcam.PaintCAM(gc, preview_plotter=simplecam.SimpleCAM.PreviewPlotter())
        cam.optimize_rapid_moves = self.options.sort_paths
        cam.tool_width = self.options.brush_size
        cam.biarc_tolerance = self.options.biarc_tolerance
        cam.biarc_max_depth = self.options.biarc_max_depth
        cam.biarc_line_flatness = self.options.line_flatness
        cam.brush_landing_angle = self.options.brush_landing_angle
        cam.brush_landing_end_height = self.options.brush_landing_end_height
        cam.brush_landing_start_height = self.options.brush_landing_start_height
        cam.brush_liftoff_angle = self.options.brush_liftoff_angle
        cam.brush_liftoff_height = self.options.brush_liftoff_height
        cam.brush_overshoot = self.options.brush_overshoot
        cam.brush_reload_before_path = self.options.brush_reload_path and self.options.brush_dwell > 0.0
        cam.brush_reload_after_interval = self.options.brushstroke_max > 0.0
        cam.brush_reload_angle = self.options.brush_reload_angle
        cam.brush_reload_dwell = self.options.brush_dwell
        cam.brush_trail_down = self.options.brush_trail_down
        cam.brush_depth = self.options.brush_depth
        cam.brush_overshoot = cam.brush_trail_down + cam.tool_width/2
        cam.use_fillets = not self.options.brush_no_fillets
        cam.close_polygons = self.options.brush_close_polygons
        #cam.fillet_radius = cam.tool_width / 2
        if self.options.brushstroke_max > 0.0:
            cam.feed_interval = self.options.brushstroke_max
            cam.brush_flip_before_reload = True
        
        cam.generate_gcode(cutpath_list, cam.brush_depth)
        return cam
        
option_info = [
    supereffect.optargs('--active-tab', type='string', dest='active_tab'),
    
    supereffect.optargs('--origin-ref', dest='origin_ref', default='paper', help='Lower left origin reference.'),
    supereffect.optargs('--sort-paths', type='inkbool', dest='sort_paths', default=True, help='Sort paths to minimize rapid moves.'),
    supereffect.optargs('--biarc-tolerance', type='float', convert_to='world', default=0.01, help='Biarc approximation fitting tolerance.'),                
    supereffect.optargs('--biarc-max-depth', type='int', default=4, help='Biarc approximation maximum curve splitting recursion depth.'),                
    supereffect.optargs('--line-flatness', type='float', convert_to='world', default=0.001, help='Curve to line flatness.'),                
    supereffect.optargs('--min-arc-radius', type='float', convert_to='world', default=0.01, help='All arcs having radius less than minimum will be considered as straight line.'),
    supereffect.optargs('--epsilon', type='float', convert_to='world', default=0.0001, help='Epsilon'),

    supereffect.optargs('--z-scale', type='float', dest='z_scale', default=1.0, help='Scale factor Z'), 
    supereffect.optargs('--z-offset', type='float', convert_to='world', default=0.0, help='Offset along Z'),
    supereffect.optargs('--x-scale', type='float', dest='x_scale', default=1.0, help='Scale factor X'), 
    supereffect.optargs('--x-offset', type='float', convert_to='world', default=0.0, help='Offset along X'),
    supereffect.optargs('--y-scale', type='float', dest='y_scale', default=1.0, help='Scale factor Y'), 
    supereffect.optargs('--y-offset', type='float', convert_to='world', default=0.0, help='Offset along Y'),
    supereffect.optargs('--a-scale', type='float', dest='a_scale', default=1.0, help='Angular scale along rotational axis'),
    supereffect.optargs('--a-offset', type='float', convert_to='world', default=0.0, help='Angular offset along rotational axis'),
   
    supereffect.optargs('--gcode-units', dest='gcode_units', default='in', help='G code output units (inch or mm).'),
    supereffect.optargs('--xy-feed', type='float', dest='xy_feed', default=10.0, help='XY axis feed rate in unit/s'),
    supereffect.optargs('--z-feed', type='float', dest='z_feed', default=10.0, help='Z axis feed rate in unit/s'),
    supereffect.optargs('--a-feed', type='float', dest='a_feed', default=60.0, help='A axis feed rate in deg/s'),
    supereffect.optargs('--z-safe', type='float', convert_to='world', default=5.0, help='Z axis safe height for rapid moves'),
    supereffect.optargs('--z-wait', type='float', dest='z_wait', default=500, help='Z axis wait (milliseconds)'),
    supereffect.optargs('--traj-mode', dest='traj_mode', default='', help='Trajectory planning mode.'),
    supereffect.optargs('--traj-tolerance', type='float', dest='traj_tolerance', default='0', help='Trajectory blending tolerance.'),
    
    supereffect.optargs('--brush-angle', type='float', convert_to='rad', default=90, help='Brush angle'),
    supereffect.optargs('--brush-overshoot', type='float', dest='brush_overshoot', convert_to='world', default=0.5, help='Brushstroke overshoot distance.'),
    supereffect.optargs('--brush-liftoff-height', type='float', convert_to='world', default=0.1, help='Brushstroke liftoff height.'),
    supereffect.optargs('--brush-liftoff-angle', type='float', convert_to='rad', default=45, help='Brushstroke liftoff angle.'),
    supereffect.optargs('--brush-landing-start-height', type='float', convert_to='world', default=0.1, help='Brushstroke landing start height.'),
    supereffect.optargs('--brush-landing-end-height', type='float', convert_to='world', default=-0.2, help='Brushstroke landing end height.'),
    supereffect.optargs('--brush-landing-angle', type='float', convert_to='rad', default=45, help='Brushstroke landing angle.'),
    supereffect.optargs('--brush-flip-stroke', type='inkbool', dest='brush_flip_stroke', default=True, help='Flip brush before every stroke.'),
    supereffect.optargs('--brush-flip-path', type='inkbool', dest='brush_flip_path', default=True, help='Flip after path.'),
    supereffect.optargs('--brush-flip-reload', type='inkbool', dest='brush_flip_reload', default=True, help='Flip after reload.'),
    supereffect.optargs('--brush-reload', type='inkbool', dest='brush_reload', default=True, help='Enable brush reload.'),
    supereffect.optargs('--brush-reload-path', type='inkbool', dest='brush_reload_path', default=True, help='Reload brush after every path.'),
    supereffect.optargs('--brushstroke-max', type='float', convert_to='world', default=0.0, help='Maximum brushstroke distance.'),
    supereffect.optargs('--brushstroke-overlap', type='float', convert_to='world', default=0.0, help='Brushstroke overlap.'),
    supereffect.optargs('--brush-dwell', type='float', dest='brush_dwell', default=0.0, help='Brush reload time (seconds).'),
    supereffect.optargs('--brush-reload-angle', type='float', convert_to='rad', default=90.0, help='Brush reload angle (degrees).'),
    supereffect.optargs('--brush-size', type='float', convert_to='world', default=1.0, help='Brush size'),
    supereffect.optargs('--brush-trail-down', type='float', convert_to='world', default=0.25, help='Brush trail down'),
    supereffect.optargs('--brush-trail-up', type='float', convert_to='world', default=0.15, help='Brush trail up'),
    supereffect.optargs('--brush-depth', type='float', convert_to='world', default=-0.15, help='Brush depth'),
    supereffect.optargs('--brush-no-fillets', type='inkbool', default=False, help='No fillets at sharp corners'),
    supereffect.optargs('--brush-fillet-radius', type='float', convert_to='world', default=0.0, help='Brush trail'),
    supereffect.optargs('--brush-fillet-radius-auto', type='inkbool', default=False, help='Auto-calc fillet radius'),
    supereffect.optargs('--brush-close-polygons', type='inkbool', default=False, help='Close polygons'),
   
    supereffect.optargs('--directory', type='string', dest='directory', default='~', help='Directory for gcode file'),
    supereffect.optargs('--filename', type='string', dest='filename', default='untitled', help='File name'), 
    supereffect.optargs('--append-suffix', type='inkbool', dest='append_suffix', default=True, help='Append auto-incremented numeric suffix to filename'), 
    supereffect.optargs('--separate-layers', type='inkbool', dest='separate_layers', default=True, help='Seaparate gcode file per layer'), 
    supereffect.optargs('--create-log', type='inkbool', dest='log_create_log', default=True, help='Create log files'),
    supereffect.optargs('--log-level', type='string', dest='log_level', default='DEBUG', help='Log level'),
    supereffect.optargs('--log-filename', type='string', dest='log_filename', default='tcnc.log', help='Full pathname of log file'),
   
    supereffect.optargs('--preview-show', type='inkbool', dest='preview_show', default=True, help='Show generated cut paths on preview layer.'),
    supereffect.optargs('--debug-layer', type='inkbool', dest='debug_layer', default=True, help='Create debug layer.'),
    supereffect.optargs('--debug-biarcs', type='inkbool', dest='debug_biarcs', default=True),
   
    # Unused lagacy options - derived from gcodetools
    # TODO: use or delete
    supereffect.optargs('--z-depth', type='float', convert_to='world', default=-0.125, help='Z full depth of cut'),
    supereffect.optargs('--z-step', type='float', convert_to='world', default=-0.125, help='Z cutting step depth'),
    supereffect.optargs('--path-to-gcode-order',type='string', dest='path_to_gcode_order', default='path by path', help='Defines cutting order path by path or layer by layer.'), 
    supereffect.optargs('--path-to-gcode-depth-function', type='string', dest='path_to_gcode_depth_function', default='zd', help='Path to gcode depth function.'),
]

if __name__ == '__main__':
    tcnc = TCnc(option_info)
    tcnc.affect()
