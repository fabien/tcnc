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
_ = gettext.gettext

import inkex
import cubicsuperpath
import simpletransform

from lib import svg
from lib import geom
from lib import gcode
from lib import simplecam
from lib import paintcam

VERSION = "0.1"

DEFAULT_LOG_LEVEL = logging.DEBUG
logger = logging.getLogger(__name__)

PREVIEW_LAYER_NAME = 'tcnc-preview'
DEBUG_LAYER_NAME = 'tcnc-debug'


class TCnc(svg.SuperEffect):
    """Inkscape plugin that converts selected SVG elements into gcode suitable for a
    four axis (XYZA) CNC machine with a tangential tool (ie a knife or a brush) as the A axis.
    """
        
    class SVGPreviewPlotter(simplecam.SimpleCAM.PreviewPlotter):
        """Plotter used by SimpleCAM to generate G code preview output."""
        _OPACITY = 0.5
        styles = {'cutpath_end_marker': 'fill-rule:evenodd;fill:#FFFFFF;stroke:#ff0000;stroke-width:1.0pt;marker-start:none',
                  'movepath_end_marker': 'fill-rule:evenodd;fill:#00ff00;stroke:#00ff00;stroke-width:1.0pt;marker-start:none',
                  'feedline': 'fill:none;stroke:#6060c0;stroke-width:%dpx;stroke-opacity:%f%s',
                  'feedrotate': 'fill:#6060c0;stroke:none;stroke-width:%dpx;fill-opacity:%f',
                  'feedarc': 'fill:none;stroke:#6060c0;stroke-width:%dpx;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:4;stroke-opacity:%f%s',
                  'moveline': 'fill:none;stroke:#00ff00;stroke-width:%dpx;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;stroke-miterlimit:4;stroke-dasharray:1.2, 5;stroke-dashoffset:0;marker-end:url(#PreviewLineEnd2)',
                  'intervalmark': 'fill:red;stroke:red',}
        
        def __init__(self, inkex, feed_line_width):
            self.inkex = inkex
            self.feed_line_width = feed_line_width
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
        
        def draw_feed_line(self, p1, p2, width=1, depth=0.0):
            """Draw a line from :p1: to :p2: on the preview layer
            that represents a linear feed.
            :width: Line (tool) width in machine units
            :depth: Current tool depth in machine units
            """
            style = self._build_style('feedline', '#PreviewLineEnd0', width, depth)
            self.inkex.create_line(p1.x, p1.y, p2.x, p2.y, style)
        
        def draw_feed_rotate(self, center, angle1, angle2, width=1.0, depth=0.0):
            """Draw a tool rotation at :center: from :angle1: to :angle2:
            on the preview layer.
            :width: Line (tool) width in machine units
            :depth: Current tool depth in machine units
            This draws a more or less hourglass shape at the point of rotation.
            """
            if width > 3:
                r = self.feed_line_width / 2
                a90 = math.pi / 2
                p1 = center + geom.P.from_polar(r, angle1 + a90)
                p2 = center + geom.P.from_polar(r, angle2 + a90)
                p3 = center + geom.P.from_polar(r, angle2 - a90)
                p4 = center + geom.P.from_polar(r, angle1 - a90)
                arc1 = geom.Arc(p1, p2, r, angle2 - angle1, center)
                arc2 = geom.Arc(p3, p4, r, -(angle2 - angle1), center)
                style = self.styles['feedrotate'] % (width, self._OPACITY)
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
        
        def draw_feed_arc(self, arc, width=1, depth=0.0):
            """Draw an arc on the preview layer that represents a circular feed.
            :arc: A geom.Arc object
            :width: Line (tool) width in machine units
            :depth: Current tool depth in machine units
            """
            sweep_flag = 0 if arc.angle < 0 else 1
            style = self._build_style('feedarc',
                                      '#PreviewLineEnd' + str(sweep_flag),
                                      width, depth)
            attrs = { 'd': 'M %5f %5f A %5f %5f 0 0 %d %5f %5f' % \
                      (arc.p1.x, arc.p1.y, arc.radius, arc.radius,
                       sweep_flag, arc.p2.x, arc.p2.y),
                     'style': style }
            self.inkex.create_path(attrs)
        
        def draw_interval_marker(self, p, depth=0.0):
            """Draw a marker glyph (a small filled circle for example)
            at the specified location :p: on the preview layer.
            :depth: Current tool depth in machine units
            """
            self.create_circle(p.x, p.y, '2pt', self.styles['intervalmark'])
    
    def effect(self):
        """Main entry point for Inkscape plugins.
        """
        # Create a log file for debugging
        if self.options.log_create_log and self.options.log_filename:
            log_path = os.path.abspath(os.path.expanduser(self.options.log_filename))
            log_level = getattr(logging, self.options.log_level, DEFAULT_LOG_LEVEL)
            logging.basicConfig(filename=log_path, filemode='w', level=log_level)
        
        self.process_options()
        if self.options.units not in ('in', 'mm'):
            inkex.errormsg(_('Document units must be either inches or mm.'))
            return
        
        # Get selected SVG elements if any
        rootnodes = self.selected.values()
        if not rootnodes:
            # Use entire document if nothing is selected
            rootnodes = self.document.getroot()
        
        # Clear out the old preview layer if any
        self.clear_layer(PREVIEW_LAYER_NAME)
        self.clear_layer(DEBUG_LAYER_NAME)

        # Extract all the recognized SVG shape elements
        skip_layers=(PREVIEW_LAYER_NAME,)
        shapelist = svg.flatten_nodetree(rootnodes, skip_layers=skip_layers)
        
        # Create a new layer that will contain the G code preview
        self.preview_layer = self.create_layer(PREVIEW_LAYER_NAME)
        
        if self.options.debug_layer:
            debug_layer = self.create_layer(DEBUG_LAYER_NAME)
            # setup geom module for debug output
            geom.DEBUG_EFFECT = self
            geom.DEBUG_LAYER = debug_layer
        
        self.current_layer = self.preview_layer        

        # Add a transform to the preview layer to flip the coordinates
        # from cartesian to SVG (flip Y axis from lower left to upper left).
        page_height = self.get_document_size()[1]
        #page_width = float(self.docroot.get('width'))
        flip_transform_attr = 'translate(0, %f) scale(1, -1)' % page_height
        self.preview_layer.set('transform', flip_transform_attr)
        flip_transform = simpletransform.parseTransform(flip_transform_attr)
        
        # Process the SVG shape/path elements
        cutpath_list = self.process_shapes(shapelist, flip_transform)
        
        # Generate and export G code
        cam = self.generate_gcode(cutpath_list)
        try:
            cam.export(self.options.filename, self.options.directory,
                      append_suffix=self.options.append_suffix)
        except IOError, e:
            inkex.errormsg(str(e))            
            
    def process_options(self):
        """Convert option units, etc..."""
        if self.options.units == 'doc':
            self.options.units = self.get_document_units()
        unit_scale = inkex.uuconv[self.options.units]
        # Perform any necessary unit conversion on plugin options
        self.convert_option_units(default_unit_scale=unit_scale)
        
    def process_shapes(self, shapelist, transform):
        """Process the SVG shape elements and generate G code.
        Return a list of cut paths as SVG path elements.
        """
        cutpath_list = []
        for node, layer_transform in shapelist:
            # Convert the shape element to a simplepath
            path = svg.convert_element_to_path(node)
            
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
                        
            # Convert cubic path segments to arcs and lines
            cutpath = []
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

            # The cutpath is a tuple ('path-id', path_list)
            cutpath_list.append((node.get('id'), cutpath))
        
        return cutpath_list
            
    def generate_gcode(self, cutpath_list):
        """Generate G code from cutpaths."""
        gc = gcode.GCode(zsafe=float(self.options.z_safe),
                         zfeed=float(self.options.z_feed),
                         xyfeed=float(self.options.xy_feed),
                         afeed=float(self.options.a_feed * 60))
        gc.set_axis_offsets('XYZA', (self.options.x_offset,
                            self.options.y_offset, self.options.z_offset,
                            self.options.a_offset + self.options.brush_angle))
        gc.set_axis_scales('XYZA', (self.options.x_scale,
                           self.options.y_scale, self.options.z_scale,
                           self.options.a_scale))
        gc.units = self.options.units
        gc.unit_scale = self.get_unit_scale(self.options.units)
        gc.trajectory_mode = self.options.traj_mode
        gc.trajectory_tolerance = self.options.traj_tolerance
        
        cam = paintcam.PaintCAM(gc, preview_plotter=TCnc.SVGPreviewPlotter(self, self.options.brush_size))
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
        cam.brush_reload_angle = self.options.brush_reload_angle
        cam.brush_reload_dwell = self.options.brush_dwell
        if self.options.brushstroke_max > 0.0:
            cam.feed_interval = self.options.brushstroke_max
            cam.brush_flip_before_reload = True
        
        cam.generate_gcode(cutpath_list)
        return cam
        
option_info = [
    svg.optargs('--active-tab', type='string', dest='active_tab'),
    
    svg.optargs('--units', dest='units', default='in', help='Document units.'),
    svg.optargs('--auto-select-paths', type='inkbool', dest='auto_select_paths', default=True, help='Select all paths if nothing is selected.'),
    svg.optargs('--biarc-tolerance', type='float', convert_to='world', default=0.01, help='Biarc approximation fitting tolerance.'),                
    svg.optargs('--biarc-max-depth', type='int', default=4, help='Biarc approximation maximum curve splitting recursion depth.'),                
    svg.optargs('--line-flatness', type='float', convert_to='world', default=0.001, help='Curve to line flatness.'),                
    svg.optargs('--min-arc-radius', type='float', convert_to='world', default=0.01, help='All arcs having radius less than minimum will be considered as straight line.'),
    svg.optargs('--sort-paths', type='inkbool', dest='sort_paths', default=True, help='Sort paths to minimize rapid moves.'),

    svg.optargs('--origin-ref', dest='origin_ref', default='paper', help='Lower left origin reference.'),
    svg.optargs('--z-scale', type='float', dest='z_scale', default=1.0, help='Scale factor Z'), 
    svg.optargs('--z-offset', type='float', convert_to='world', default=0.0, help='Offset along Z'),
    svg.optargs('--x-scale', type='float', dest='x_scale', default=1.0, help='Scale factor X'), 
    svg.optargs('--x-offset', type='float', convert_to='world', default=0.0, help='Offset along X'),
    svg.optargs('--y-scale', type='float', dest='y_scale', default=1.0, help='Scale factor Y'), 
    svg.optargs('--y-offset', type='float', convert_to='world', default=0.0, help='Offset along Y'),
    svg.optargs('--a-scale', type='float', dest='a_scale', default=1.0, help='Angular scale along rotational axis'),
    svg.optargs('--a-offset', type='float', convert_to='world', default=0.0, help='Angular offset along rotational axis'),
   
    svg.optargs('--xy-feed', type='float', dest='xy_feed', default=10.0, help='XY axis feed rate in unit/s'),
    svg.optargs('--z-feed', type='float', dest='z_feed', default=10.0, help='Z axis feed rate in unit/s'),
    svg.optargs('--a-feed', type='float', dest='a_feed', default=60.0, help='A axis feed rate in deg/s'),
    svg.optargs('--z-safe', type='float', convert_to='world', default=5.0, help='Z axis safe height for rapid moves'),
    svg.optargs('--z-wait', type='float', dest='z_wait', default=500, help='Z axis wait (milliseconds)'),
    svg.optargs('--traj-mode', dest='traj_mode', default='G64', help='Trajectory planning mode.'),
    svg.optargs('--traj-tolerance', type='float', dest='traj_tolerance', default='0', help='Trajectory blending tolerance.'),
    
    svg.optargs('--brush-angle', type='float', convert_to='rad', default=90, help='Brush angle'),
    svg.optargs('--brush-overshoot-enabled', type='inkbool', dest='brush_overshoot_enabled', default=False, help='Enable brush overshoot.'),
    svg.optargs('--brush-overshoot', type='float', dest='brush_overshoot', default=0.5, help='Brushstroke overshoot distance.'),
    svg.optargs('--brush-liftoff-height', type='float', convert_to='world', default=45, help='Brushstroke liftoff height.'),
    svg.optargs('--brush-liftoff-angle', type='float', convert_to='rad', default=45, help='Brushstroke liftoff angle.'),
    svg.optargs('--brush-landing-start-height', type='float', convert_to='world', default=45, help='Brushstroke landing start height.'),
    svg.optargs('--brush-landing-end-height', type='float', convert_to='world', default=45, help='Brushstroke landing end height.'),
    svg.optargs('--brush-landing-angle', type='float', convert_to='rad', default=45, help='Brushstroke landing angle.'),
    svg.optargs('--brush-flip-stroke', type='inkbool', dest='brush_flip_stroke', default=True, help='Flip brush before every stroke.'),
    svg.optargs('--brush-flip-path', type='inkbool', dest='brush_flip_path', default=True, help='Flip after path.'),
    svg.optargs('--brush-flip-reload', type='inkbool', dest='brush_flip_reload', default=True, help='Flip after reload.'),
    svg.optargs('--brush-reload', type='inkbool', dest='brush_reload', default=True, help='Enable brush reload.'),
    svg.optargs('--brush-reload-path', type='inkbool', dest='brush_reload_path', default=True, help='Reload brush after every path.'),
    svg.optargs('--brushstroke-max', type='float', convert_to='world', default=0.0, help='Maximum brushstroke distance.'),
    svg.optargs('--brushstroke-overlap', type='float', convert_to='world', default=0.0, help='Brushstroke overlap.'),
    svg.optargs('--brush-dwell', type='float', dest='brush_dwell', default=0.0, help='Brush reload time (seconds).'),
    svg.optargs('--brush-reload-angle', type='float', convert_to='rad', default=90.0, help='Brush reload angle (degrees).'),
    svg.optargs('--brush-size', type='float', convert_to='world', default=1.0, help='Brush size'),
   
    svg.optargs('--directory', type='string', dest='directory', default='~', help='Directory for gcode file'),
    svg.optargs('--filename', type='string', dest='filename', default='untitled', help='File name'), 
    svg.optargs('--append-suffix', type='inkbool', dest='append_suffix', default=True, help='Append auto-incremented numeric suffix to filename'), 
    svg.optargs('--separate-layers', type='inkbool', dest='separate_layers', default=True, help='Seaparate gcode file per layer'), 
    svg.optargs('--create-log', type='inkbool', dest='log_create_log', default=True, help='Create log files'),
    svg.optargs('--log-level', type='string', dest='log_level', default='DEBUG', help='Log level'),
    svg.optargs('--log-filename', type='string', dest='log_filename', default='tcnc.log', help='Full pathname of log file'),
   
    svg.optargs('--preview-show', type='inkbool', dest='preview_show', default=True, help='Show generated cut paths on preview layer.'),
    svg.optargs('--debug-layer', type='inkbool', dest='debug_layer', default=True, help='Create debug layer.'),
    svg.optargs('--debug-biarcs', type='inkbool', dest='debug_biarcs', default=True),
   
    # Unused lagacy options - derived from gcodetools
    # TODO: use or delete
    svg.optargs('--z-depth', type='float', convert_to='world', default=-0.125, help='Z full depth of cut'),
    svg.optargs('--z-step', type='float', convert_to='world', default=-0.125, help='Z cutting step depth'),
    svg.optargs('--path-to-gcode-order',type='string', dest='path_to_gcode_order', default='path by path', help='Defines cutting order path by path or layer by layer.'), 
    svg.optargs('--path-to-gcode-depth-function', type='string', dest='path_to_gcode_depth_function', default='zd', help='Path to gcode depth function.'),
]

tcnc = TCnc(option_info)
tcnc.affect()
