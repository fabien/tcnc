'''An Inkscape extension that will output G-code from
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
'''
import os
import fnmatch
import logging
import math
import gettext
_ = gettext.gettext

import inkex
import simpletransform

from lib import svg
from lib import geom
from lib import gcode

from lib.geom import P, Line

VERSION = "0.1"

DEFAULT_LOG_LEVEL = logging.DEBUG
logger = logging.getLogger(__name__)

PREVIEW_LAYER_NAME = 'wash-preview'

class Wash(svg.SuperEffect):
    '''Inkscape plugin that converts selected SVG elements into gcode suitable for a
    four axis (XYZA) CNC machine with a tangential tool (ie a knife or a brush) as the A axis.
    '''
    styles = {
              'simple': 'fill:none;stroke:#cccc99;stroke-width:0.25pt;stroke-opacity:1',
              'feedline': 'fill:none;stroke:#c000c0;stroke-width:0.75pt;stroke-opacity:1;marker-end:url(#PreviewLineEnd0)',
              'feedline1': 'fill:none;stroke:#ff00ff;stroke-width:1pt;stroke-opacity:1;',
              'cutpath_end_marker': 'fill-rule:evenodd;fill:#FFFFFF;stroke:#ff0000;stroke-width:1.0pt;marker-start:none',
              'movepath_end_marker': 'fill-rule:evenodd;fill:#00ff00;stroke:#00ff00;stroke-width:1.0pt;marker-start:none',
              'tangent_tool': 'fill:none;stroke:#00cc00;stroke-width:2.0pt;stroke-opacity:.5',
              'moveline': 'fill:none;stroke:#00ff00;stroke-width:1pt;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;stroke-miterlimit:4;stroke-dasharray:1.2, 5;stroke-dashoffset:0;marker-end:url(#PreviewLineEnd2)',
              }
        
    last_point = (0.0, 0.0)
            
    def effect(self):
        '''Main entry point for Inkscape plugins.
        '''
        # Create a log file for debugging
        if self.options.log_create_log and self.options.log_filename:
            log_path = os.path.abspath(os.path.expanduser(self.options.log_filename))
            log_level = getattr(logging, self.options.log_level, DEFAULT_LOG_LEVEL)
            logging.basicConfig(filename=log_path, filemode='w', level=log_level)
        
        if self.options.units == 'doc':
            self.units = self.get_document_units()
        else:
            self.units = self.options.units
            
        # Convert units
        unit_scale = inkex.uuconv[self.units]
        self.options.line_spacing *= unit_scale
        self.options.brush_size *= unit_scale
        self.options.brush_overlap *= unit_scale
        self.options.brush_overshoot *= unit_scale
        self.options.margin_right *= unit_scale
        self.options.margin_left *= unit_scale
        self.options.margin_top *= unit_scale
        self.options.margin_bottom *= unit_scale
        self.options.brush_reload_angle = math.pi
        self.options.brush_depth = -1.0 * unit_scale
               
        self.styles['brushline'] = \
            'fill:none;stroke:#c0c0ff;stroke-width:%dpx;stroke-opacity:.5' \
            % self.options.brush_size

        # Create the SVG for the line end markers used in preview layer
        self.create_inkscape_markers()
            
        # Clear out the old preview layer if any
        self.clear_layer(PREVIEW_LAYER_NAME)

        # Create a new layer that will contain the G code preview
        self.preview_layer = self.create_layer(PREVIEW_LAYER_NAME)       
        self.current_layer = self.preview_layer        

        # Add a transform to the preview layer to flip the coordinates
        # from cartesian to SVG (flip Y axis from lower left to upper left).
        page_height = self.get_document_size()[1]
        flip_transform_attr = 'translate(0, %f) scale(1, -1)' % page_height
        self.preview_layer.set('transform', flip_transform_attr)
#        flip_transform = simpletransform.parseTransform(flip_transform_attr)
        
        if self.options.gen_gcode:
            self.gc = gcode.GCode(xyfeed=self.options.brush_speed,
                                  zwait=self.options.brush_start_delay,
                                  zsafe=1.0 * unit_scale,
                                  afeed=self.options.brush_rotate_speed * 60.0,)
            self.gc.units = self.units
            self.gc.set_unit_scale(1.0 / inkex.uuconv[self.units])
            self.gc.default_header(units=self.units,
                                   description=('Generated by TCNC/Wash Inkscape extension version %s' % VERSION,))
            # Move off the paper and reset the Z axis
            self.gc.home(axes='X')
            self.gc.tool_down()
            self.gc.tool_up()
            
        last_point = P(0,0)
        if self.options.line_orient == 'h':
            self.generate_horizontal_lines(last_point)
        else:
            self.generate_vertical_lines(last_point)
        
        if self.options.gen_gcode:
            # Add a default G code footer
            self.gc.home(axes='XYA')
            self.gc.default_footer()
            # Export the G code to a file
            self.gc.export(self.options.directory, self.options.filename,
                           self.options.append_suffix)
        
    def generate_horizontal_lines(self, last_point):
        if self.options.brush_size > 0.0:
            line_spacing = self.options.brush_size - self.options.brush_overlap
        else:
            line_spacing = self.options.line_spacing
        min_x = self.options.margin_left
        max_x = self.get_document_width() - self.options.margin_right
        min_y = self.options.margin_bottom
        max_y = self.get_document_height() - self.options.margin_top
        if self.options.line_top2bottom:
            line_offset = max_y - (self.options.brush_size / 2)
            if self.options.brush_overlap > 0:
                line_offset += self.options.brush_overlap
            line_spacing = -line_spacing
        else:
            line_offset = min_y + (self.options.brush_size / 2)
            if self.options.brush_overlap > 0:
                line_offset -= self.options.brush_overlap
        self.render_horizontal_lines(last_point, min_x, max_x, min_y, max_y,
                                    line_offset, line_spacing,
                                    self.options.line_left2right,
                                    self.options.line_alt)
#        if self.options.line_alt:
#            self.render_horizontal_lines(last_point, min_x, max_x, min_y, max_y,
#                                         line_offset, line_spacing * 2,
#                                         self.options.line_left2right)
#            line_offset += line_spacing
#            self.render_horizontal_lines(last_point, min_x, max_x, min_y, max_y,
#                                         line_offset, line_spacing * 2,
#                                         not self.options.line_left2right)
#        else:
#            self.render_horizontal_lines(last_point, min_x, max_x, min_y, max_y,
#                                        line_offset, line_spacing,
#                                        self.options.line_left2right)

    def render_horizontal_lines(self, last_point, min_x, max_x, min_y, max_y,
                                line_offset, line_spacing, left2right, alternate=False):
        line_count = 0
        line_start = self.options.line_start - 1
        line_skip = self.options.line_skip + 1
        while line_offset >= min_y and line_offset <= max_y:
            if left2right:
                p1 = P(min_x, line_offset)
                p2 = P(max_x + self.options.brush_overshoot, line_offset)
                p2b = P(max_x, line_offset)
            else:
                p1 = P(max_x, line_offset)
                p2 = P(min_x - self.options.brush_overshoot, line_offset)
                p2b = P(min_x, line_offset)
            if alternate:
                left2right = not left2right
            if line_count >= (line_start - 1) \
            and (line_skip == 1 or (line_count - line_start) % line_skip == 0):
                feedline = Line(p1, p2)
                brushline = Line(p1, p2b)
                last_point = self.render_line(last_point, feedline, brushline)
            line_offset += line_spacing
            line_count += 1
        return last_point
            
    def generate_vertical_lines(self, last_point):
        if self.options.brush_size > 0.0:
            line_spacing = self.options.brush_size - self.options.brush_overlap
        else:
            line_spacing = self.options.line_spacing
        min_x = self.options.margin_left
        max_x = self.get_document_width() - self.options.margin_right
        min_y = self.options.margin_bottom
        max_y = self.get_document_height() - self.options.margin_top
        if self.options.line_left2right:
            line_offset = min_x + (self.options.brush_size / 2)
            if self.options.brush_overlap > 0:
                line_offset -= self.options.brush_overlap
        else:
            line_offset = max_x - (self.options.brush_size / 2)
            line_spacing = -line_spacing
            if self.options.brush_overlap > 0:
                line_offset += self.options.brush_overlap
        self.render_vertical_lines(last_point, min_x, max_x, min_y, max_y,
                                    line_offset, line_spacing,
                                    self.options.line_top2bottom,
                                    self.options.line_alt)
            
    def render_vertical_lines(self, last_point, min_x, max_x, min_y, max_y,
                                line_offset, line_spacing, top2bottom, alternate=False):
        line_count = 0
        line_start = self.options.line_start - 1
        line_skip = self.options.line_skip + 1
        while line_offset >= min_x and line_offset <= max_x:
            if top2bottom:
                p1 = P(line_offset, max_y)
                p2 = P(line_offset, min_y - self.options.brush_overshoot)
                p2b = P(line_offset, min_y)
            else:
                p1 = P(line_offset, min_y)
                p2 = P(line_offset, max_y + self.options.brush_overshoot)
                p2b = P(line_offset, max_y)
            if alternate:
                top2bottom = not top2bottom
            if line_count >= (line_start - 1) \
            and (line_skip == 1 or (line_count - line_start) % line_skip == 0):
                feedline = Line(p1, p2)
                brushline = Line(p1, p2b)
                last_point = self.render_line(last_point, feedline, brushline)
            line_offset += line_spacing
            line_count += 1
        return last_point
            
    def render_line(self, last_point, feedline, brushline):
        ''''''
        if self.options.gen_gcode:
            self.gc.tool_up()
            # Rotate the brush to the reload position
            self.gc.rapid_rotate(self.options.brush_reload_angle,
                                 comment='Rotate brush to reload')
            # Move to beginning of line
            self.gc.rapid_move(feedline.p1.x, feedline.p1.y,
                               comment='move to beginning of line')
            # Pause to add more paint to the brush
            if self.options.brush_pause:
                self.gc.pause(comment='Pause for brush reload')
            elif self.options.brush_dwell > 0.0:
                self.gc.dwell(self.options.brush_dwell * 1000)
            # Rotate brush to direction of line travel
            self.gc.rapid_rotate(feedline.angle() + (math.pi/2),
                                 comment='Rotate brush to line direction')
            # Lower the brush and paint to end of line
            self.gc.tool_down(self.options.brush_depth, comment='lower brush to paper')
            self.gc.feed(feedline.p2.x, feedline.p2.y, comment='draw line')
            # Flip the brush 180deg if specified
            if self.options.brush_flip_reload:
                self.flip_brush(self.gc)

        # Draw a dot on the preview layer to show where a brush reload will occur
        self.draw_preview_dot(feedline.p1, color='#0000ff')
        # Render the preview lines
        moveline = Line(last_point, feedline.p1)
        self.draw_preview_line(moveline, self.styles['moveline'])
        self.draw_preview_line(feedline, self.styles['feedline'])
        self.draw_preview_line(brushline, self.styles['brushline'])
        return feedline.p2 # Last feed point
            
    def create_inkscape_markers(self):
        '''Create Inkscape line end marker glyphs and insert them into the document.
        '''
        self.create_simple_marker('PreviewLineEnd0', 'M 5.77,0.0 L -2.88,5.0 L -2.88,-5.0 L 5.77,0.0 z',
                          self.styles['cutpath_end_marker'], 'scale(0.5) translate(-4.5,0)')
        self.create_simple_marker('PreviewLineEnd1', 'M 5.77,0.0 L -2.88,5.0 L -2.88,-5.0 L 5.77,0.0 z',
                          self.styles['cutpath_end_marker'], 'scale(-0.5) translate(-4.5,0)')
        self.create_simple_marker('PreviewLineEnd2', 'M 5.77,0.0 L -2.88,5.0 L -2.88,-5.0 L 5.77,0.0 z',
                          self.styles['movepath_end_marker'], 'scale(0.5) translate(-4.5,0)')


    def draw_preview_line(self, line, style):
        '''Draw an SVG line path on to the preview layer'''
        self.create_line(line.p1.x, line.p1.y, line.p2.x, line.p2.y, style)

    def draw_preview_dot(self, point, size='medium', color='#000000'):
        '''Draw an SVG dot (a small filled circle) at the specified location'''
        radius = {'small': '2pt', 'medium': '5pt', 'large': '10pt',}
        style = 'fill:%s;stroke:%s' % (color, color)
        self.create_circle(point.x, point.y, radius[size], style)
        
    brush_flipped = -1
    def flip_brush(self, gc):
        # Offset brush rotation by 180deg.
        self.brush_flipped *= -1 # Toggle rotation direction
        gc.offset['A'] += self.brush_flipped * 180

    def generate_brush_reload_gcode(self, stop_point, current_angle=None):
        '''Generate the G code for reloading the brush tool.'''
        # Rotate the brush to a position that makes it easy to add paint
        if current_angle is None:
            self.gc.rapid_rotate(self.options.brush_reload_angle, comment='rotate brush to reload')
        else:
            # Calculate smallest angle offset to get to reload angle
            # from current angle
            da = math.fmod(current_angle, 2*math.pi)
            if math.fabs(da) < math.pi:
                a = current_angle - da
            elif da < 0:
                a = current_angle - (2*math.pi + da)
            else:
                a = current_angle + (2*math.pi - da)
            reload_angle = self.options.brush_reload_angle + a
            self.gc.rapid_rotate(reload_angle)
    
option_info = (
    ('--active-tab', '', 'store', 'string', 'active_tab', '', ''),
    
    ('--units', '', 'store', 'string', 'units', 'in', 'Document units.'),
   
    ('--line-orient', '', 'store', 'string', 'line_orient', 'h', 'Line orientation'),
    ('--line-spacing', '', 'store', 'float', 'line_spacing', '1', 'Line spacing'),
    ('--line-angle', '', 'store', 'float', 'line_angle', '1', 'Line angle'),
    ('--line-left2right', '', 'store', 'inkbool', 'line_left2right', True, 'Draw left to right'),
    ('--line-top2bottom', '', 'store', 'inkbool', 'line_top2bottom', True, 'Draw top to bottom'),
    ('--line-alt', '', 'store', 'inkbool', 'line_alt', False, 'Draw alternate lines'),
    ('--line-alt-rev', '', 'store', 'inkbool', 'line_alt_rev', False, 'Draw alternate lines in reverse'),
    ('--line-skip', '', 'store', 'int', 'line_skip', '0', 'Skip lines'),
    ('--line-start', '', 'store', 'int', 'line_start', '0', 'Start lines at'),

    ('--brush-size', '', 'store', 'float', 'brush_size', '1', 'Brush size'),
    ('--brush-overlap', '', 'store', 'float', 'brush_overlap', '0', 'Brush overlap'),
    ('--brush-overshoot', '', 'store', 'float', 'brush_overshoot', '0', 'Brushstroke overshoot distance.'),
    ('--brush-speed', '', 'store', 'float', 'brush_speed', '100', 'Brush speed.'),
    ('--brush-rotate-speed', '', 'store', 'float', 'brush_rotate_speed', '100', 'Brush rotation speed.'),
    ('--brush-start-delay', '', 'store', 'float', 'brush_start_delay', '500', 'Brush start delay.'),
    ('--brush-flip-reload', '', 'store', 'inkbool', 'brush_flip_reload', True, 'Flip after reload.'),
    ('--brush-dwell', '', 'store', 'float', 'brush_dwell', '0', 'Brush reload time (seconds).'),
    ('--brush-pause', '', 'store', 'inkbool', 'brush_pause', False, 'Pause to reload (overrides reload time).'),
    
    ('--brush-flip-path', '', 'store', 'inkbool', 'brush_flip_path', True, 'Flip after path.'),
    ('--brush-reload', '', 'store', 'inkbool', 'brush_reload', True, 'Enable brush reload.'),
    ('--brush-reload-path', '', 'store', 'inkbool', 'brush_reload_path', True, 'Force brush reload every path.'),
#    ('--brushstroke-max', '', 'store', 'float', 'brushstroke_max', '10', 'Maximum brushstroke distance.'),
#    ('--brushstroke-overlap', '', 'store', 'float', 'brushstroke_overlap', '0', 'Brushstroke overlap.'),
    ('--brush-reload-angle', '', 'store', 'float', 'brush_reload_angle', 90.0, 'Brush reload angle (degrees).'),

    ('--margin-left', '', 'store', 'float', 'margin_left', '1', 'Left margin'),
    ('--margin-right', '', 'store', 'float', 'margin_right', '1', 'Right margin'),
    ('--margin-top', '', 'store', 'float', 'margin_top', '1', 'Top margin'),
    ('--margin-bottom', '', 'store', 'float', 'margin_bottom', '1', 'Bottom margin'),

    ('--xy-feed', '', 'store', 'float', 'xy_feed', '10.0', 'XY axis feed rate in unit/s'),
    ('--z-feed', '', 'store', 'float', 'z_feed', '10.0', 'Z axis feed rate in unit/s'),
    ('--a-feed', '', 'store', 'float', 'a_feed', '60.0', 'A axis feed rate in deg/s'),
    ('--z-safe', '', 'store', 'float', 'z_safe', '5.0', 'Z axis safe height for rapid moves'),
    ('--z-wait', '', 'store', 'float', 'z_wait', '500', 'Z axis wait (milliseconds)'),
   
    ('--gen-gcode', '', 'store', 'inkbool', 'gen_gcode', True, 'Generate G code'), 
    ('--directory', '', 'store', 'string', 'directory', '~', 'Directory for gcode file'),
    ('--filename', '', 'store', 'string', 'filename', 'wash', 'File name'), 
    ('--append-suffix', '', 'store', 'inkbool', 'append_suffix', True, 'Append auto-incremented numeric suffix to filename'), 
    ('--create-log', '', 'store', 'inkbool', 'log_create_log', True, 'Create log files'),
    ('--log-level', '', 'store', 'string', 'log_level', 'DEBUG', 'Log level'),
    ('--log-filename', '', 'store', 'string', 'log_filename', 'wash.log', 'Full pathname of log file'),
   
    ('--preview-show', '', 'store', 'inkbool', 'preview_show', True, 'Show generated cut paths on preview layer.'),
)

wash = Wash(option_info)
wash.affect()
