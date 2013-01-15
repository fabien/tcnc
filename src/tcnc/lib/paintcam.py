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
import logging
import math

import geom
import simplecam

VERSION = "0.1"

DEFAULT_LOG_LEVEL = logging.DEBUG
logger = logging.getLogger(__name__)

class PaintCAM(simplecam.SimpleCAM):
    """2D CAM library that converts selected SVG elements into gcode
    suitable for a straightforward three axis machine with an optional
    fourth angular axis (A) that rotates about the Z axis. The fourth axis
    position is always tangential to the movement along the X and Y axes.
    This is usually called a tangential tool (ie a knife or a brush).
    """
    
    brush_overshoot = 0.0
    brush_liftoff_angle = 0.0
    brush_liftoff_height = 1.0
    brush_landing_angle = 0.0
    brush_landing_start_height = 1.0
    brush_landing_end_height = 0.0
    brush_flip_after_every_path = False
    brush_flip_before_reload = False
    brush_reload_angle = 0.0
    brush_reload_dwell = 0.0
    
    def __init__(self, gcode, **kwargs):
        """
        :gcode: a GCode instance
        """
        super(PaintCAM, self).__init__(gcode, **kwargs)
        self._first_liftoff = True
        # Keep track of brush flip state
        self._brush_flip_toggle = -1
    
    def generate_tool_landing_gcode(self, tool_direction, depth, start_path=False):
        """Generate the G code for the brush landing."""
        if start_path:
            self.generate_brush_reload_gcode(tool_direction)
        delta = geom.P(math.cos(tool_direction), math.sin(tool_direction))
        landing_distance = (self.brush_landing_start_height -
                            self.brush_landing_end_height)
        xydist = landing_distance / math.tan(self.brush_landing_angle)
        startp = self.last_point - (delta * xydist)
        self.gc.comment('Brush landing')
        self.gc.rapid_move(startp.x, startp.y)
        self.gc.tool_down(self.brush_landing_start_height, comment=None)
        self.gc.feed(self.last_point.x, self.last_point.y,
                     self.brush_landing_end_height)
    
    def generate_tool_liftoff_gcode(self, tool_direction, end_path=False):
        """Generate the G code and preview for the brush liftoff.
        Also generates the overshoot vector.
        The liftoff vector is determined by the current brush direction,
        the lift-off angle, and the height.
        Once the brush is lifted off then it is brought up to the safe height.
        """
        delta = geom.P(math.cos(tool_direction), math.sin(tool_direction))
        overshoot_endp = self.last_point + (delta * self.brush_overshoot)
        self.gc.feed(overshoot_endp.x, overshoot_endp.y)
        self.feed_distance += self.brush_overshoot
        liftoff_dist = self.brush_liftoff_height / math.tan(self.brush_liftoff_angle)
        liftoff_endp = self.last_point + (delta * liftoff_dist)
        self.gc.comment('Brush lift-off')
        self.gc.feed(liftoff_endp.x, liftoff_endp.y, self.brush_liftoff_height)
        self.gc.tool_up(comment=None)
        self.preview_plotter.draw_feed_line(self.last_point, overshoot_endp,
                                            self.tool_width)
        if end_path and self.brush_flip_after_every_path:
            self.flip_brush()
            
    def generate_feed_interval_gcode(self, tool_direction, depth):
        """Generate G code to reload brush after the feed travel has passed
        a threshold SimpleCAM.:feed_interval:
        :tool_direction: current tool angle.
        :depth: current cut depth
        """
        self.generate_tool_liftoff_gcode(tool_direction, start_path=False)
        if self.brush_flip_before_reload:
            self.flip_brush()
        self.generate_brush_reload_gcode(tool_direction)
        self.generate_tool_landing_gcode(tool_direction, depth, start_path=False)
            
    def generate_brush_reload_gcode(self, current_angle):
        """Generate the G code for reloading the brush tool."""
        self.gc.comment('Rotate brush and pause for reload')
        # Rotate the brush to a position that makes it easy to add paint
        da = math.fmod(current_angle, 2*math.pi)
        if math.fabs(da) < math.pi:
            a = current_angle - da
        elif da < 0:
            a = current_angle - (2*math.pi + da)
        else:
            a = current_angle + (2*math.pi - da)
        self.gc.rapid_rotate(self.brush_reload_angle + a)
        # Pause to add more paint to the brush
        self.gc.dwell(self.brush_reload_dwell * 1000)
        # TODO: backup brush to overlap previous stroke
        # Rotate brush back to drawing angle
        self.gc.rapid_rotate(current_angle)

    def flip_brush(self):
        """Offset tangential tool rotation by 180deg.
        This useful for brush-type or double sided tools to even out wear.
        """
        self._tool_flip_toggle *= -1 # Toggle rotation direction
        self.gc.offset['A'] += self._tool_flip_toggle * 180
        
