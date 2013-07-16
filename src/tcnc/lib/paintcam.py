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
import math
import logging

import geom
import paths
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
    
    brush_depth = 0.0
    brush_overshoot = 0.0
    brush_liftoff_angle = 0.0
    brush_liftoff_height = 0.0
    brush_landing_angle = 0.0
    brush_landing_start_height = 0.0
    brush_landing_end_height = 0.0
    brush_flip_after_every_path = False
    brush_flip_before_reload = False
    brush_reload_before_path = False
    brush_reload_after_interval = False
    brush_reload_angle = 0.0
    brush_reload_dwell = 0.0
    # Brush trail is the distance the brush bristles trail the handle center
    # when in contact with the painting surface. Usually a function of tool
    # depth, bristle length, and bristle stiffness.
    brush_trail_down = 0.0
    fillet_radius = None
#    min_fillet_angle = math.pi * .5
#    use_fillets = True
    min_corner_angle = math.pi * .3
    close_polygons = True
    allow_corner_reversals = False
    
    def __init__(self, gcode, **kwargs):
        """
        :gcode: a GCode instance
        """
        super(PaintCAM, self).__init__(gcode, **kwargs)
#        self._first_liftoff = True
        # Keep track of brush flip state
#        self._brush_flip_toggle = -1
        # Uses spindle on/off command to lower/raise pneumatics.
        self.gc.tool_down_prefix = 'M3 S1'
        self.gc.tool_up_suffix = 'M5'
    
    def preprocess_cutpaths(self, cutpath_list):
        """Preprocess the cutpath list before it is converted to G code.
        This attempts to create a fillet where two lines meet.
        """
        cutpath_list = super(PaintCAM, self).preprocess_cutpaths(cutpath_list)
        if self.fillet_radius is None:
            self.fillet_radius = math.hypot(self.tool_width/2, self.brush_trail_down)
        # Do this in two passes in case the corner processor needs
        # to break up the path
        new_cutpath_list1 = []
        for cutpath in cutpath_list:
            new_cutpaths = self._cutpath_process_corners(cutpath)
            new_cutpath_list1.extend(new_cutpaths)
        new_cutpath_list2 = []
        for cutpath in new_cutpath_list1:
            new_cutpath = self._cutpath_process_pathends(cutpath)
            new_cutpath_list2.append(new_cutpath)
#        new_cutpath_list.sort(key=lambda cp: len(cp))
        return new_cutpath_list2
    
    def _cutpath_process_pathends(self, cutpath):
        # Brush landing
        segment = cutpath[0]
        if not isinstance(segment, geom.Line) and not isinstance(segment, geom.Arc):
            logging.debug('segment type=%s' % segment.__name__)
        assert(isinstance(segment, geom.Line) or isinstance(segment, geom.Arc))
        d = max(self.brush_trail_down, 0.01)
        if segment.length() > d:
            seg1, seg2 = segment.subdivide(d / segment.length())
            cutpath[0] = seg1
            cutpath.insert(1, seg2)
        cutpath[0].inline_z = self.brush_depth        
        # Brush liftoff
        segment = cutpath[-1]
        assert(isinstance(segment, geom.Line) or isinstance(segment, geom.Arc))
        brush_direction = getattr(segment, 'inline_end_angle', segment.end_tangent_angle())
        # Default overshoot distance works fine for 90d angles
        overshoot_dist = self.brush_trail_down#self.tool_width/2
        # Closed path?
        if geom.float_eq(segment.p2, cutpath[0].p1):
            # tangent to reverse brush direction
            t1 = geom.P.from_polar(1.0, brush_direction + math.pi)
            corner_angle = abs(cutpath[0].p1.angle2(cutpath[0].p1 + t1, cutpath[0].p2))
            if not geom.float_eq(corner_angle, math.pi):
#                overshoot_dist = self.brush_trail_down
#            else:
                # Calculate overshoot based on corner angle
                d = (math.pi - (corner_angle % math.pi)) / math.pi
                overshoot_dist *= d
                logger.debug('overshoot dist: %.4f, %.4f, %.4f' % (math.degrees(corner_angle), d, overshoot_dist))
#        logger.debug('brush_direction= %.4f [%.4f]' % (math.degrees(brush_direction), math.degrees(segment.end_tangent_angle())))
        delta = geom.P(math.cos(brush_direction), math.sin(brush_direction))
        overshoot_endp = segment.p2 + (delta * overshoot_dist)
        overshoot_line = geom.Line(segment.p2, overshoot_endp)
        liftoff_dist = self.brush_trail_down
        liftoff_endp = overshoot_endp + (delta * liftoff_dist)
        liftoff_line = geom.Line(overshoot_endp, liftoff_endp)
        liftoff_line.inline_z = 0.0
        cutpath.append(overshoot_line)
        cutpath.append(liftoff_line)
        return cutpath

    def _corner_angle(self, seg1, seg2):
        t1 = geom.P.from_polar(1.0, seg1.end_tangent_angle() + math.pi)
        return seg1.p1.angle2(seg1.p1 + t1, seg1.p2)
    
    def _cutpath_process_corners(self, cutpath):
        if len(cutpath) <= 1:
            return [cutpath,]
        cutpath_list = []
        new_cutpath = paths.Path(name=getattr(cutpath, 'name'))
        firstline = cutpath[0]
        line1 = firstline
        for line2 in cutpath[1:]:
            # TODO: also process curves using tangent lines
            if isinstance(line1, geom.Line) and isinstance(line2, geom.Line):
                corner_angle = line1.p2.angle2(line1.p1, line2.p2)
                new_segments = ()
                if self.allow_corner_reversals and abs(corner_angle) <= self.min_corner_angle:
                    new_segments = self._corner_reversal(line1, line2, corner_angle)
                elif abs(corner_angle) > self.min_corner_angle and not geom.float_eq(abs(corner_angle), math.pi):
                    new_segments = self._corner_turn(line1, line2, corner_angle)
                if new_segments:
                    new_cutpath.extend(new_segments[:-1])
                    line2 = new_segments[-1]
                else:
                    new_cutpath.append(line1)
                    if not geom.float_eq(abs(corner_angle), math.pi):
                        # Split this path if the corner can't be processed
                        cutpath_list.append(new_cutpath)
                        new_cutpath = paths.Path()
            else:
                new_cutpath.append(line1)
            line1 = line2
        new_cutpath.append(line2)
        # process last segment of polygon
        if self.close_polygons and new_cutpath.is_closed() \
        and isinstance(new_cutpath[-1], geom.Line) \
        and isinstance(firstline, geom.Line):
#        and geom.float_eq(new_cutpath[0].p1, new_cutpath[-1].p2):
            line1 = new_cutpath[-1]
            line2 = firstline
            corner_angle = line1.p2.angle2(line1.p1, line2.p2)
            new_segments = ()
            if self.allow_corner_reversals and abs(corner_angle) <= self.min_corner_angle:
                new_segments = self._corner_reversal(line1, line2, corner_angle)
            elif abs(corner_angle) > self.min_corner_angle and not geom.float_eq(abs(corner_angle), math.pi):
                new_segments = self._corner_turn(line1, line2, corner_angle)
            if new_segments:
                new_cutpath[0] = new_segments[-1]
                del new_cutpath[-1]
                new_cutpath.extend(new_segments[0:-1])
        if len(new_cutpath) > 0:
            cutpath_list.append(new_cutpath)
        return cutpath_list

    def _corner_reversal(self, line1, line2, corner_angle):
        # For highly acute angles add a segment the size of the brush
        # trail and reverse brush direction.
        xline1_angle = line1.angle() + line1.p2.angle2(line1.p1, line2.p2)/2
        xline1 = geom.Line.from_polar(line1.p2, self.brush_trail_down, xline1_angle)
        xline1.inline_z = self.brush_depth / 8
        xline1.inline_delta_a = corner_angle / 2
        xline1.inline_ignore_angle = True
        xline2 = xline1.reversed()
        xline2.inline_z = self.brush_depth
        xline2.inline_delta_a = corner_angle / 2
        xline2.inline_flip_tool = True
#            line2.inline_ignore_angle = True
        return (line1, xline1, xline2, line2)

    def _corner_reversal2(self, line1, line2, corner_angle):
        # For highly acute angles add a segment the size of the brush
        # trail and reverse brush direction.
        xline1_angle = line1.angle() + line1.p2.angle2(line1.p1, line2.p2)/2
        xline1 = geom.Line.from_polar(line1.p2, self.brush_trail_down, xline1_angle)
        xline1.inline_z = self.brush_depth / 8
        xline1.inline_delta_a = corner_angle / 2
        xline2 = xline1.reversed()
        xline2.inline_z = self.brush_depth
#            line2.inline_ignore_angle = True
        return (line1, xline1, xline2, line2)
            
    def _corner_turn(self, line1, line2, corner_angle):
        """Try to replace the two segments with shorter segments
        connected with a fillet arc that compensates for brush trail.
        Returns a list (segment1, arc, segment2).
        Returns an empty list if the segments could not be connected with
        a fillet."""
        # Determine the fillet arc center which is the
        # intersection of two interior offset lines
#        fillet_overlap = self.tool_width/3
        lineside = line1.which_side(line2.p2)
        offset = (self.tool_width/2) * lineside
        offset_line1 = line1.offset(offset)
        offset_line2 = line2.offset(offset)
        fillet_center = offset_line1.intersection(offset_line2)
#        logger.debug('fillet_center: x=%.4f, y=%.4f' % (fillet_center.x, fillet_center.y))
#        self.preview_plotter.draw_interval_marker(fillet_center)
        u1 = line1.projection(fillet_center)
        u2 = line2.projection(fillet_center)
        ux1 = self.brush_trail_down / line1.length()
        ux2 = self.brush_trail_down / line2.length()
        u1 += ux1
        u2 += ux2
        if (0.0 < u1 < (1.0+ux1)) and (0.0 < u2 < 1.0):#(1.0+ux2)):
            p1 = line1.p1 + ((line1.p2 - line1.p1) * u1)
            p2 = line2.p1 + ((line2.p2 - line2.p1) * u2)
            fillet_angle = fillet_center.angle2(p1, p2)
            fillet_arc = geom.Arc(p1, p2, self.fillet_radius, fillet_angle, fillet_center)
            xline1 = geom.Line(line1.p1, fillet_arc.p1)
            xline2 = geom.Line(fillet_arc.p2, line2.p2)
            fillet_arc.inline_ignore_angle = True
            fillet_arc.inline_end_angle = xline2.start_tangent_angle()
            return (xline1, fillet_arc, xline2)
        else:
            logger.debug('corner turn too tight')
            logger.debug('corner_angle= %.4f' % (math.degrees(corner_angle),))
            logger.debug('epsilon= %f, angle=%f' % (geom.EPSILON, corner_angle))
            logger.debug('u1=%.4f, ux1=%.4f' % (u1, ux1))
            logger.debug('u2=%.4f, ux2=%.4f' % (u2, ux2))
            return ()

    def generate_tool_landing_gcode(self, tool_direction, depth, start_path=False):
        """Generate the G code for the brush landing."""
        if start_path and self.brush_reload_before_path:
            self.generate_brush_reload_gcode(tool_direction)
#        # Calculate the point to backup the brush to where the descent starts
#        delta = geom.P(math.cos(tool_direction), math.sin(tool_direction))
#        xydist = abs(self.brush_depth) / math.tan(self.brush_landing_angle)
#        startp = self.last_point - (delta * xydist)
#        self.gc.comment('Brush landing')
#        # Move to beginning of descent
#        self.gc.rapid_move(startp.x, startp.y)
        # Bring the brush straight down to the start of the angled descent
        self.gc.tool_down(self.brush_landing_start_height, comment=None)
#        # Move the brush at the desired downward angle towards the paper
#        self.gc.feed(self.last_point.x, self.last_point.y, self.brush_depth)
    
    def generate_tool_liftoff_gcode(self, tool_direction, end_path=False):
        """Generate the G code and preview for the brush liftoff.
        Also generates the overshoot vector.
        The liftoff vector is determined by the current brush direction,
        the lift-off angle, and the height.
        Once the brush is lifted off then it is brought up to the safe height.
        """
        #logger.debug('tool direction= %f' % tool_direction)
#        delta = geom.P(math.cos(tool_direction), math.sin(tool_direction))
#        overshoot_endp = self.last_point + (delta * self.brush_overshoot)
#        self.gc.feed(overshoot_endp.x, overshoot_endp.y, comment='brush overshoot')
#        self.preview_plotter.draw_feed_line(self.last_point, overshoot_endp,
#                                            self.tool_width)
#        self.last_point = overshoot_endp
#        self.feed_distance += self.brush_overshoot
#        liftoff_dist = abs(self.brush_depth) / math.tan(self.brush_liftoff_angle)
#        liftoff_endp = self.last_point + (delta * liftoff_dist)
#        self.gc.comment('Brush lift-off')
#        self.gc.feed(liftoff_endp.x, liftoff_endp.y, self.brush_liftoff_height)
        self.gc.tool_up(comment=None)
        if end_path and self.brush_flip_after_every_path:
            self.flip_tool()
            
    def generate_feed_interval_gcode(self, tool_direction, depth):
        """Generate G code to reload brush after the feed travel has passed
        a threshold SimpleCAM.:feed_interval:
        :tool_direction: current tool angle.
        :depth: current cut depth
        """
        self.generate_tool_liftoff_gcode(tool_direction, start_path=False)
        if self.brush_flip_before_reload:
            self.flip_tool()
        if self.brush_reload_after_interval:
            self.generate_brush_reload_gcode(tool_direction)
        self.generate_tool_landing_gcode(tool_direction, depth, end_path=False)
            
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

#    def flip_brush(self):
#        """Offset tangential tool rotation by 180deg.
#        This useful for brush-type or double sided tools to even out wear.
#        """
#        self._tool_flip_toggle *= -1 # Toggle rotation direction
#        self.gc.offset['A'] += self._tool_flip_toggle * 180
        
