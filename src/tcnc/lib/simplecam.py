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
import fnmatch
import logging
import math

import geom
import gcode

VERSION = "0.1"

DEFAULT_LOG_LEVEL = logging.DEBUG
logger = logging.getLogger(__name__)

class SimpleCAM(object):
    """Simple 2D CAM library that converts simple line/arc geometry into gcode
    suitable for a straightforward three axis machine with an optional
    fourth angular axis (A) that rotates about the Z axis. The fourth axis
    position is always tangential to the movement along the X and Y axes.
    This is usually called a tangential tool (ie a knife or a brush).
    """
    
    class PreviewPlotter(object):
        """Base class for plotting CAM preview geometry.
        This should be subclassed by users of SimpleCAM.
        This default class does not draw anything.
        """
        def draw_move_line(self, p1, p2):
            """Draw a line from :p1: to :p2: on the preview layer
            that represents a rapid move.
            """
            pass
        
        def draw_feed_line(self, p1, p2, width=1.0, depth=0.0):
            """Draw a line from :p1: to :p2: on the preview layer
            that represents a linear feed.
            :width: Line (tool) width in machine units
            :depth: Current tool depth in machine units
            """
            pass
        
        def draw_feed_rotate(self, angle1, angle2, width=1.0, depth=0.0):
            """Draw a tool rotation from :angle1: to :angle2: on the
            preview layer.
            :width: Line (tool) width in machine units
            :depth: Current tool depth in machine units
            """
            pass
        
        def draw_feed_arc(self, arc, width=1.0, depth=0.0):
            """Draw an arc on the preview layer that represents a circular feed.
            :arc: A geom.Arc object
            :width: Line (tool) width in machine units
            :depth: Current tool depth in machine units
            """
            pass
        
        def draw_interval_marker(self, p, depth=0.0):
            """Draw a marker glyph (a small filled circle for example)
            at the specified location :p: on the preview layer.
            :depth: Current tool depth in machine units
            """
            pass
            
    # Sort paths to optimize/minimize rapid moves
    optimize_rapid_moves = True
    # Maximum cutting depth per pass
    cut_depth = 0.0
    # Final cutting depth
    final_cut_depth = 0.0
    # Tool width in machine units
    tool_width = 1.0
    # Feed distance to travel (in machine units) before outputting special
    # interval G code. 0 means no special output.
    feed_interval = 0.0
    
    # Biarc approximation tolerance
    biarc_tolerance = 0.01
    # Maximum bezier curve subdivision recursion
    biarc_max_depth = 4
    # Flatness of curve to convert to line
    biarc_line_flatness = 0.001
    
    def __init__(self, gcode, preview_plotter=None):
        """
        :gcode: a GCode instance
        :preview_plotter: an instance of PreviewPlotter.
        Default is a null plotter (no preview output)
        """
        self.gc = gcode
        if preview_plotter is None:
            preview_plotter = SimpleCAM.PreviewPlotter()
        self.preview_plotter = preview_plotter
        # Cumulative tool cutting distance
        self.feed_distance = 0.0
        # Last tool location on XY plane
        self.last_point = geom.P(0.0, 0.0)
        # Current tool angle
        self.current_angle = 0.0
        # Feed interval travel
        self.feed_interval_travel = 0.0
    
    def generate_gcode(self, cutpath_list):
        """Generate G code from cutpaths.
        :cutpath_list: A list of tuples of the form
        (<path_id>, <cutpath>) where :path_id: is a string
        and :cutpath: is a list of geom.Line or geom.Arc objects
        """
        if self.optimize_rapid_moves:
            cutpath_list = self.sort_cutpaths(cutpath_list)
        self.gc.default_header()
        self.gc.tool_up()
        self.gc.spindle_on()
            
        for path_id, cutpath in cutpath_list:
            self.gc.addline()
            if path_id:
                self.gc.comment('Path "%s"' % path_id)
            self.generate_cutpath_gcode(cutpath)
            
        self.gc.spindle_off()
        # Home the X, Y, and A axes.
        self.gc.home(axes='XYA')
        self.gc.default_footer()
        
    def generate_cutpath_gcode(self, cutpath, depth=0):
        """Generate G code for the given cutpath.
        This method also creates SVG output that serves as a preview to the
        G code output. The SVG output is added to the specified group or layer.
        The cutpath is specified by a list of curve segments which are either
        line segments or circular arc segments.
        :depth: The depth of cut for this pass
        """
        if not cutpath:
            return
        
        self.current_angle = cutpath[0].start_tangent()
        # Rapidly move to the beginning of the cut path
        self.gc.rapid_move(cutpath[0].p1.x, cutpath[0].p1.y, a=self.current_angle)
        # Draw the initial rapid move line on the SVG preview layer
        self.preview_plotter.draw_move_line(self.last_point, cutpath[0].p1)
        # Generate the G code to lower the tool to cut depth
        self.generate_tool_landing_gcode(self.current_angle, depth, start_path=True)
        self.last_point = cutpath[0].p1
        
        self.gc.comment('Start cutting path')
        # Create G-code for each segment of the cutpath
        for segment in cutpath:
            # Convert any CubicBeziers into Arcs and Lines
            # using biarc approximation.
            if isinstance(segment, geom.CubicBezier):
                biarcs = segment.biarc_approximation(tolerance=self.biarc_tolerance,
                                                     max_depth=self.biarc_max_depth,
                                                     line_flatness=self.biarc_line_flatness)
                for biarc_segment in biarcs:
                    self._generate_segment_gcode(biarc_segment, depth)
            else:
                self._generate_segment_gcode(segment, depth)

        # Generate G code to lift the tool to safe height
        self.generate_tool_liftoff_gcode(self.current_angle, end_path=True)
        # Do a fast unwind if angle is > 360deg
        if self.current_angle > (math.pi * 2):
            self.gc.rehome_rotational_axis()
        return

    def generate_feed_interval_gcode(self, tool_direction, depth):
        """Generate G code for some purpose after the feed travel has passed
        a threshold SimpleCAM.:feed_interval:
        This does nothing by default since it is a hook for subclasses.
        :tool_direction: current tool angle.
        :depth: current cut depth
        """
        pass
    
    def generate_tool_landing_gcode(self, tool_direction, depth, start_path=False):
        """Generate the G code and preview for the tool landing sequence.
        By default this just lowers the tool along the Z axis to the cutting depth
        :tool_direction: current tool angle.
        :depth: current cut depth
        :start_path: True if the current position is at the start of a path
        """
        self.gc.tool_down(depth)

    def generate_tool_liftoff_gcode(self, tool_direction, end_path=False):
        """Generate the G code and preview for the tool liftoff sequence.
        By default this is just a matter of bring the tool up to a safe height.
        :tool_direction: current tool angle.
        :start_path: True if the current position is at the end of a path
        """
        self.gc.tool_up()
            
    def sort_cutpaths(self, cutpath_list):
        """Sort the cutpaths to minimize tool movements.
        :cutpath_list: A list of tuples of the form
        (<path-id>, <path-geometry)
        """
        # TODO: use a better sort method...
        # Just sort the paths from bottom to top, left to right.
        # Only the first point of the path is used as a sort key...
        cutpath_list.sort(key=lambda cp: cp[1][0].p1.x)
        cutpath_list.sort(key=lambda cp: cp[1][0].p1.y)
        return cutpath_list
            
    def export(self, filename, directory='.', append_suffix=False):
        """Export the generated g code to a specified file.
        :filename: Base file name.
        :directory: Directory in which to create the file.
        :append_suffix: Append an auto-incrementing numeric suffix to the
        file name if True. Default is False.
        """
        filedir = os.path.expanduser(directory)
        filename = os.path.basename(filename)
        path = os.path.abspath(os.path.join(filedir, filename))
        if append_suffix:
            file_root, file_ext = os.path.splitext(filename)
            # Get a list of existing files that match the numeric suffix.
            # They should already be sorted.
            filter_str = '%s_[0-9]*%s' % (file_root, file_ext)
            files = fnmatch.filter(os.listdir(filedir), filter_str)
            if len(files) > 0:
                # Get the suffix from the last one and add one to it.
                # This seems overly complicated but it takes care of the case
                # where the user deletes a file in the middle of the
                # sequence which guarantees the newest file is always last.
                last_file = files[-1]
                file_root, file_ext = os.path.splitext(last_file)
                try:
                    suffix = int(file_root[-4:]) + 1
                except Exception:
                    suffix = 0
                filename = file_root[:-4] + ('%04d' % suffix) + file_ext
                path = os.path.join(filedir, filename)
            else:
                path = os.path.join(filedir, file_root + '_0000' + file_ext)
        with open(path, 'w') as f:
            f.write(str(self.gc))

    def _calc_rotation(self, current_angle, new_angle):
        """Calculate the relative rotation amount in radians."""
        # Normalize the angles to 0-360
        prev_angle = math.fmod(current_angle, 2*math.pi)
        if prev_angle < 0:
            prev_angle += 2*math.pi
        new_angle = math.fmod(new_angle + 2*math.pi, 2*math.pi)
        rotation_angle = new_angle - prev_angle
        if prev_angle < math.pi and new_angle > (prev_angle + math.pi):
            rotation_angle -= 2*math.pi
        elif prev_angle > math.pi and new_angle < (prev_angle - math.pi):
            rotation_angle += 2*math.pi
        return rotation_angle
        
    def _generate_segment_gcode(self, segment, depth):
        """Generate G code for Line and Arc path segments."""
        # If enabled and the tool has traveled a specified interval distance
        # then allow a subclass to generate some G code at this point.
        if self.feed_interval > 0.0 \
        and self.feed_interval_travel >= self.feed_interval:
            self.generate_feed_interval_gcode(self.current_angle, depth)
            self.preview_plotter.draw_interval_marker(self.last_point, depth)
            self.feed_interval_travel = 0.0
            
        seglen = segment.length()
        self.feed_distance += seglen
        self.feed_interval_travel += seglen
        previous_angle = self.current_angle
        self.current_angle += self._calc_rotation(self.current_angle,
                                                  segment.start_tangent())
        if seglen > geom.EPSILON:
            if isinstance(segment, geom.Line):
                self.gc.feed_rotate(self.current_angle)
                self.preview_plotter.draw_feed_rotate(self.last_point,
                                                      previous_angle,
                                                      self.current_angle,
                                                      self.tool_width, depth)
                self.gc.feed(segment.p2.x, segment.p2.y)
                self.preview_plotter.draw_feed_line(self.last_point,
                                                    segment.p2,
                                                    self.tool_width, depth)
            elif isinstance(segment, geom.Arc):
                arcv = segment.center - segment.p1
                self.gc.feed_arc(segment.is_clockwise(),
                                 segment.p2.x, segment.p2.y,
                                 arcv.x, arcv.y,
                                 self.current_angle + segment.angle)
                self.current_angle += segment.angle
                self.preview_plotter.draw_feed_arc(segment,
                                                   self.tool_width, depth)
        self.last_point = segment.p2
    
