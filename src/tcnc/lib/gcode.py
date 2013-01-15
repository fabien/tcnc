"""A G-code generator class that is suitable for a four axis
machine with X, Y, and Z axis along with angular A
axis that rotates around the Z axis. It is more general but
that's the machine I have and the code might reflect that.

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

import sys
import math
from datetime import datetime
import logging

logger = logging.getLogger(__name__)

# Default tolerance for floating point comparison
EPSILON = 0.00001

# Current machine target
TARGET = 'linuxcnc'
# Target machine info - machine name, version
TARGET_INFO = {'linuxcnc': ('LinuxCNC', '2.4+'),}

DEFAULT_TRAJECTORY_MODE = 'G64'
DEFAULT_TRAJECTORY_TOLERANCE = 0.0
DEFAULT_FEED_RATE = 10.0

CLOCKWISE = CW = 1
COUNTERCLOCKWISE = CCW = 0

class GCode(object):
    """GCode generation class that describes a basic two axis (XY),
    three axis (XYZ), or four axis (XYZA)
    machine. The G code output is compatible with LinuxCNC.
    Note: Angles are always specified in radians but output as degrees.
    Axis values are always specified in user/world coordinates and output
    as machine units (ie inches or millimeters).
    Set GCode.unit_scale to the appropriate unit scale/conversion value.
    """
    
    # Private constants
    # Order in which G code parameters are written to output file
    _GCODE_ORDERED_PARAMS = 'XYZUVWABCIJKRDHLPQSF'
    # Non-modal G codes (LinuxCNC.)
    _GCODE_NONMODAL_GROUP = ('G4', 'G04', 'G10', 'G28','G30', 'G53', 'G92')
    
    header_comments = None
    # Tolerance for float comparisons
    tolerance = EPSILON
    # Tolerance for angle comparisons
    atolerance = EPSILON
    # Default X,Y axis feed rate
    xyfeed = DEFAULT_FEED_RATE
    # Safe height of Z axis (for tool up)
    zsafe = 10.0
    # Default Z axis feed rate
    zfeed = DEFAULT_FEED_RATE
    # Delay time in millis for tool-down
    zwait = 0.0
    # Default rotational axis feed rate
    afeed = DEFAULT_FEED_RATE
    # Angles < 360 ?
    wrap_angles = False
    # Show comments if True
    show_comments = True
    # Show line numbers if True
    show_line_numbers = False
    # Trajectory planning mode and blend tolerance
    trajectory_mode = DEFAULT_TRAJECTORY_MODE
    trajectory_tolerance = DEFAULT_TRAJECTORY_TOLERANCE
    # Default spindle speed
    spindle_speed = 0
    # Axis scale factors
    axis_scale = {'X': 1.0, 'Y': 1.0, 'Z': 1.0, 'A': 1.0}
    # Axis offsets
    axis_offset = {'X': 0.0, 'Y': 0.0, 'Z': 0.0, 'A': 0.0}
    # Units can be 'in' or 'mm'
    units = 'in'
    # User to machine unit scale
    unit_scale = 1.0
    
    def __init__(self, xyfeed, zsafe=10.0, zfeed=None, afeed=None,
                 tolerance=EPSILON, atolerance=EPSILON):
        self.zsafe = zsafe
        self.zfeed = zfeed if zfeed is not None else xyfeed
        self.xyfeed = xyfeed
        self.afeed = afeed if afeed is not None else xyfeed
        # Current angle of the rotational axis in radians
        self.current_angle = 0.0
        # True if the tool is on the cutting surface
        self.is_tool_down = False
        # Current line number
        self.line_number = 0
        # Last value for G code parameters
        self.last_val = {}
        for param in self._GCODE_ORDERED_PARAMS:
            self.last_val[param] = sys.float_info.max
        # G code output buffer. See _add() method for implementation details
        self._gcodebuf = None
        
    def set_axis_offsets(self, axes, offsets):
        """Set the offset for the specified axes.
        Axis offsets are always specified in machine units.
        Angular offsets are always in degrees.
        """
        for axis, offset in zip(axes, offsets):
            self.axis_offset[axis] = offset
            
    def set_axis_scales(self, axes, scales):
        """Set the scaling factors for the specified axes."""
        for axis, scale in zip(axes, scales):
            self.axis_scale[axis] = scale
    
    def addline(self, s='', comment=None):
        """Add a line to the G code output.
        A newline character is always appended, even if the string is empty.
        :s: The string to add (optional)
        :comment: An inline comment (optional)
        """
        self.line_number += 1
        linebuf = ''
        if self.show_line_numbers:
            linebuf += 'N%d' % self.line_number
            if s:
                linebuf += ' '
        linebuf += s
        if self.show_comments and comment:
            if s:
                linebuf += ' '
            linebuf += '(' + comment + ')'
        linebuf += '\n'
        self._add(linebuf)
        
    def comment(self, s):
        """Add a G code comment line."""
        if self.show_comments:
            self._add('(' + s + ')\n')
            
    def default_header(self, comment=None):
        """Output a standard default G code file header.
        :comment: A header comment or a list of comments (optional).
        """
        self.addline('%')
        self.comment('--------------------------------------------------------')
        if comment is None:
            comment = self.header_comments
        if comment:
            if isinstance(comment, tuple) \
            or isinstance(comment, list):
                for comment_line in comment:
                    self.comment(comment_line)
            else:
                self.comment(comment)
        self.comment('Creation date: %s' % datetime.today().isoformat(' '))
        self.comment('Target machine: %s, version %s' % \
                     (TARGET_INFO[TARGET][0], TARGET_INFO[TARGET][1]))
        self.comment('--------------------------------------------------------')
        self.addline('G17', 'XY plane')
        if self.units == 'mm':
            self.addline('G21', 'Units are in millimeters')
        else:
            self.addline('G20', 'Units are in inches')
        self.addline('G90', 'Use absolute positioning')
        self.addline('G40', 'Cancel tool diameter compensation')
        self.addline('G49', 'Cancel tool length compensation')
        if self.trajectory_mode == 'G64P':
            self.addline('G64 P%.4f' % self.trajectory_tolerance,
                         'Trajectory plan: blend with tolerance')
        elif self.trajectory_mode:
            self.addline(self.trajectory_mode, 'Trajectory planning mode')
    
    def default_footer(self):
        """Output a default G code file footer."""
        self.addline()
        self.addline('M2', 'End program.')
        self.addline('%')
        
    def home(self, axes='XY', absolute=True):
        """Rapid move to home position.
        Uses absolute home location, ignoring touch-off offsets.
        :axes: Which axes to home (default is 'XY')
        :absolute: Use absolute G53 coordinates (default is True)
        """
        if self.is_tool_down:
            self.tool_up()
        gcode = 'G53 G0' if absolute else 'G0'
        for axis in axes:
            gcode += ' %s0.0' % axis.upper()
        self.addline(gcode, comment='Home')
    
    def pause(self, conditional=False, comment='Pause'):
        """Pause the G code interpreter.
        Normally, pressing the start button in LinuxCNC/Axis
        will restart the interpreter.
        :conditional: use conditional stop if True.
        """
        mcode = 'M1' if conditional else 'M0'
        self._gcode_line(mcode, comment=comment)
    
    def dwell(self, milliseconds, comment=None):
        """Output a dwell command which pauses the tool for the specified
        number of milliseconds.
        """
        # LinuxCNC interprets P as seconds whereas pretty much everything else
        # (ie Fanuc) interprets the parameter as milliseconds...
        if milliseconds > 0:
            if TARGET == 'linuxcnc':
                seconds = milliseconds / 1000.0
                if comment is None:
                    comment = 'Pause tool for %.4f seconds' % seconds
                self.addline('G04 P%.4f' % seconds, comment=comment)
            else:
                if comment is None:
                    comment = 'Pause tool for %d milliseconds' % milliseconds
                self.addline('G04 P%d' % milliseconds, comment=comment)
    
    def tool_up(self, comment='Tool up'):
        """Moves tool to a safe Z axis height."""
        self._gcode_line('G00', Z=self.zsafe, force_value='Z', comment=comment)
        self.is_tool_down = False
    
    def tool_down(self, z, wait=None, comment='Tool down'):
        """Moves tool on Z axis to specified cutting height.
        :wait: the number of milliseconds to wait for the
        tool to actually get to the specified depth.
        This is mainly useful for
        pneumatically controlled up/down axes where the actuator may take
        a few milliseconds to extend."""
        if not self.is_tool_down:
            self._gcode_line('G01', Z=z, F=self.zfeed, comment=comment)
            if wait is None:
                wait = self.zwait
            if wait > 0:
                self.dwell(wait)
            self.is_tool_down = True
            
    def spindle_on(self, speed=None, direction=CLOCKWISE, comment='Spindle on'):
        """Turn the spindle on and set the spindle speed.
        If <speed> is None (default) the last specified speed will be used.
        """
        if speed is None:
            speed = self.spindle_speed
        else:
            self.spindle_speed = speed
        mcode = 'M3' if direction == CLOCKWISE else 'M4'
        self.addline('%s S%d' % (mcode, speed), comment=comment)
        
    def spindle_off(self, comment='Spindle off'):
        """Turn off the spindle."""
        self.addline('M5')
            
    def rehome_rotational_axis(self):
        """This will re-home the rotational axis to zero with minimal
        angular movement (<180).
        Useful when drawing large spirals with a tangent knife to minimize
        unwinding.
        """
        if self.is_tool_down:
            self.tool_up()
        self.comment('Rehome rotational axis to zero to minimize long unwinds')
        da = math.fmod(self.current_angle, 2*math.pi)
        if math.fabs(da) < math.pi:
            a = self.current_angle - da
        elif da < 0:
            a = self.current_angle - (2*math.pi + da)
        else:
            a = self.current_angle + (2*math.pi - da)
        self.current_angle = 0
        self._gcode_line('G00', A=a, comment='Rotate to home position')
        self._gcode_line('G92', A=0.0, comment='Set axis origin to zero')
        # This might be EMC2-specific:
        self._gcode_line('G92.1', comment='Reset variables 5211-5219 to zero')
    
    def rapid_move(self, x, y, z=None, a=None, comment=None):
        """Perform a rapid move to the specified location."""
        if not self.is_tool_down:
            self._gcode_line('G00', X=x, Y=y, Z=z, A=a, comment=comment)
    
    def rapid_rotate(self, a, comment=None):
        """Perform a rapid rotation to the specified angle (in radians)."""
        if not self.is_tool_down:
            self._gcode_line('G00', A=a, comment=comment)
    
    def feed(self, x, y, z=None, a=None, feed=None, comment=None):
        """Perform a linear tool feed (cut) to the specified location."""
        self._gcode_line('G01', X=x, Y=y, Z=z, A=a,
                         F=(feed if feed is not None else self.xyfeed),
                         comment=comment)
    
    def feed_arc(self, clockwise, x, y, arc_x, arc_y, a, z=None,
                 feed=None, comment=None):
        """Perform an arc feed (cut). Angle specified in radians."""
        gc = 'G02' if clockwise else 'G03'
        self._gcode_line(gc, X=x, Y=y, Z=z, I=arc_x, J=arc_y, A=a,
                         F=(feed if feed is not None else self.xyfeed),
                         comment=comment)
        
    def feed_rotate(self, angle, feed=None, comment=None):
        """Rotate the A axis (around the Z axis) to the specified angle
        (in radians).
        """
        self._gcode_line('G01', A=angle,
                         F=(feed if feed is not None else self.afeed),
                         comment=comment)
    
    def _gcode_line(self, gcode, **kwargs):
        """Output a line of gcode for G00, G01, G02, G03 commands.
        :comment: an optional inline comment for show_comments mode
        :force_value: a string containing the modal parameter names whose
        values will be output regardless of whether their values have changed
        """
        gcode_params = dict()
        force_value = kwargs.get('force_value', '')
        for k in ('X', 'Y', 'Z', 'I', 'J', 'R', 'A', 'F'):
            v = kwargs.get(k)
            if v is not None:
                if k in 'ABC':
                    # Use angle tolerance for comparing angle values
                    tolerance = self.atolerance
                    if self.wrap_angles:
                        v = math.fmod(v, 2*math.pi)
                else:
                    # Otherwise use float tolerance
                    tolerance = self.tolerance
                value_has_changed = abs(v - self.last_val[k]) > tolerance
                gcode_is_nonmodal = gcode in self._GCODE_NONMODAL_GROUP
                if k in force_value or value_has_changed or gcode_is_nonmodal:
                    self.last_val[k] = v
                    # Apply any axis transforms 
                    v *= self.axis_scale.get(k, 1.0)
                    v += self.axis_offset.get(k, 0.0)
                    if k in 'ABC':
                        self.current_angle = v
                        v = math.degrees(v)
                    elif k in 'XYZIJ':
                        # Extra height check for safety
                        if k == 'Z':
                            self.is_tool_down = (v < self.zsafe)
                        # Apply unit scale (user/world to machine units)
                        v *= self.unit_scale
                    gcode_params[k] = v
        if len(gcode_params) > 0:
            # Suppress feedrate-only lines
            if len(gcode_params.keys()) > 1 or 'F' not in gcode_params:
                gcode_line = '' + gcode
                for k in self._GCODE_ORDERED_PARAMS:
                    v = gcode_params.get(k)
                    if v is not None:
                        gcode_line += ' %s%.5f' % (k, v)
                self.addline(gcode_line, comment=kwargs.get('comment'))
        elif gcode not in ('G01','G02','G03'):
            # Single G code with no parameters
            self.addline(gcode, comment=kwargs.get('comment'))
            
    def _add(self, s):
        """Append the string to the gcode output buffer."""
        if self._gcodebuf is None:
            self._gcodebuf = []
        self._gcodebuf.append(s)

    def __str__(self):
        return ''.join(self._gcodebuf)
    
