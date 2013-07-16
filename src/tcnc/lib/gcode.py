"""A G-code generator class that is suitable for a four axis
machine with X, Y, and Z axis along with angular A
axis that rotates around the Z axis. It is much more general but
that's the machine I have and the code might reflect that.

Copyright (C) 2012, 2013 Claude Zervas, claude@utlco.com

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
from datetime import datetime

logger = logging.getLogger(__name__)

_DEFAULT_TOLERANCE = 1e-9
"""Default tolerance for floating point comparison"""


class PreviewPlotter(object):
    """Default base class that can be subclassed by users of the GCode
    class to provide a graphical preview of the G-code output."""
    def plot_move(self, endp):
        """Plot G00 - rapid move from current tool location to to ``endp``.
        
        :param endp: Endpoint of move as a 4-tuple (x, y, z, a).
        """
        pass
    
    def plot_feed(self, endp):
        """Plot G01 - linear feed from current tool location to ``endp``.
        
        :param endp: Endpoint of feed as a 4-tuple (x, y, z, a).
        """
        pass
    
    def plot_arc(self, endp, clockwise, arc_x, arc_y):
        """Plot G02/G03 - arc feed from current tool location to to ``endp``.
        
        :param endp: Endpoint of feed as a 4-tuple (x, y, z, a).
        :param clockwise: True if the arc moves in a clockwise direction.
        :param arc_x: Center of arc - X coordinate relative to ``x``
        :param arc_y: Center of arc - Y coordinate relative to ``y``
        """
        pass
    

class GCode(object):
    """GCode generation class that describes a basic two axis (XY),
    three axis (XYZ), or four axis (XYZA)
    machine. The G code output is compatible with LinuxCNC.
    
    Angles are always specified in radians but output as degrees.
    
    Axis values are always specified in user/world coordinates and output
    as machine units (ie inches or millimeters) using ``GCode.unit_scale``
    as the scaling factor.
    """
    
    #: Current machine target
    TARGET = 'linuxcnc'
    #: Target machine info - machine name, version
    _TARGET_INFO = {'linuxcnc': ('LinuxCNC', '2.4+'),}
    
    # Private constants
    # Order in which G code parameters are specified in a line of G code
    _GCODE_ORDERED_PARAMS = 'XYZUVWABCIJKRDHLPQSF'
    # Non-modal G codes (LinuxCNC.)
    _GCODE_NONMODAL_GROUP = ('G4', 'G04', 'G10', 'G28','G30', 'G53', 'G92')
    # G codes where a feed rate is required
    _GCODE_FEED = ('G01', 'G02', 'G03')
    # G codes that change the position of the tool
    _GCODE_MOTION = ('G0', 'G00', 'G1', 'G01', 'G2', 'G02', 'G3', 'G03')
    # G codes that are suppressed if the parameters remain unchanged
    _GCODE_MODAL_MOTION = ('G00', 'G01', 'G02', 'G03')
        
    #: Values can be changed to accomodate machines that expect
    #: different axis names (ie. using C instead of A or UVW instead of XYZ)
    axis_map = {'X': 'X', 'Y': 'Y', 'Z': 'Z', 'A': 'A'}
    
    # :Comments added to header part
    header_comments = None
    # :Tolerance for float comparisons
    tolerance = _DEFAULT_TOLERANCE
    #: Tolerance for angle comparisons
    atolerance = _DEFAULT_TOLERANCE
    #: X,Y axis feed rate
    xyfeed = 0.0
    #: Safe height of Z axis (for tool up)
    zsafe = 10.0
    #: Z axis feed rate
    zfeed = 0.0
    #: Delay time in millis for tool-down
    zwait = 0.0
    #: Delay time in milliseconds for spindle on.
    wait_spindle_on = 0.0
    #: Delay time in milliseconds for spindle off
    wait_spindle_off = 0.0
    #: Rotational axis feed rate
    afeed = 0.0
    #: Angles < 360 ?
    wrap_angles = False
    #: Show comments if True
    show_comments = True
    #: Show line numbers if True
    show_line_numbers = False
    #: Trajectory planning mode and blend tolerance
    trajectory_mode = None
    trajectory_tolerance = 0
    #: Default spindle speed
    spindle_speed = 0
    #: Axis scale factors
    axis_scale = {'X': 1.0, 'Y': 1.0, 'Z': 1.0, 'A': 1.0}
    #: Axis offsets
    axis_offset = {'X': 0.0, 'Y': 0.0, 'Z': 0.0, 'A': 0.0}
    #: Units can be 'in' or 'mm'
    units = 'in'
    #: User to machine unit scale
    unit_scale = 1.0
    #: G code to output just before the tool down Z axis move
    tool_down_prefix = None
    #: G code to output just after the tool up Z axis move
    tool_up_suffix = None
    
    
    def __init__(self, xyfeed, zsafe, zfeed=None, afeed=None,
                 tolerance=_DEFAULT_TOLERANCE, atolerance=_DEFAULT_TOLERANCE, plotter=None):
        """
        :param xyfeed: Default feed rate along X and Y axes
        :param zsafe: The safe height of the Z axis for rapid XY moves
        :param zfeed: Default feed rate along Z axis (defaults to :xyfeed:)
        :param afeed: Default feed rate along A axis
        :param tolerance: The tolerance for floating point comparisons
        :param atolerance: The tolerance for angle comparisons
        :param plotter: Preview plotter
        """
        self.zsafe = zsafe
        self.zfeed = zfeed if zfeed is not None else xyfeed
        self.xyfeed = xyfeed
        self.afeed = afeed if afeed is not None else xyfeed
        # Current angle of the rotational axis in radians
        self.current_angle = 0.0
        # Current line number
        self.line_number = 0
        # The current preview plotter
        self.preview_plotter = plotter if plotter is not None else PreviewPlotter()
        # Last value for G code parameters
        self._last_val = {}
        for param in self._GCODE_ORDERED_PARAMS:
            self._last_val[param] = None
        # G code output buffer; a list of G code lines
        # See _add() method for implementation details
        self._gcodebuf = None
        # True if the tool is above the Z axis safe height for rapid moves
        self._is_tool_up = False
        
    def set_axis_offsets(self, axes, offsets):
        """Set the offset for the specified axes.
        
        Axis offsets are always specified in machine units.
        Angular offsets are always in degrees.
        
        :param axes: A string containing the axes (ie 'XYZ')
        """
        for axis, offset in zip(axes, offsets):
            self.axis_offset[axis] = offset
            
    def set_axis_scales(self, axes, scales):
        """Set the scaling factors for the specified axes.
        """
        for axis, scale in zip(axes, scales):
            self.axis_scale[axis] = scale
    
    def comment(self, s=None):
        """Add a G code comment line.
        
        Encloses the comment string in parentheses.
        Outputs a newline if the comment string is None (default).
        
        :param s: The comment string.
        """
        if s is None:
            self._add('\n')
        elif self.show_comments:
            self._add('(' + s + ')\n')
            
    def default_header(self, comment=None):
        """Output a pretty standard G code file header.
        
        :param comment: A header comment or a list of comments (optional).
        """
        self.comment('--------------------------------------------------------')
        if comment is None:
            comment = self.header_comments
        if comment is not None:
            if isinstance(comment, tuple) \
            or isinstance(comment, list):
                for comment_line in comment:
                    self.comment(comment_line)
            else:
                self.comment(comment)
        self.comment('Creation date: %s' % datetime.today().isoformat(' '))
        self.comment('Target machine: %s, version %s' % \
                     (self._TARGET_INFO[self.TARGET][0],
                      self._TARGET_INFO[self.TARGET][1]))
        self.comment('--------------------------------------------------------')
        self._addline('G17', 'XY plane')
        if self.units == 'mm':
            self._addline('G21', 'Units are in millimeters')
        else:
            self._addline('G20', 'Units are in inches')
        self._addline('G90', 'Use absolute positioning')
        self._addline('G40', 'Cancel tool diameter compensation')
        self._addline('G49', 'Cancel tool length compensation')
        if self.trajectory_mode == 'G64P':
            self._addline('G64 P%.4f' % self.trajectory_tolerance,
                         'Trajectory plan: blend with tolerance')
        elif self.trajectory_mode:
            self._addline(self.trajectory_mode, 'Trajectory planning mode')
        self.feed_rate(self.xyfeed)
    
    def default_footer(self):
        """Output a G code file footer."""
        self._add('\n')
        self._addline('M2', 'End program.')
        
    def feed_rate(self, feed_rate):
        """Output the specified feed rate (if changed).
        
        :param feed_rate: The feed rate in units (inches or mm) per minute.
        """
        if self._last_val['F'] is None \
        or abs(feed_rate - self._last_val['F']) > _DEFAULT_TOLERANCE:
            self._addline('F%.4f' % feed_rate)
            self._last_val['F'] = feed_rate
            
    def pause(self, conditional=False, comment='Pause'):
        """Pause the G code interpreter.
        Normally, pressing the start button in LinuxCNC/Axis
        will restart the interpreter.
        
        :param conditional: use conditional stop if True.
        :param comment: Optional comment string.
        """
        mcode = 'M1' if conditional else 'M0'
        self._gcode_line(mcode, comment=comment)
    
    def dwell(self, milliseconds, comment=None):
        """Output a dwell command which pauses the tool for the specified
        number of milliseconds.
        
        :param milliseconds: Number of milliseconds to pause.
        :param comment: Optional comment string.
        """
        # LinuxCNC interprets P as seconds whereas pretty much everything else
        # (ie Fanuc) interprets the parameter as milliseconds...
        if milliseconds > 0:
            if self.TARGET == 'linuxcnc':
                seconds = milliseconds / 1000.0
                if comment is None:
                    comment = 'Pause tool for %.4f seconds' % seconds
                self._addline('G04 P%.4f' % seconds, comment=comment)
            else:
                if comment is None:
                    comment = 'Pause tool for %d milliseconds' % milliseconds
                self._addline('G04 P%d' % milliseconds, comment=comment)
    
    def tool_up(self, wait=None, comment=None):
        """Moves tool to a safe Z axis height.
        
        :param wait: the number of milliseconds to wait for the tool to retract.\
        Uses GCode.zwait value by default if None specified.\
        This parameter is mainly useful for pneumatically\
        controlled up/down axes where the actuator may take\
        a few milliseconds to extend/retract.
        :param comment: Optional comment string.
        """
        # Note: self._is_tool_up is purposely not checked here to insure
        # that the tool is forced to a safe height regardless of internal state
        self._gcode_line('G00', Z=self.zsafe, force_value='Z', comment=comment)
        if self.tool_up_suffix:
            self._addline(self.tool_up_suffix)
        if wait is None:
            wait = self.zwait
        if wait > 0:
            self.dwell(wait)
        self._is_tool_up = True
    
    def tool_down(self, z, feed=None, wait=None, comment=None):
        """Moves tool on Z axis down to specified cutting height.
        
        :param z: Height of Z axis to move to.
        :param wait: the number of milliseconds to wait for the tool to\
        actually get to the specified depth. Uses GCode.zwait value by default\
        if None specified. This parameter is mainly useful for pneumatically\
        controlled up/down axes where the actuator may take\
        a few milliseconds to extend/retract.
        :param comment: Optional comment string.
        """
        if self._is_tool_up:
            if feed is None:
                feed = self.zfeed
            if self.tool_down_prefix:
                self._addline(self.tool_down_prefix)
            self._gcode_line('G01', Z=z, F=feed, comment=comment)
            self._is_tool_up = False
            if wait is None:
                wait = self.zwait
            if wait > 0:
                self.dwell(wait)
            
    def spindle_on(self, speed=None, clockwise=True, wait=None, comment='Spindle on'):
        """Turn on the spindle.
        
        :param speed: Spindle speed in RPM\
            (if None use last specified speed default).
        :param clockwise: Spindle direction (clockwise by default).
        :param wait: Number of milliseconds to wait for the spindle to reach\
        full speed. Uses ``GCode.wait_spindle_on`` value by default.
        :param comment: Optional comment string.
        """
        if speed is None:
            speed = self.spindle_speed
        else:
            self.spindle_speed = speed
        mcode = 'M3' if clockwise else 'M4'
        self._addline('%s S%d' % (mcode, speed), comment=comment)
        if wait is None:
            wait = self.wait_spindle_on
        self.dwell(wait)
        
    def spindle_off(self, wait=None, comment='Spindle off'):
        """Turn off the spindle.
        
        :param wait: the number of milliseconds to wait for the spindle\
        to stop. Uses ``GCode.wait_spindle_off`` value by default.
        :param comment: Optional comment string.
        """
        self._addline('M5')
        if wait is None:
            wait = self.wait_spindle_off
        self.dwell(wait)
            
    def rehome_rotational_axis(self):
        """This will re-home the rotational axis to zero with minimal
        angular movement (<180).
        
        Useful when drawing large spirals with a tangent knife to minimize
        unwinding.

        :param comment: Optional comment string.
        """
        if not self._is_tool_up:
            self.tool_up()
        self.comment('Rehome rotational axis to zero to minimize long unwinds')
        a = math.fmod(self.current_angle, math.pi)
        if abs(a) > self.atolerance:
            self.current_angle -= a
            self._gcode_line('G00', A=self.current_angle,
                             comment='Rotate to home position')
            self._gcode_line('G92', A=0.0, comment='Set axis origin to zero')
            # This might be EMC2-specific:
            self._gcode_line('G92.1', comment='Reset variables 5211-5219 to zero')
            self.current_angle = 0.0
    
    def rapid_move(self, x=None, y=None, z=None, a=None, comment=None):
        """Perform a rapid move to the specified location.

        At least one axis should be specified.
        If the Z axis is specified it must be greater than the safe move height
        ``GCode.zsafe``.

        :param x: X axis value (optional)
        :param y: Y axis value (optional)
        :param z: Z axis value (optional)
        :param a: A axis value (optional)
        :param comment: Optional comment string.
        """
        if not self._is_tool_up:
            self.tool_up()
        z = max(self.zsafe, z)
        self._gcode_line('G00', X=x, Y=y, Z=z, A=a, comment=comment)
        self.preview_plotter.plot_move(self._endp(x, y, z, a))
        
    def feed(self, x=None, y=None, z=None, a=None, feed=None, comment=None):
        """Perform a linear tool feed to the specified location.
        
        At least one axis should be specified.
        
        :param x: X axis value (optional)
        :param y: Y axis value (optional)
        :param z: Z axis value (optional)
        :param a: A axis value (optional)
        :param feed: Feed rate (optional - last used feed rate by default)
        :param comment: Optional comment string.
        """
        # Determine default feed rate appropriate for the move
        if feed is None:
            if x is not None or y is not None:
                feed = self.xyfeed
            elif z is not None:
                feed = self.zfeed
            elif a is not None:
                feed = self.afeed
            else:
                # No feed rate, no axis specified - nothing to feed
                return
        self._gcode_line('G01', X=x, Y=y, Z=z, A=a, F=feed, comment=comment)
        self.preview_plotter.plot_feed(self._endp(x, y, z, a))
    
    def feed_arc(self, clockwise, x, y, arc_x, arc_y, a=None, z=None,
                 feed=None, comment=None):
        """Perform an arc feed.

        :param x: X value of arc end point
        :param y: Y value of arc end point
        :param arc_x: Center of arc relative to ``x``
        :param arc_y: Center of arc relative to ``y``
        :param clockwise: True if the arc moves in a clockwise direction.
        :param a: A axis value at endpoint (in radians)
        :param feed: Feed rate (optional - last used XY feed rate by default)
        :param comment: Optional comment string.
        """
        gc = 'G02' if clockwise else 'G03'
        self._gcode_line(gc, X=x, Y=y, Z=z, I=arc_x, J=arc_y, A=a,
                         F=(feed if feed is not None else self.xyfeed),
                         force_value='IJ',
                         comment=comment)
        self.preview_plotter.plot_arc(self._endp(x, y, z, a),
                                      clockwise, arc_x, arc_y)
        
    def add_gcode(self, code, params='', comment=None):
        """Add a line of G code.
        
        The line will be numbered if the line numbering option is enabled.
        
        Not to used with codes such as G0, G1, G2, G3 that can change the
        position of the tool. Use appropriate move or feed methods instead.
        
        :param code: The G code to output
        :param params: Optional parameters
        :param comment: Optional comment string.
        """
        code = code.upper()
        if code in self._GCODE_MOTION:
            raise Exception('%s will cause position errors.' % code)
        self._addline(code + ' ' + params, comment)
        
    def __str__(self):
        """Return the complete G code file as a string."""
        return ''.join(('%\n', ''.join(self._gcodebuf), '%\n'))
    
    def _gcode_line(self, gcode, **kwargs):
        """Output a line of gcode for G00, G01, G02, G03 commands.

        :param force_value: a string containing the modal parameter names whose
        values will be output regardless of whether their values have changed.
        By default if the specified value of a modal parameter has not changed
        since its last value then it will not be output.
        :param comment: Optional comment string.
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
#                gcode_params[k] = self._process_modal_value(k, v)
                value_has_changed = self._last_val[k] is None or \
                                    abs(v - self._last_val[k]) > tolerance
                gcode_is_nonmodal = gcode in self._GCODE_NONMODAL_GROUP
                if k in force_value or value_has_changed or gcode_is_nonmodal:
                    self._last_val[k] = v
                    # Apply any axis transforms 
                    v *= self.axis_scale.get(k, 1.0)
                    v += self.axis_offset.get(k, 0.0)
                    if k in 'ABC':
                        self.current_angle = v
                        v = math.degrees(v)
                    elif k in 'XYZIJ':
                        # Extra height check for safety
                        if k == 'Z':
                            self._is_tool_up = (v >= self.zsafe)
                        # Apply unit scale (user/world to machine units)
                        v *= self.unit_scale
                    gcode_params[k] = v
#        if gcode not in self._GCODE_NONMODAL_GROUP:
#            self._trim_unchanged_modal_params(gcode, gcode_params, tolerance)
        if len(gcode_params) > 0:
            # Suppress feedrate-only lines
            if len(gcode_params.keys()) > 1 or 'F' not in gcode_params:
#                if gcode in self._GCODE_FEED:
#                    # Just in case a feed rate wasn't specified use a default
#                    f = gcode_params.get('F', self.xyfeed)
#                    if abs(f - self._last_val['F']) > self.tolerance:
#                        gcode_params['F'] = f
                gcode_line = '' + gcode
                for k in self._GCODE_ORDERED_PARAMS:
                    v = gcode_params.get(k)
                    if v is not None:
                        k = self.axis_map.get(k, k)
                        gcode_line += ' %s%.5f' % (k, v)
                self._addline(gcode_line, comment=kwargs.get('comment'))
        # Note: this check will suppress output of gcodes with unchanged params
        elif gcode not in self._GCODE_MODAL_MOTION:
            # Single G code with no parameters
            self._addline(gcode, comment=kwargs.get('comment'))
    
#    def _trim_unchanged_modal_params(self, gcode, params, tolerance):
#        for k in params:
#            v = params[k]
#            value_has_changed = abs(v - self._last_val[k]) > tolerance
#            if not value_has_changed:
#                # Make sure arc commands aren't left with nothing
#                del params[k]
#            else:
#                self._last_val[k] = v
        
#    def _process_modal_value(self, k, v):
##        self._last_val[k] = v
#        # Apply any axis transforms 
#        v *= self.axis_scale.get(k, 1.0)
#        v += self.axis_offset.get(k, 0.0)
#        if k in 'ABC':
#            self.current_angle = v
#            v = math.degrees(v)
#        elif k in 'XYZIJ':
#            # Extra height check for safety
#            if k == 'Z':
#                self.is_tool_down = (v < self.zsafe)
#            # Apply unit scale (user/world to machine units)
#            v *= self.unit_scale
#        return v        
             
    def _addline(self, s='', comment=None):
        """Add a numbered line to the G code output.
        A newline character is always appended, even if the string is empty.
        :s: The string to add (optional)
        :comment: An inline comment (optional)
        """
        self.line_number += 1
        linebuf = ''
        if self.show_line_numbers and s:
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
        
    def _add(self, s):
        """Append the string to the gcode output buffer."""
        if self._gcodebuf is None:
            self._gcodebuf = []
        self._gcodebuf.append(s)

    def _endp(self, x, y, z, a):
        """Return the end point of the current trajectory.
        Used for preview plotting."""
        return (x if x is not None else self._last_val['X'],
                y if y is not None else self._last_val['Y'],
                z if z is not None else self._last_val['Z'],
                a if a is not None else self._last_val['A'])

