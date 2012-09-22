'''A G-code generator class that is suitable for a four axis
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
'''
import os
import sys
import fnmatch
import math
from datetime import datetime
import logging
logger = logging.getLogger(__name__)

EPSILON = 0.00001   # Default tolerance for floating point comparison

# Current machine target
TARGET = 'linuxcnc'
# Target machine info - machine name, version
TARGET_INFO = {'linuxcnc': ('LinuxCNC', '2.4+'),}

_MAXFLOAT = sys.float_info.max

DEFAULT_TRAJECTORY_TOLERANCE = 0.0

class GCode(object):
    '''
    Note: Angles are always specified in radians but output as degrees.
    '''
    # Order in which G code parameters are presented
    _gcode_ordered_params = 'XYZUVWABCIJKRDHLPQSF'
    
    # Non-modal G codes (LinuxCNC.)
    _gcode_nonmodal_group = ('G4', 'G04', 'G10', 'G28','G30', 'G53', 'G92')
    
    def __init__(self, xyfeed, zsafe=10.0, zfeed=None, zwait=0.0, afeed=None,
                 tolerance=EPSILON, atolerance=EPSILON,
                 wrap_angles=False, aaxis_name='A'):
        # Tolerance for float comparisons
        self.tolerance = tolerance
        # Tolerance for angle comparisons
        self.atolerance = atolerance
        # Trajectory plan tolerance (G64 P value). -1 ==> Exact path mode
        self.trajectory_tolerance = DEFAULT_TRAJECTORY_TOLERANCE
        # Safe height of Z axis (for tool up)
        self.zsafe = zsafe
        # Default Z axis feed rate
        self.zfeed = zfeed if zfeed is not None else xyfeed
        # Delay time in millis for tool-down
        self.zwait = zwait
        # Default X,Y axis feed rate
        self.xyfeed = xyfeed
        # Default rotational axis feed rate
        self.afeed = afeed if afeed is not None else xyfeed
        # Angles < 360 ?
        self.wrap_angles = wrap_angles
        # Rotational axis name (A, B, or C). Default is 'A'
        self.aaxis_name = aaxis_name
        # Show comments if True
        self.show_comments = True
        # Show line numbers if True
        self.show_line_numbers = False
        
        # Current angle of the rotational axis in radians
        self.current_angle = 0.0
        self.axis_scale = {'X': 1.0, 'Y': 1.0, 'Z': 1.0, 'A': 1.0}
        self.axis_offset = {'X': 0.0, 'Y': 0.0, 'Z': 0.0, 'A': 0.0}
        # User to machine unit scale
        self.unit_scale = 1.0
        self.is_tool_down = False
        self.line_number = 0
        self.last_val = {}
        for param in self._gcode_ordered_params:
            self.last_val[param] = sys.float_info.max
        # G code output buffer. See _add() method for implementation details
        self._gcodebuf = None
        
    def set_axis_offset(self, axis, offset):
        '''Set the offset for the specified axis.'''
        self.axis_offset[axis] = offset
            
    def set_axis_scale(self, axis, scale):
        '''Set the scale for the specified axis.'''
        self.axis_scale[axis] = scale
    
    def set_unit_scale(self, unit_scale):
        '''Set document unit scale factor for X,Y, and Z axis.
        Default is 1.0.
        '''
        self.unit_scale = unit_scale
    
    def show_line_numbers(self, show):
        '''Turn line numbers on or off.'''
        self.show_line_numbers = show
    
    def show_comments(self, show):
        '''Turn comments on or off.'''
        self.show_comments = show
    
    def addline(self, s='', comment=None):
        '''Add a line to the G code output.
        This will append a newline character.
        :comment: an optional inline comment.
        '''
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
        '''Add a G code comment line.'''
        if self.show_comments:
            self._add('(' + s + ')\n')
            
    def default_header(self, units='in', description=None):
        '''Output a standard default G code file header.
        :units: 'mm' or 'in'. Default is 'in' (inches)
        :description: a comment or a list of comments
        '''
        self.addline('%')
        self.comment('--------------------------------------------------------')
        if isinstance(description, tuple) or isinstance(description, list):
            for line in description:
                self.comment(line)
        else:
            self.comment(description)
        self.comment('Creation date: %s' % str(datetime.today()))
        self.comment('Target machine: %s, version %s' % \
                     (TARGET_INFO[TARGET][0], TARGET_INFO[TARGET][1]))
        self.comment('--------------------------------------------------------')
        self.addline('G17', 'XY plane')
        if units == 'mm':
            self.addline('G21', 'Units are in millimeters')
        else:
            self.addline('G20', 'Units are in inches')
        self.addline('G90', 'Use absolute positioning')
        self.addline('G40', 'Cancel tool diameter compensation')
        self.addline('G49', 'Cancel tool length compensation')
        if self.trajectory_tolerance < 0.0:
            self.addline('G61', 'Trajectory plan: exact path')
        elif self.trajectory_tolerance > 0.0:
            self.addline('G64 P%.4f' % self.trajectory_tolerance,
                         'Trajectory plan: blend with tolerance')
        else:
            self.addline('G64', 'Trajectory plan: blend without tolerance')
    
    def default_footer(self):
        '''Output a default G code file footer.'''
        self.addline()
        self.addline('M2', 'End program.')
        self.addline('%')
        
    def home(self, axes='XY', absolute=True):
        '''Rapid move to home position.
        Uses absolute home location, ignoring touch-off offsets.
        '''
        if self.is_tool_down:
            self.tool_up()
        gcode = 'G53 G0' if absolute else 'G0'
        for axis in axes:
            gcode += ' %s0.0' % axis.upper()
        self.addline(gcode, comment='Home')
    
    def pause(self, conditional=False):
        '''Pause the G code interpreter.
        Normally, pressing the start button in LinuxCNC/Axis
        will restart the interpreter.
        :conditional: use conditional stop if True.
        '''
        mcode = 'M1' if conditional else 'M0'
        self._gcode_line(mcode, comment='Pause')
    
    def dwell(self, milliseconds):
        '''Output a dwell command which pauses the tool for the specified
        number of milliseconds'''
        # LinuxCNC interprets P as seconds whereas pretty much everything else
        # (ie Fanuc) interprets the parameter as milliseconds...
        if milliseconds > 0:
            if TARGET == 'linuxcnc':
                seconds = milliseconds / 1000.0
                self.addline('G04 P%.4f' % seconds,
                             'Pause tool for %.4f seconds' % seconds)
            else:
                self.addline('G04 P%d' % milliseconds,
                             'Pause tool for %d milliseconds' % milliseconds)
    
    def tool_up(self, comment='Tool up'):
        '''Moves tool to a safe Z axis height.'''
        self._gcode_line('G00', Z=self.zsafe, comment=comment)
        self.is_tool_down = False
    
    def tool_down(self, z, wait=None, comment='Tool down'):
        '''Moves tool on Z axis to specified cutting height.
        :wait: the number of milliseconds to wait for the
        tool to actually get to the specified depth. This is useful for
        pneumatically controlled up/down axes where the actuator may take
        longer than the controller expects.'''
        if not self.is_tool_down:
            self._gcode_line('G01', Z=z, F=self.zfeed, comment=comment)
            if wait is None:
                wait = self.zwait
            if wait > 0:
                self.dwell(wait)
            self.is_tool_down = True
            
    def rehome_rotational_axis(self):
        '''This will re-home the rotational axis to zero with minimal
        angular movement (<180).
        Useful when drawing large spirals with a tangent knife to minimize
        unwinding.
        '''
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
        '''Perform a rapid move to the specified location.'''
        if not self.is_tool_down:
            self._gcode_line('G00', X=x, Y=y, Z=z, A=a, comment=comment)
    
    def rapid_rotate(self, a, comment=None):
        '''Perform a rapid rotation to the specified angle (in radians).'''
        if not self.is_tool_down:
            self._gcode_line('G00', A=a, comment=comment)
    
    def feed(self, x, y, z=None, a=None, feed=None, comment=None):
        '''Perform a linear tool feed (cut) to the specified location.'''
        self._gcode_line('G01', X=x, Y=y, Z=z, A=a,
                         F=(feed if feed is not None else self.xyfeed),
                         comment=comment)
    
    def feed_arc(self, clockwise, x, y, z, arc_x, arc_y, a,
                 feed=None, comment=None):
        '''Perform an arc feed (cut). Angle specified in radians.'''
        gc = 'G02' if clockwise else 'G03'
        self._gcode_line(gc, X=x, Y=y, Z=z, I=arc_x, J=arc_y, A=a,
                         F=(feed if feed is not None else self.xyfeed),
                         comment=comment)
        
    def feed_rotate(self, angle, feed=None, comment=None):
        '''Rotate the A axis (around the Z axis) to the specified angle
        (in radians).
        '''
        self._gcode_line('G01', A=angle,
                         F=(feed if feed is not None else self.afeed),
                         comment=comment)
    
    def export(self, directory, filename, append_suffix=False):
        '''Export the generated g code to a specified file.
        :directory: Directory in which to create the file.
        :filename: Base file name.
        :append_suffix: Append an auto-incrementing numeric suffix to the
        file name if True. Default is False.
        '''
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
            f.write(str(self))
            
    def _gcode_line(self, gcode, **kwargs):
        '''Output a line of gcode for G00, G01, G02, G03 commands.
        :comment: an optional inline comment for show_comments mode
        '''
        gcode_params = dict()
        for k in ('X', 'Y', 'Z', 'I', 'J', 'R', 'A', 'F'):
            v = kwargs.get(k)
            if v is not None:
                if k == 'A':
                    # Use angle tolerance for comparing angle values
                    tolerance = self.atolerance
                    if self.wrap_angles:
                        v = math.fmod(v, 2*math.pi)
                else:
                    # Otherwise use float tolerance
                    tolerance = self.tolerance
                value_has_changed = abs(v - self.last_val[k]) > tolerance
                gcode_is_nonmodal = gcode in self._gcode_nonmodal_group
                if value_has_changed or gcode_is_nonmodal:
                    self.last_val[k] = v
                    if k == 'A':
                        self.current_angle = v
                        v = math.degrees(v)
                    elif k in 'XYZIJ':
                        # Extra height check for safety
                        if k == 'Z':
                            self.is_tool_down = (v < self.zsafe)
                        # Apply unit scale (ie mm or inch)
                        v *= self.unit_scale
                    # Apply any axis transforms
                    v *= self.axis_scale.get(k, 1.0)
                    v += self.axis_offset.get(k, 0.0)
                    # Rename rotational axis if specified
                    k = self.aaxis_name if k == 'A' else k
                    gcode_params[k] = v
        if len(gcode_params) > 0:
            # Suppress feedrate-only lines
            if len(gcode_params.keys()) > 1 or 'F' not in gcode_params:
                gcode_line = '' + gcode
                for k in self._gcode_ordered_params:
                    v = gcode_params.get(k)
                    if v is not None:
                        gcode_line += ' %s%.5f' % (k, v)
                self.addline(gcode_line, comment=kwargs.get('comment'))
        elif gcode not in ('G01','G02','G03'):
            # Single G code with no parameters
            self.addline(gcode, comment=kwargs.get('comment'))
            
    def _add(self, s):
        '''Append the string to the gcode output buffer.'''
        if self._gcodebuf is None:
            self._gcodebuf = []
        self._gcodebuf.append(s)

    def __str__(self):
        return ''.join(self._gcodebuf)
    
