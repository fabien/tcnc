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
import math
from datetime import datetime

EPSILON = 0.00001   # Default tolerance for floating point comparison

# Current machine target
TARGET = 'linuxcnc'
# Target machine info - machine name, version
TARGET_INFO = {'linuxcnc': ('LinuxCNC', '2.4'),}

class GCode(object):
    '''
    '''
    def __init__(self, zsafe, zfeed, xyfeed, afeed=None, tolerance=EPSILON, atolerance=EPSILON, use_degrees=True, wrap_angles=False, aaxis_name='A'):
        self.tolerance = tolerance        # Tolerance for float comparisons
        self.atolerance = atolerance    # Tolerance for angle comparisons
        self.use_degrees = use_degrees    # Use degrees instead of radians
        self.zsafe = zsafe                # Safe height of Z axis (for tool up)
        self.zfeed = zfeed                # Default Z axis feed rate
        self.xyfeed = xyfeed            # Default X,Y axis feed rate
        self.afeed = afeed if afeed is not None else xyfeed # Default rotational axis feed rate
        self.wrap_angles = wrap_angles    # Angles < 360 ?
        self.aaxis_name = aaxis_name    # Rotational axis name (A, B, or C)
        self.current_angle = 0.0        # Current angle of the rotational axis in radians
        self.last_val = {'X': None, 'Y': None, 'Z': None, 'A': None, 'I': None, 'J': None, 'R': None, 'F': None}
        self.scale = {'X': 1.0, 'Y': 1.0, 'Z': 1.0, 'A': 1.0, 'I': 1.0, 'J': 1.0, 'R': 1.0, 'F': 1.0}
        self.unit_scale = 1.0
        self.offset = {'X': 0.0, 'Y': 0.0, 'Z': 0.0, 'A': 0.0, 'I': 0.0, 'J': 0.0, 'R': 0.0, 'F': 0.0}
        self.gcode = ''
        self.is_tool_down = False
        self.verbose = True
        
    def __str__(self):
        return self.gcode
    
    def set_unit_scale(self, unit_scale):
        '''Set document unit scale factor for X,Y, and Z axis.
        Default is 1.0.
        '''
        self.unit_scale = unit_scale
    
    def addline(self, s='', ilc=None):
        '''Add a line to the G code output.
        This will append a newline character.
        :lc: an optional inline comment.
        '''
        if self.verbose and ilc:
            self.gcode += s + ' (' + ilc + ')\n'
        else:
            self.gcode += s + '\n'
        
    def comment(self, s):
        '''Add a G code comment.'''
        if self.verbose:
            self.addline('(%s)' % s)
            
    def default_header(self, units='in', description=None):
        '''Output a standard default G code file header.
        :units: 'mm' or 'in'. Default is 'in' (inches)
        :description: a list of comment lines
        '''
        self.addline('%')
        self.comment('----------------------------------------------------------')
        for line in description:
            self.comment(line)
        self.comment('Creation date: %s' % str(datetime.today()))
        self.comment('Target machine: %s, version %s' % \
                     (TARGET_INFO[TARGET][0], TARGET_INFO[TARGET][1]))
        self.comment('----------------------------------------------------------')
        if units == 'mm':
            self.addline('G21', 'Units are in millimeters')
        else:
            self.addline('G20', 'Units are in inches')
        self.addline('G90', 'Absolute positioning')
        self.addline('G40', 'Turn off tool diameter compensation')
    
    def default_footer(self):
        '''Output a default G code file footer.'''
        self.addline('%')
    
    def dwell(self, milliseconds):
        '''Output a dwell command which pauses the tool for the specified
        number of milliseconds'''
        # LinuxCNC interprets P as seconds whereas pretty much everything else
        # (ie Fanuc) interprets the parameter as milliseconds...
        # Also, LinuxCNC does not accept X or U as the parameter
        if milliseconds > 0:
            if TARGET == 'linuxcnc':
                seconds = milliseconds / 1000
                self.addline('G04 P%.4f' % seconds,
                             'Pause tool for %.4f seconds' % seconds)
            else:
                self.addline('G04 P%f' % milliseconds,
                             'Pause tool for %f milliseconds' % milliseconds)
    
    def tool_up(self):
        '''Moves tool to a safe Z axis height.'''
        self._gcode_line('G00', Z=self.zsafe, ilc='Tool up')
        self.is_tool_down = False
    
    def tool_down(self, z):
        '''Moves tool on Z axis to specified cutting height.'''
        if not self.is_tool_down:
            self._gcode_line('G01', Z=z, F=self.zfeed, ilc='Tool down')
            self.is_tool_down = True
            
    def rehome_rotational_axis(self):
        '''This will re-home the rotational axis to zero with minimal angular movement (<180).
        Useful when drawing large spirals with a tangent knife to minimize unwinding.
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
        self._gcode_line('G00', A=a, ilc='Rotate to home position')
        self._gcode_line('G92', A=0.0, ilc='Set axis origin to zero')
        # This might be EMC2-specific:
        self._gcode_line('G92.1', ilc='Reset variables 5211-5219 to zero')
    
    def rapid_move(self, x, y, z=None, a=None):
        '''Perform a rapid move to the specified location.'''
        if not self.is_tool_down:
            self._gcode_line('G00', X=x, Y=y, Z=z, A=a, ilc='Rapid move')
    
    def feed(self, x, y, z, feed=None):
        '''Perform a linear tool feed (cut) to the specified location.'''
        self._gcode_line('G01', X=x, Y=y, Z=z, F=(feed if feed is not None else self.xyfeed))
    
    def feed_arc(self, clockwise, x, y, z, arc_x, arc_y, a, feed=None):
        '''Perform an arc feed (cut).'''
        gc = 'G02' if clockwise else 'G03'
        self._gcode_line(gc, X=x, Y=y, Z=z, I=arc_x, J=arc_y, A=a, F=(feed if feed is not None else self.xyfeed))
        
    def feed_rotate(self, angle, feed=None):
        '''Rotate the A axis (around the Z axis) to the specified angle.'''
        self._gcode_line('G01', A=angle, F=(feed if feed is not None else self.afeed))
    
    def _gcode_line(self, gc, **kwargs):
        '''Output a line of gcode for G00, G01, G02, G03 commands.
        :ilc: an optional inline comment for verbose mode
        '''
        gcparams = ''
        for k in ('X', 'Y', 'Z', 'I', 'J', 'R', 'A', 'F'):
            v = kwargs.get(k)
            if v is not None:
                if k == 'A':
                    tolerance = self.atolerance        # Use angle tolerance for comparing angle values
                    if self.wrap_angles:
                        v = math.fmod(v, 2*math.pi)
                else:
                    tolerance = self.tolerance        # Otherwise use float tolerance
                if gc not in ('G01','G02','G03') or self.last_val[k] is None or abs(v - self.last_val[k]) > tolerance:
                    self.last_val[k] = v
                    if k == 'A':
                        self.current_angle = v
                        if self.use_degrees:
                            v = math.degrees(v)
                    elif k in 'XYZIJ':
                        # Apply unit scale (ie mm or inch)
                        v *= self.unit_scale
                    v = v * self.scale[k] + self.offset[k]    # Apply any axis transforms
                    k = self.aaxis_name if k == 'A' else k    # Rename rotational axis if specified
                    gcparams += ' %s%.5f' % (k, v)
        if len(gcparams) > 0:
            self.addline('%s%s' % (gc, gcparams), ilc=kwargs.get('ilc'))
        elif gc not in ('G01','G02','G03'):
            self.addline(gc, ilc=kwargs.get('ilc'))
    
