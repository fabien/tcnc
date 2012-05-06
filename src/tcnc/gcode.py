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

EPSILON = 0.00001   # Default tolerance for floating point comparison

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
        self.offset = {'X': 0.0, 'Y': 0.0, 'Z': 0.0, 'A': 0.0, 'I': 0.0, 'J': 0.0, 'R': 0.0, 'F': 0.0}
        self.gcode = ''
        self.is_tool_down = False
        
    def __str__(self):
        return self.gcode
    
    def comment(self, s):
        self.gcode += '(%s)\n' % s
    
    def tool_up(self):
        self._gcode_line('G00', Z=self.zsafe)
        self.is_tool_down = False
    
    def tool_down(self, z):
        if not self.is_tool_down:
            self._gcode_line('G01', Z=z, F=self.zfeed)
            self.is_tool_down = True
            
    def rehome_rotational_axis(self):
        '''This will re-home the rotational axis to zero with minimal angular movement (<180).
        Useful when drawing large spirals with a tangent knife to minimize unwinding.
        '''
        self.comment('rehome rotational axis to zero')
        da = math.fmod(self.current_angle, 2*math.pi)
        if math.fabs(da) < math.pi:
            a = self.current_angle - da
        elif da < 0:
            a = self.current_angle - (2*math.pi + da)
        else:
            a = self.current_angle + (2*math.pi - da)
        self.current_angle = 0
        self._gcode_line('G01', A=a)    # Rotate to home position (0 deg)
        self._gcode_line('G92', A=0.0)    # Move 'A' axis origin to zero
        self._gcode_line('G92.1')        # Reset parameters 5211-5219 to zero (this might be EMC2-specific)
    
    def rapid_move(self, x, y, z=None, a=None):
        if not self.is_tool_down:
            self._gcode_line('G00', X=x, Y=y, Z=z, A=a)
    
    def feed(self, x, y, z, feed=None):
        self._gcode_line('G01', X=x, Y=y, Z=z, F=(feed if feed is not None else self.xyfeed))
    
    def feed_arc(self, clockwise, x, y, z, arc_x, arc_y, a, feed=None):
        gc = 'G02' if clockwise else 'G03'
        self._gcode_line(gc, X=x, Y=y, Z=z, I=arc_x, J=arc_y, A=a, F=(feed if feed is not None else self.xyfeed))
        
    def feed_rotate(self, angle, feed=None):
        self._gcode_line('G01', A=angle, F=(feed if feed is not None else self.afeed))
    
    def append_gcode(self, s):
        self.gcode += s + '\n'
        
    def _gcode_line(self, gc, **kwargs):
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
                    v = v * self.scale[k] + self.offset[k]    # Apply any axis transforms
                    k = self.aaxis_name if k == 'A' else k    # Rename rotational axis if specified
                    gcparams += ' %s%.5f' % (k, v)
        if len(gcparams) > 0:
            self.gcode += '%s %s\n' % (gc, gcparams)
        elif gc not in ('G01','G02','G03'):
            self.gcode += '%s\n' % gc
    
