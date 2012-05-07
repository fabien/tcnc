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

import inkex
import cubicsuperpath
import simpletransform

import svg
import geom
import gcode

from geom import P, CubicBezier

tcnc_current_version = "0.1"

DEFAULT_LOG_LEVEL = logging.DEBUG
logger = logging.getLogger(__name__)

PREVIEW_LAYER_NAME = 'tcnc-preview'
DEBUG_LAYER_NAME = 'tcnc-debug'

option_info = (
    ('--active-tab',        'store', 'string',         'active_tab', '',            'Defines which tab is active'),
    
    ('--unit',                'store', 'string',         'unit', 'in',                'Units'),
    ('--auto-select-paths',    'store', 'inkbool',        'auto_select_paths', True,    'Select all paths if nothing is selected.'),
    ('--biarc-tolerance',    'store', 'float',         'biarc_tolerance', '1',        'Biarc approximation fitting tolerance.'),                
    ('--biarc-maxdepth',    'store', 'float',         'biarc_max_depth', '4',        'Biarc approximation maximum curve splitting recursion depth.'),                
    ('--line-flatness',    'store', 'float',         'line_flatness', '0.5',        'Curve to line flatness.'),                
    ('--min-arc-radius',        'store', 'float',         'min_arc_radius', '.01',        'All arcs having radius less than minimum will be considered as straight line'),

    ('--origin-ref',            'store', 'string',         'origin_ref', 'paper',        'Lower left origin reference.'),
    ('--z-scale',            'store', 'float',         'z_scale', '1.0',            'Scale factor Z'),    
    ('--z-offset',            'store', 'float',         'z_offset', '0.0',            'Offset along Z'),
    ('--x-scale',            'store', 'float',         'x_scale', '1.0',            'Scale factor X'),                
    ('--x-offset',            'store', 'float',         'x_offset', '0.0',            'Offset along X'),
    ('--y-scale',            'store', 'float',         'y_scale', '1.0',            'Scale factor Y'),                
    ('--y-offset',            'store', 'float',         'y_offset', '0.0',            'Offset along Y'),
    ('--a-offset',            'store', 'float',         'a_offset', '0.0',            'Angular offset along rotational axis'),

    ('--xy-feed',            'store', 'float',         'xy_feed', '10.0',            'XY axis feed rate in unit/s'),
    ('--z-feed',                'store', 'float',         'z_feed', '10.0',            'Z axis feed rate in unit/s'),
    ('--a-feed',                'store', 'float',         'a_feed', '60.0',            'A axis feed rate in deg/s'),
    ('--z-safe',                'store', 'float',         'z_safe', '5.0',            'Z axis safe height for rapid moves'),

    ('--directory',         'store', 'string',         'directory', '~',            'Directory for gcode file'),
    ('--filename',          'store', 'string',         'filename', '-1.0',            'File name'),            
    ('--append-suffix',        'store', 'inkbool',        'append_suffix', True,        'Append auto-incremented numeric suffix to filename'),            
    ('--create-log',            'store', 'inkbool',     'log_create_log', True,    'Create log files'),
    ('--log-level',            'store', 'string',         'log_level', 'DEBUG',        'Log level'),
    ('--log-filename',        'store', 'string',         'log_filename', 'tcnc.log',            'Full pathname of log file'),

    ('--preview-show',            'store', 'inkbool',     'preview_show', True,    'Show generated cut paths on preview layer.'),
    ('--debug-layer',            'store', 'inkbool',     'debug_layer', True,    'Create debug layer.'),
    ('--debug-biarcs',            'store', 'inkbool',     'debug_biarcs', True,    ''),

    ('--z-depth',            'store', 'float',         'z_depth', '-0.125',        'Z full depth of cut'),
    ('--z-step',                'store', 'float',         'z_step', '-0.125',            'Z cutting step depth'),
    ('--path-to-gcode-order','store', 'string',         'path_to_gcode_order', 'path by path',    'Defines cutting order path by path or layer by layer.'),        
    ('--path-to-gcode-depth-function','store', 'string',     'path_to_gcode_depth_function', 'zd',    'Path to gcode depth function.'),
    ('--biarc-max-split-depth',        'store', 'int',         'biarc_max_split_depth', '4',            'Defines maximum depth of splitting while approximating using biarcs.'),
)

class TCnc(svg.SuperEffect):
    '''Inkscape plugin that converts selected SVG elements into gcode suitable for a
    four axis (XYZA) CNC machine with a tangential tool (ie a knife or a brush) as the A axis.
    '''
    styles = {
              'simple': 'fill:none;stroke:#cccc99;stroke-width:0.25pt;stroke-opacity:1',
              'cutline': 'fill:none;stroke:#c000c0;stroke-width:0.5pt;stroke-opacity:1;marker-end:url(#PreviewLineEnd0)',
              'cutline_tiny': 'fill:none;stroke:#0000ff;stroke-width:1.0pt;stroke-opacity:1',
              'cutarc0': 'fill:none;stroke:#ff0000;stroke-width:0.5pt;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:4;stroke-opacity:1;marker-end:url(#PreviewLineEnd0)',
              'cutarc1': 'fill:none;stroke:#ff0000;stroke-width:0.5pt;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:4;stroke-opacity:1;marker-end:url(#PreviewLineEnd1)',
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
        
        self.docroot = self.document.getroot()
        self.create_inkscape_markers()
            
        # Get selected SVG elements if any
        rootnodes = self.selected.values()
        if not rootnodes:
            # Use entire document if nothing is selected
            rootnodes = self.docroot
        
        # Clear out the old preview layer if any
        self.clear_layer(PREVIEW_LAYER_NAME)
        self.clear_layer(DEBUG_LAYER_NAME)

        # Extract all the recognized SVG shape elements
        skip_layers=(PREVIEW_LAYER_NAME,)
        nodelist = svg.flatten_nodetree(rootnodes, skip_layers=skip_layers)
        
        # Create a new layer that will contain the G code preview
        self.preview_layer = self.create_layer(PREVIEW_LAYER_NAME)
        #self.set_current_layer(PREVIEW_LAYER_NAME)
        
        if self.options.debug_layer:
            debug_layer = self.create_layer(DEBUG_LAYER_NAME)
            # setup geom module for debug output
            geom.DEBUG_EFFECT = self
            geom.DEBUG_LAYER = debug_layer
        
        self.current_layer = self.preview_layer        

        # Process the SVG elements and generate G code
        self.process_shapes(nodelist)
            
            
    def create_inkscape_markers(self):
        '''Create Inkscape line end marker glyphs and insert them into the document.
        '''
        self.create_simple_marker('PreviewLineEnd0', 'M 5.77,0.0 L -2.88,5.0 L -2.88,-5.0 L 5.77,0.0 z',
                          self.styles['cutpath_end_marker'], 'scale(0.4) translate(-4.5,0)')
        self.create_simple_marker('PreviewLineEnd1', 'M 5.77,0.0 L -2.88,5.0 L -2.88,-5.0 L 5.77,0.0 z',
                          self.styles['cutpath_end_marker'], 'scale(-0.4) translate(-4.5,0)')
        self.create_simple_marker('PreviewLineEnd2', 'M 5.77,0.0 L -2.88,5.0 L -2.88,-5.0 L 5.77,0.0 z',
                          self.styles['movepath_end_marker'], 'scale(0.5) translate(-4.5,0)')


    def process_shapes(self, nodelist):
        '''Process the SVG shape elements and generate G code.
        Return a list of cut paths as SVG path elements.
        '''
        biarc_tolerance = float(self.options.biarc_tolerance)
        biarc_max_depth = float(self.options.biarc_max_depth)
        line_flatness = float(self.options.line_flatness)
        logging.debug('Biarc tolerance & depth: %f %d' % (biarc_tolerance, biarc_max_depth))
        cutpath_list = []
        for node in nodelist:
            # Convert the shape element to a simplepath
            path = svg.convert_element_to_path(node)
            
            # Convert the simplepath to a 'cubicsuperpath' which is
            # just a list of cubic bezier curves.
            # This seems to be how most Inkscape plugins do things...
            # TODO: This is really an unnecessary step since
            # we could deal with SVG shapes directly...
            csp = cubicsuperpath.CubicSuperPath(path)
            
            # Apply the SVG element transform to the path segments so that
            # we are working in absolute coordinates
            transform = node.get('transform')
            if transform:
                logging.debug('transform: ' + transform)
                mat = simpletransform.parseTransform(transform)
                simpletransform.applyTransformToPath(mat, csp)
            
            # Convert cubic path segments to arcs and lines
            for subcsp in csp:
                logging.debug('converting path segments to biarcs')
                cutpath = []
                for i in range(1,len(subcsp)):
                    sp1 = subcsp[i-1]
                    sp2 = subcsp[i]
                    curve = CubicBezier(P(sp1[1]), P(sp1[2]), P(sp2[0]), P(sp2[1]))
#                    curve.SVG_plot()
                    biarcs = curve.biarc_approximation(tolerance=biarc_tolerance,
                                                       max_depth=biarc_max_depth,
                                                       line_flatness=line_flatness)
                    cutpath.extend(biarcs)
                    
                cutpath_list.append(cutpath)
        
        # Generate and export G code
        gcodebuf = ''
        for cutpath in cutpath_list:
            gcodebuf += self.generate_gcode(cutpath)
        self.export_gcode(gcodebuf)
#        logging.info('gcode:\n' + gcodebuf + '\nend gcode\n')
                             
        return
    
    def draw_preview_line(self, x1, y1, x2, y2, style_id=None):
        '''Draw an SVG line path on to the preview layer'''
        if style_id is None:
            style_id = 'simple'     
        self.create_line(x1, y1, x2, y2, self.styles[style_id])
    
    def draw_preview_Line(self, line, style_id=None):
        '''Draw an SVG line path on to the preview layer'''
        self.draw_preview_line(line.p1.x, line.p1.y, line.p2.x, line.p2.y, style_id)
        
    def draw_preview_arc(self, x1, y1, r, sweep_flag, x2, y2, style_id):
        '''Draw an SVG arc on to the preview layer'''
        attrs = { 'd': 'M %5f %5f A %5f %5f 0.0 0 %d %5f %5f' % (x1, y1, r, r, sweep_flag, x2, y2),
                 'style': self.styles[style_id + str(sweep_flag)] }
        self.create_path(attrs)
        
    def draw_preview_path(self, d, style_id):
        '''Draw an SVG path on the preview layer'''
        attrs = { 'd': d, 'style': self.styles[style_id]}
        self.create_path(attrs)

    def draw_dot(self, x, y, size='small', color='#000000'):
        '''Draw an SVG dot (a small filled circle) at the specified location'''
        radius = {'small': '2pt', 'medium': '5pt', 'large': '10pt',}
        style = 'fill:%s;stroke:%s' % (color, color)
        self.create_circle(x, y, radius[size], style)

    def generate_gcode(self, cutpath, depth=0):
        """Generate G code for the given cutpath.
        This method also creates SVG output that serves as a preview to the
        G code output. The SVG output is added to the specified group or layer.
        The cutpath is specified by a list of curve segments which are either
        line segments or circular arc segments.
        """
        def calc_rotation(current_angle, new_angle):
            '''Calculate the relative rotation amount in radians'''
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
            
        if len(cutpath) == 0:
            return ''                
        
        gc = gcode.GCode(
                    zsafe=float(self.options.z_safe),
                    zfeed=float(self.options.z_feed),
                    xyfeed=float(self.options.xy_feed),
                    afeed=float(self.options.a_feed),
                    atolerance=geom.EPSILON,
                    )
        
        current_angle = 0.0
        
        # Move to the start point of the first path segment.
        # Move the tool rapidly to the end point
        gc.tool_up()
        gc.rapid_move(cutpath[0].p1.x, cutpath[0].p1.y, a=current_angle)
        
        # Draw the initial rapid move line on the SVG preview layer
        self.draw_preview_line(self.last_point[0], self.last_point[1],
                               cutpath[0].p1.x, cutpath[0].p1.y, 'moveline')            

        # Create G-code for each segment of the cutpath
        for segment in cutpath:
            start_x = segment.p1.x
            start_y = segment.p1.y
            end_x = segment.p2.x
            end_y = segment.p2.y
            start_z = 0 + depth
            end_z = 0 + depth

            self.last_point = (end_x, end_y)
            
            if isinstance(segment, geom.Line):
                # Calculate tool tangent angle for line
                a = math.atan2(end_y - start_y, end_x - start_x)
                current_angle += calc_rotation(current_angle, a)
                gc.feed_rotate(current_angle)
                
                # Don't bother outputting code if the line is too small
                if math.hypot(end_x - start_x, end_y - start_y) > geom.EPSILON:
                    gc.tool_down(start_z)    
                    gc.feed(end_x, end_y, end_z)
                    
                    # Add the line to the SVG preview layer
                    self.draw_preview_line(start_x, start_y, end_x, end_y, 'cutline')
                
            elif isinstance(segment, geom.Arc):
                # Calculate starting tool tangent angle
                center_x = segment.center.x
                center_y = segment.center.y
                arc_angle = segment.angle
                if arc_angle < 0: # CW ?
                    a = math.atan2(center_y - start_y, center_x - start_x) + math.pi/2
                    current_angle += calc_rotation(current_angle, a)
                else: # CCW
                    a = math.atan2(start_y - center_y, start_x - center_x) + math.pi/2
                    current_angle += calc_rotation(current_angle, a)
                gc.feed_rotate(current_angle)    # Rotate to starting arc tangent
                current_angle += arc_angle       # Endpoint tangent
                    
                gc.tool_down(start_z)    
    
                arc_radius = math.hypot(start_x - center_x, start_y - center_y)
                if arc_radius > self.options.min_arc_radius:
                    gc.feed_arc((arc_angle<0), end_x, end_y, end_z, (center_x-start_x), (center_y-start_y), current_angle)
                    
                    # Add the arc to the SVG preview layer
                    sweep_flag = 0 if arc_angle < 0 else 1
                    self.draw_preview_arc(start_x, start_y, arc_radius, sweep_flag, end_x, end_y, 'cutarc')
                    
                elif math.hypot(end_x - start_x, end_y - start_y) > geom.EPSILON:
                    # Arc radius was too small to draw so calculate tangent for line
                    a = math.atan2(end_y - start_y, end_x - start_x)
                    current_angle += calc_rotation(current_angle, a)
                    gc.feed_rotate(current_angle)
                    # Draw a very short line
                    gc.feed(end_x, end_y, end_z)
                    
                    # Add the line to the SVG preview layer
                    self.draw_preview_line(start_x, start_y, end_x, end_y, 'cutline')

        # Post-path G code
        gc.tool_up()
        gc.rehome_rotational_axis()
        return gc.gcode
    
    def export_gcode(self, gcodebuf):
        '''Export the generated g code to a specified file.'''
        filedir = os.path.expanduser(self.options.directory)
        filename = os.path.basename(self.options.filename)
        path = os.path.abspath(os.path.join(filedir, filename))
        if self.options.append_suffix:
            file_root, file_ext = os.path.splitext(filename)
            # Get a list of existing files that match the numeric suffix.
            # They should already be sorted.
            filter_str = '%s_[0-9]*%s' % (file_root, file_ext)
            files = fnmatch.filter(os.listdir(filedir), filter_str)
            logging.debug('filter: %s [%s]' % (filter_str, str(files)))
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
#        logging.debug('path: %s' % path)
        
        try:
            with open(path, 'w') as f:
                f.write(gcodebuf)
        except Exception:
            inkex.errormsg("Can't write to file: %s" % path)
        

tcnc = TCnc(option_info)
tcnc.affect()