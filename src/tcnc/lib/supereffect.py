"""An extension of the Inkscape extension class.
Also includes some utility methods that are handy for generating
SVG output.

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
import gettext
_ = gettext.gettext

import inkex
import simpletransform
import simplepath
import simplestyle

logger = logging.getLogger(__name__)

SUPPORTED_SHAPES = ('path', 'rect', 'line', 'circle',
                   'ellipse', 'polyline', 'polygon')

def optargs(*args, **kwargs):
    """Bundle optarg params.
    
    Handles automatic unit conversions such as degrees to radians, world
    to document, etc.
    """
    options = [args, kwargs, {}]
    # Hoist the conversion specs from the optargs
    if 'convert_from' in kwargs:
        options[2]['convert_from'] = kwargs['convert_from']
        del kwargs['convert_from']
    if 'convert_to' in kwargs:
        options[2]['convert_to'] = kwargs['convert_to']
        del kwargs['convert_to']
    return options
    
class SuperEffect(inkex.Effect):
    """Sort of a beefed up version of inkex.Effect that includes
    more handy methods for generating SVG and traversing the inkscape document.
    """
    doc_size = None
    doc_units = None
    option_conversions = {}
        
    def __init__(self, option_info=None, *args, **kwargs):
        """
        :param option_info: a list of optarg arguments specifying plugin options.
        """
        #super(SuperEffect, self).__init__(*args, **kwargs)
        # Unfortunately inkex.EFfect is declared using old syntax so use old-style init
        inkex.Effect.__init__(self, *args, **kwargs)
        
        # Derived class options will override defaults.
        if option_info:
            for optargs in option_info:
                opt = self.OptionParser.add_option(*optargs[0], **optargs[1])
                # Perform unit conversion if specified as the last two items
                # in the option specification tuple.
                if len(optargs[2]) > 0:
                    from_unit = optargs[2].get('convert_from')
                    to_unit = optargs[2].get('convert_to', 'world')
                    if from_unit == 'doc':
                        from_unit = self.get_document_units()
                    self.option_conversions[opt.dest] = (from_unit, to_unit)

    def convert_option_units(self, default_unit_scale=1.0):
        """Perform any unit conversion on options.
        
        :param default_unit_scale: the default scale factor if
            `convert_from` is not specified in the optargs.
        """
        for option_name in self.options.__dict__:
            if option_name in self.option_conversions:
                from_value = getattr(self.options, option_name)
                to_value = from_value # default
                from_unit, to_unit = self.option_conversions[option_name]
                if from_unit == 'deg' or to_unit == 'rad':
                    to_value = math.radians(from_value)
                elif from_unit == 'rad' or to_unit == 'deg':
                    to_value = math.degrees(from_value)
                else:
                    if from_unit == 'doc':
                        from_unit = self.get_document_units()
                    if from_unit != 'world':
                        if not from_unit:
                            to_value = from_value * default_unit_scale
                        else:
                            to_value = from_value * inkex.uuconv[from_unit]
                    if to_unit != 'world':
                        to_value /= inkex.uuconv[to_unit]
                setattr(self.options, option_name, to_value)
                    
    def get_document_units(self):
        """Return the Inkscape document units (in, mm, etc.).
        """
        if self.doc_units is None:
#            basedoc = self.document.xpath('//sodipodi:namedview[@id="base"]', namespaces=inkex.NSS)
            basedoc = self.getNamedView()
            if basedoc is not None:
                self.doc_units = basedoc.get(inkex.addNS('document-units', 'inkscape'))
                if self.doc_units is None:
                    self.doc_units = basedoc.get('units')
        return self.doc_units
    
    def get_unit_scale(self, units=None):
        """Return the scale factor to scale document units to SVG/world units.
        
        :units: Unit to convert from. Use document units if None (default).
        """
        if units is None:
            units = self.get_document_units()
        return inkex.uuconv[units]
    
    def get_document_size(self):
        """Return width and height of document in SVG units as a tuple (W, H).
        """
        if not self.doc_size:
            x = inkex.unittouu(self.document.getroot().get('width'))
            y = inkex.unittouu(self.document.getroot().get('height'))
            self.doc_size = (x, y)
        return self.doc_size
    
    def get_document_width(self):
        """Return width of document in SVG units."""
        return self.get_document_size()[0]
                
    def get_document_height(self):
        """Return height of document in SVG units."""
        return self.get_document_size()[1]
                
    def find_layer(self, layer_name):
        """Find and return the layer whose layer name (label) is <layer_name>.
        
        Returns None if none found.
        """
#        layer = self.layer_cache.get(layer_name)
#        if layer is None:
        layers = self.document.xpath('//svg:g[@inkscape:label="%s"]' %
                                 layer_name, namespaces=inkex.NSS)
        # If there is more than one layer with the same name just return
        # the first one...
        layer = layers[0] if len(layers) > 0 else None
        return layer
    
    def clear_layer(self, layer_name):
        """Delete the contents of the specified layer.
        Does nothing if the layer doesn't exist.
        """
        layer = self.find_layer(layer_name)
        if layer is not None:
            del layer[:]
            
    def create_layer(self, layer_basename, clear=False, incr_suffix=False):
        """Create an SVG layer or return an existing layer.
        
        :param parentnode: should be document root node
        :param clear: if a layer of the same name already exists then erase it first
            if True otherwise just return it.
        :param incr_suffix: if a layer of the same name already exists then add an
            auto-incrementing numeric suffix to the name (overrides <clear>)
        """
        layer_name = layer_basename
        layer = self.find_layer(layer_name)
        if incr_suffix and layer is not None and len(layer) > 0:
            suffix_n = 0
            while layer is not None and len(layer) > 0:
                layer_name = '%s_%04d' % (layer_basename, suffix_n)
                suffix_n += 1
                layer = self.find_layer(layer_name)
        if layer is None:
            layer = inkex.etree.SubElement(self.document.getroot(), 'g')
            layer.set(inkex.addNS('groupmode', 'inkscape'), 'layer')
            layer.set(inkex.addNS('label', 'inkscape'), layer_name)
#            self.layer_cache[layer_name] = layer
        elif clear:
            del layer[:]
        return layer

    def set_current_layer(self, layer_name):
        """Set the current layer to the specified layer.
        
        If the specified layer doesn't exist then make the default layer
        the current layer.
        """
        layer = self.find_layer(layer_name)
        if layer is None:
            # find the default layer...
            # TODO:
            pass
        if layer is not None:
            self.current_layer = layer
        
    def get_current_layer(self, default_layer):
        """Get the current layer.
        
        If there isn't a current layer then return the specified default layer.
        """
        return self.current_layer if self.current_layer is not None else default_layer
    
    def get_parent_layer(self, node):
        """Return the layer that the node resides in.
        Returns None if the node is not in a layer.
        """
        # This assumes that layers reside at the document root level
        # and that Inkscape still doesn't support sub-layers
        for layer in self.document.getroot():
            if node_is_layer(layer) and node in layer.getchildren():
                return layer
        return None            

    def create_path(self, attrs, style=None, layer=None):
        """Create an SVG path element in the current layer."""
        if layer is None:
            layer = self.current_layer
        if layer is not None:
            if style:
                attrs['style'] = style
            return inkex.etree.SubElement(layer, inkex.addNS('path', 'svg'), attrs)

    def create_circle(self, cx, cy, radius, style=None, layer=None):
        """Create an SVG circle in the current layer."""
        if layer is None:
            layer = self.current_layer
        if layer is not None:
            attrs = {'r': str(radius), 'cx': str(cx), 'cy': str(cy)}
            if style:
                attrs['style'] = style
            return inkex.etree.SubElement(layer, inkex.addNS('circle', 'svg'), attrs)
    
    def create_line(self, x1, y1, x2, y2, style=None, layer=None):
        """Shortcut to create an SVG path consisting of one line segment."""
        attrs = {'d': 'M %5f %5f L %5f %5f' % (x1, y1, x2, y2)}
        return self.create_path(attrs, style, layer)

    def create_arc(self, x1, y1, x2, y2, radius, angle, style=None, layer=None):
        """Create an SVG arc in the current layer."""
        sweep_flag = 0 if angle < 0 else 1
        attrs = { 'd': 'M %5f %5f A %5f %5f 0 0 %d %5f %5f' % \
                  (x1, y1, radius, radius,
                   sweep_flag, x2, y2),}
        return self.create_path(attrs, style, layer)

    def create_polygon(self, vertices, style=None, layer=None):
        d = 'M %f %f L' % (vertices[0].x, vertices[0].y)
        for p in vertices[1:]:
            d = d + ' %f,%f' % (p.x, p.y)
        if vertices[0] != vertices[-1]:
            d = d + ' %f,%f' % (vertices[0].x, vertices[0].y)
        attrs = {'d': d}
        return self.create_path(attrs, style, layer)

    def create_simple_marker(self, marker_id, d, style, transform):
        """Create an Inkscape line end marker glyph.
        :param parentnode: should be document root node
        """
        docroot = self.document.getroot()
        defs = docroot.find(inkex.addNS('defs', 'svg'))
        if defs is None:
            defs = inkex.etree.SubElement(docroot, inkex.addNS('defs', 'svg'))
        marker = inkex.etree.SubElement(defs, inkex.addNS('marker','svg'),
                        {'id': marker_id, 'orient': 'auto', 'refX':  '0.0',
                         'refY': '0.0', 'style': 'overflow:visible',})
        inkex.etree.SubElement( marker, inkex.addNS('path','svg'), 
                        { 'd': d, 'style': style, 'transform': transform, })
        return marker
    

def node_is_layer(node):
    """Return True if the node is an Inkscape layer root node."""
    return node.get(inkex.addNS('groupmode', 'inkscape')) == 'layer'

def flatten_nodetree(nodetree, mtransform=[[1.0, 0.0, 0.0],[0.0, 1.0, 0.0]],
                     parent_visibility='visible', ignored_layers=None, _nodelist=None):
    """Recursively traverse an SVG node tree and flatten it to a list of 
    tuples containing an SVG shape element and its accumulated transform.
    
    This does a depth-first traversal of <g> and <use> elements.
    Anything besides paths, rectangles, circles, ellipses, lines, polygons,
    and polylines are ignored.
    
    Invisible elements are ignored.
    
    :param mtransform: a transform matrix to add to each node's transforms
    :param ignored_layers: a list of names of layers that should be skipped
    
    Returns a list of SVG elements.
    """
    if _nodelist is None:
        _nodelist = []
    for node in nodetree:
        # Ignore invisible nodes
        v = node.get('visibility', parent_visibility)
        if v == 'inherit':
            v = parent_visibility
        if v == 'hidden' or v == 'collapse':
            continue
#        logger.debug('flattening %s' % node.get('id'))
        # first apply the current transform matrix to this node's tranform
        transform = simpletransform.parseTransform(node.get('transform'))
#        logger.debug('transform=%s' % node.get('transform'))
        mtransform_node = simpletransform.composeTransform(mtransform, transform)

        if node.tag == inkex.addNS('g', 'svg') or node.tag == 'g':
            if ignored_layers and node_is_layer(node):
                layer_name = node.get(inkex.addNS('label', 'inkscape'))
                layer_display = simplestyle.parseStyle(node.get('style')).get('display')
                if layer_name in ignored_layers or layer_display == 'none':
                    continue
            # Recursively traverse groups
#            logger.debug('flattening layer %s' % layer_name)
            flatten_nodetree(node, mtransform_node, parent_visibility=v,
                             _nodelist=_nodelist)

        elif node.tag == inkex.addNS('use', 'svg') or node.tag == 'use':
            # A <use> element refers to another SVG element via an xlink:href="#blah"
            # attribute.  We will handle the element by doing an XPath search through
            # the document, looking for the element with the matching id="blah"
            # attribute.  We then recursively process that element after applying
            # any necessary (x,y) translation.
            #
            # Notes:
            #  1. We ignore the height and width attributes as they do not apply to
            #     path-like elements, and
            #  2. Even if the use element has visibility="hidden", SVG still calls
            #     for processing the referenced element.  The referenced element is
            #     hidden only if its visibility is "inherit" or "hidden".
            refid = node.get( inkex.addNS( 'href', 'xlink' ) )
            if refid:
                # [1:] to ignore leading '#' in reference
                path = '//*[@id="%s"]' % refid[1:]
                refnode = node.xpath( path )
                if refnode:
                    x = float( node.get( 'x', '0' ) )
                    y = float( node.get( 'y', '0' ) )
                    # Note: the transform has already been applied
                    if ( x != 0 ) or (y != 0 ):
                        mtransform_node = simpletransform.composeTransform(mtransform_node, simpletransform.parseTransform('translate(%f,%f)' % (x,y)))
                    v = node.get( 'visibility', v )
                    flatten_nodetree(refnode, mtransform_node,
                                     parent_visibility=v, _nodelist=_nodelist)

        elif get_node_tag(node) in SUPPORTED_SHAPES:
            # Set this node's transform to the combined transform
            #node.set('transform', simpletransform.formatTransform(mtransform_node))
#            logger.debug('adding node: %s' % node.get('id'))
            _nodelist.append((node, mtransform_node))

        else:
            # silently ignore text, images, etc.
            pass

    return _nodelist


def get_node_tag(node):
    """Get the node tag stripped of it's namespace part if any"""
    return node.tag.rpartition('}')[2]


def convert_rect_to_path(node):
    """Convert an SVG rect shape element to a simplepath.
    
    Convert this:
       <rect x='X' y='Y' width='W' height='H'/>
    to this:
       'M X1 Y1 L X1 Y2 L X2 Y2 L X2 Y1 Z'
    """
    x1 = float(node.get('x', 0))
    y1 = float(node.get('y', 0))
    x2 = x1 + float(node.get('width', 0))
    y2 = y1 + float(node.get('height', 0))
    return [['M', [x1, y1]], ['L', [x1, y2]], ['L', [x2, y2]],
            ['L', [x2, y1]], ['L', [x1, y1]]]
    #return 'M %f %f L %f %f L %f %f L %f %f Z' % (x1, y1, x1, y2, x2, y2, x2, y1)


def convert_line_to_path(node):
    """Convert an SVG line shape element to a simplepath.
    
    Convert this:
       <line x1='X1' y1='Y1' x2='X2' y2='Y2/>
    to this:
       'MX1 Y1 LX2 Y2'
    """
    x1 = float(node.get('x1', 0))
    y1 = float(node.get('y1', 0))
    x2 = float(node.get('x2', 0))
    y2 = float(node.get('y2', 0))
    return [['M', [x1, y1]], ['L', [x2, y2]]]
    #return 'M %s %s L %s %s' % (node.get('x1'), node.get('y1'), node.get('x2'), node.get('y2'))


def convert_circle_to_path(node, subdivision=1):
    """Convert an SVG circle shape element to a simplepath.
    The circle is divided into <subdivision> number of circular arcs.
    
    Convert this:
       <circle r='RX' cx='X' cy='Y'/>
    to this:
       'MX1,CY A RX,0 0 1 0 X2,CY [A ...]'
    """
    r = float(node.get('r', 0))
    cx = float(node.get('cx', 0))
    cy = float(node.get('cy', 0))
    x = cx - r
    if subdivision == 1:
        d = [['M', [x, cy]], ['A', [r, r, 0, 1, 0, x, cy]]]
        #d = 'M %f %f A %f,%f 0 1 0 %f,%f' % (x, cy, r, r, x, cy)
    else:
        d = [['M', [x, cy]],]
        #d = 'M %f %f' % (x, cy)
        y = cy
        arc = math.pi / subdivision
        for i in range(1, subdivision):
            x = math.cos(arc * i) * r
            y = math.sin(arc * i) * r
            d.append(['A', [r, r, 0, 1, 0, x, y]])
            #d += ' A %f,%f 0 1 0 %f,%f' % (r, r, x, y)
    return d


def convert_ellipse_to_path(node):
    """Convert an SVG ellipse shape element to a path.
    The ellipse is divided into two 180deg elliptical arcs.
    
    Convert this:
       <ellipse rx='RX' ry='RY' cx='X' cy='Y'/>
    to this:
       'M X1,CY A RX,RY 0 1 0 X2,CY A RX,RY 0 1 0 X1,CY'
    """
    rx = float(node.get('rx', 0))
    ry = float(node.get('ry', 0))
    cx = float(node.get('cx', 0))
    cy = float(node.get('cy', 0))
    x1 = cx - rx
    x2 = cx + rx
    a1 = [rx, ry, 0, 1, 0, x2, cy]
    a2 = [rx, ry, 0, 1, 0, x1, cy]
    return [ ['M', [x1, cy]], ['A', a1], ['A', a2] ]
#    return 'M %f %f ' % ( x1, cy ) + \
#        'A %f,%f 0 1 0 %f,%f ' % ( rx, ry, x2, cy ) + \
#        'A %f,%f 0 1 0 %f,%f' % ( rx, ry, x1, cy )


def convert_polyline_to_path(node):
    """Convert an SVG line shape element to a path.
    
    Convert this:
       <polyline points='x1,y1 x2,y2 x3,y3 [...]'/>
    to this:
       'M x1 y1 L x2 y2 L x3 y3 [...]'/>
    """
    points = node.get('points', '').split()
    point = points[0].split(',')
    d = [ [ 'M', [float(point[0]), float(point[1])] ], ]
    #d = 'M ' + ''.join(points[0].split(',')) # remove comma separator
    for i in range(1, len(points)):
        point = points[i].split(',')
        d.append( ['L', [float(point[0]), float(point[1])]] )
        #d += ' L ' + ''.join(points[i].split(','))
    return d


def convert_polygon_to_path(node):
    """Convert an SVG line shape element to a path.
    
    Convert this:
       <polygon points='x1,y1 x2,y2 x3,y3 [...]'/>
    to this:
       'M x1 y1 L x2 y2 L x3 y3 [...]'/>
    """
    d = convert_polyline_to_path(node)
    #d += ' Z' # close path for polygons
    d.append(['L', d[0][1]])
    return d


def convert_element_to_path(node):
    """Convert an SVG element into a simplepath.
    
    This handles paths, rectangles, circles, ellipses, lines, and polylines.
    Anything else raises and exception.
    """    
    node_tag = get_node_tag(node) # node tag stripped of namespace part
        
    if node_tag == 'path':
        return simplepath.parsePath(node.get('d'))
    elif node_tag == 'rect':
        return convert_rect_to_path(node)
    elif node_tag == 'line':
        return convert_line_to_path(node)
    elif node_tag == 'circle':
        return convert_circle_to_path(node)
    elif node_tag == 'ellipse':
        return convert_ellipse_to_path(node)    
    elif node_tag == 'polyline':
        return convert_polyline_to_path(node)    
    elif node_tag == 'polygon':
        return convert_polygon_to_path(node)
    elif node_tag == 'text':
        raise Exception(_('Unable to convert text; please convert text to paths first.'))
    elif node_tag == 'image':
        raise Exception(_('Unable to convert bitmap images; please convert them to line art first.'))
    else:
        raise Exception(_('Unable to convert this SVG element to a path: <%s>') % (node.tag))


def create_path_element(node, d=None, style=None, transform=None):
    """Create an SVG path element.
    
    If no path is specified then copy it from <node> if <node> is
    a path or convert <node> shape (ie rect, line, etc) into a path.
    Copies style and transform attributes from <node> if none given.
    """
    if d is None:
        d = node.get('d', simplepath.formatPath(convert_element_to_path(node)) )
    newpath = inkex.etree.Element(inkex.addNS('path', 'svg'))
    newpath.set('d', d);
    if not style:
        style = node.get('style')
    if style:
        newpath.set('style', style)
    if not transform:
        transform = node.get('transform')
    if transform:
        newpath.set('transform', transform)
    return newpath
