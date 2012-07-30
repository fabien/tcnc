'''
Python port of quasi.c which was originally written by Eric Weeks
link: http://www.physics.emory.edu/~weeks/software/quasi.html
email: weeks@physics.emory.edu

Mostly unchanged except to make it a little more pythonic...
'''

import math
from optparse import OptionParser
from Carbon.Aliases import false

def quasi(zfill, midon, symmetry, maxmax, plotter):
	vx = []
	vy = []
	mm = []
	b = []
	halfmax = maxmax / 2

	# Initialize vectors
	for t in range(symmetry):
		phi = ((t * 2.0) / symmetry) * math.pi
		x = math.cos(phi)
		y = math.sin(phi)
		m = y / x
		vx.append(x)
		vy.append(y)
		mm.append(m)
		vr = []
		for r in range(maxmax):
			y1 = y * (t * 0.1132) - x * (r - halfmax)  # offset 
			x1 = x * (t * 0.2137) + y * (r - halfmax)
			vr.append(y1 - m * x1)		# intercept
		b.append(vr)

	# t is 1st direction, r is 2nd.  look for intersection between pairs
	# of lines in these two directions. (will be x0,y0) 

	index = [0,] * symmetry
	fill_color = 0.2
	stroke_color = 0.0
	themax = maxmax - 1
	themin = themax / 2
	minmin = 0.0
	midsix = 0
	
	while minmin <= float(themax):
		rad1 = minmin * minmin
		minmin += 0.4
		rad2 = minmin * minmin
		for n in range(1, themax):
			for m in range(1, themax):
				rad = float((n-themin) * (n-themin) + (m-themin) * (m-themin))
				if rad >= rad1 and rad < rad2:
					for t in range(symmetry-1):
						for r in range(t + 1, symmetry):
							x0 = (b[t][n] - b[r][m]) / (mm[r] - mm[t])
							y0 = mm[t] * x0 + b[t][n]
							flag = False
							for i in range(symmetry):
								if i != t and i != r:
									dx = -x0 * vy[i] + (y0 - b[i][0]) * vx[i]
									index[i] = int(-dx)
									if index[i] > maxmax-3 or index[i] < 1:
										flag = True
							if not flag:
								index[t] = n - 1
								index[r] = m - 1
								x0 = 0.0
								y0 = 0.0
								for i in range(symmetry):
									x0 += vx[i] * index[i]
									y0 += vy[i] * index[i]
								vertices = [(x0, y0)]
								x0 += vx[t]
								y0 += vy[t]
								vertices.append((x0, y0))
								x0 += vx[r]
								y0 += vy[r]
								vertices.append((x0, y0))
								x0 -= vx[t]
								y0 -= vy[t]
								vertices.append((x0, y0))
								x0 -= vx[r]
								y0 -= vy[r]
								vertices.append((x0, y0))
								if midon[0] > 0:
									stroke_color = 0.8 # faint lines
								# color of tile unless zfill==1
								fill_color += .05
								if fill_color > 1.0:
									fill_color = 0.2
								if zfill:
									fill_color = 0.0
									for i in range(symmetry):
										fill_color += index[i]
									while fill_color > (symmetry - 1.0) / 2.0:
										fill_color -= (symmetry - 1.0) / 2.0
									fill_color = fill_color / ((symmetry - 1.0) / 2.0) * 0.8 + 0.1
									fill_color += math.fabs(vx[t] * vx[r] + vy[t] * vy[r]) # dot product
									if fill_color > 1.0:
										fill_color -= 1.0
								plotter.set_fill_color(fill_color)
								plotter.set_stroke_color(stroke_color)
								plotter.plot_polygon(vertices)
								if midon[0] > 0:
									midx1 = x0 + vx[t]*0.5
									midy1 = y0 + vy[t]*0.5
									midx2 = x0 + vx[t] + vx[r]*0.5
									midy2 = y0 + vy[t] + vy[r]*0.5
									midx3 = x0 + vx[r] + vx[t]*0.5
									midy3 = y0 + vy[r] + vy[t]*0.5
									midx4 = x0 + vx[r]*0.5
									midy4 = y0 + vy[r]*0.5
									dx1 = midx1 - midx2
									dy1 = midy1 - midy2
									dist1 = dx1 * dx1 + dy1 * dy1
									dx2 = midx2 - midx3
									dy2 = midy2 - midy3
									dist2 = dx2 * dx2 + dy2 * dy2
									if dist1 * dist2 < 0.1:
										segtype = midon[0]
									else:
										segtype = midon[1]
									if segtype == 1 or segtype == 2:
										if dist1 > dist2:
											segtype = 3 - segtype
									elif segtype == 5:
										midsix = 1 - midsix
										segtype = midsix + 1
									elif segtype == 6:
										midsix += 1
										if midsix > 2:
											midsix = 0
										segtype = midsix + 1
									stroke_color = 0.0 # dark lines
									plotter.set_stroke_color(stroke_color)
									if segtype == 3:
										# X's 
										plotter.plot_segment(midx1, midy1, midx3, midy3)
										plotter.plot_segment(midx2, midy2, midx4, midy4)
									elif segtype == 1:
										plotter.plot_segment(midx1, midy1, midx2, midy2)
										plotter.plot_segment(midx3, midy3, midx4, midy4)
									elif segtype == 2:
										plotter.plot_segment(midx1, midy1, midx4, midy4)
										plotter.plot_segment(midx2, midy2, midx3, midy3)
									elif segtype == 4:
										# boxes
										plotter.plot_segment(midx1, midy1, midx2, midy2)
										plotter.plot_segment(midx3, midy3, midx4, midy4)
										plotter.plot_segment(midx1, midy1, midx4, midy4)
										plotter.plot_segment(midx2, midy2, midx3, midy3)



class QuasiPlotter(object):
	''''''
	fill_color = 1.0
	stroke_color = 0.0
	use_color = false
	
	def set_fill_color(self, fill_color):
		self.fill_color = fill_color

	def set_stroke_color(self, stroke_color):
		self.stroke_color = stroke_color
		
	def get_fill_color(self):
		return self.fill_color
	
	def get_stroke_color(self):
		return self.stroke_color
	
	def get_fill_color_rgb(self):
		if self.use_color:
			phi = self.fill_color * 2.0 * math.pi
			theta = (math.pi / 2.0) * math.sin(3.0 * phi) + math.pi / 2.0
			r = 0.5 * math.sin(theta) * math.cos(phi) + 0.5
			g = 0.5 * math.sin(theta) * math.sin(phi) + 0.5
			b = 0.5 * math.cos(theta) + 0.5
			return (r, g, b)
		else:
			return (self.fill_color, self.fill_color, self.fill_color)
		
	def plot_polygon(self, vertices):
		''''''
		pass

	def plot_segment(self, x1, y1, x2, y2):
		''''''
		pass

class PSQuasiPlotter(QuasiPlotter):
	'''Postscript plotter'''
	clipped = False
	scale = 1.0
	offsetx = 0.0
	offsety = 0.0
	xcenter = 0.0
	ycenter = 0.0
	rotate = False
	window = 1.0
	stroke_width = 0.015
	
	def __init__(self, *args, **kwargs):
		super(PSQuasiPlotter, self).__init__(*args, **kwargs)		
	
	def plot_polygon(self, vertices):
		''''''
		self._psplot(vertices[0][0], vertices[0][1], 0)
		for p in vertices[1:-1]:
			self._psplot(p[0], p[1], 1)
		self._psplot(vertices[-1][0], vertices[-1][1], 2)

	def plot_segment(self, x1, y1, x2, y2):
		''''''
		self._psplot(x1, y1, 0)
		self._psplot(x2, y2, 2)

	def ps_header(self):
		'''Output header'''
		# PostScript Header (taken from CGLE output)
		print("%!PS-Adobe-1.0 ")
		print("%%BoundingBox: -1 -1 766.354 567.929 ")
		print("%%EndComments ")
		print("%%EndProlog ")
		print("gsave ")
		print(" ")
		print("/f {findfont exch scalefont setfont} bind def ")
		print("/s {show} bind def ")
		print("/ps {true charpath} bind def ")
		print("/l {lineto} bind def ")
		print("/m {newpath moveto} bind def ")
		print("/sg {setgray} bind def")
		print("/a {stroke} bind def")
		print("/cp {closepath} bind def")
		print("/g {gsave} bind def")
		print("/h {grestore} bind def")
		print("matrix currentmatrix /originmat exch def ")
		print("/umatrix {originmat matrix concatmatrix setmatrix} def ")
		print(" ")
		print("% Flipping coord system ")
		print("[8.35928e-09 28.3465 -28.3465 8.35928e-09 609.449 28.6299] umatrix ")
		print("[] 0 setdash ")
		print("0 0 0 setrgbcolor ")
		print("0 0 m ")
		print("%.6f setlinewidth " % self.stroke_width)
	
	def ps_footer(self):
		'''Output footer'''
		print("showpage grestore ")
		print("%%Trailer")
	
	def getdx(self, x, center):
		dx = (x - center) / self.window
		dx = 0.5 * (dx + 1.0)
		return dx

	def _psplot(self, x, y, plotflag):
		swap = 0.0
		cmx = 0.0
		cmy = 0.0		# x,y in centimeters 
		# plotflag:  0 = start line; 1 = lineto; 2 = endpoint 
	
		dx = self.getdx(x, self.xcenter)
		dy = self.getdx(y, self.ycenter)
		if self.rotate:
			swap = dx
			dx = dy
			dy = swap
	
		if ((dx<1.3) and (dy<1.0) and (dx>0) and (dy>0)):		# in window 
			cmx = (dx * self.scale) +  self.offsetx
			cmy = (dy * self.scale)  +  self.offsety
			if plotflag < 1:
				print("%.2f %.2f m" % (cmx, cmy))
			else:
				if self.clipped:
					print("%.2f %.2f m" % (cmx,cmy))
				print("%.2f %.2f l" % (cmx,cmy))
				if plotflag == 2:
					print("cp")
					if self.fillon:
						print("g")
						if self.use_color:
							print("%.2f %.2f %.2f sc" % self.get_fill_color_rgb())
						else:
							print("%.1f sg" % self.get_fill_color())
						print("fill")
						print("h")
					else:
						print("%.1f sg" % self.get_stroke_color())
					print("a")
			self.clipped = False
		else:
			self.clipped = True


def main():
	parser = OptionParser()

	parser.add_option('-s', type='float', default='15.0', dest='scale', help='size of output (width) in cm [15]')
	parser.add_option('-f', type='int', default='0', dest='fillon', help='fill polygons')
	parser.add_option('-z', type='int', default='0', dest='zfill', help='fill color is according to polygon type')
	parser.add_option('-F', type='int', default='0', dest='rotate', help='flip picture 90 degrees')
	parser.add_option('-m', type='float', default='1.0', dest='magnifya', help='magnification factor')
	parser.add_option('-M', type='float', default='0.0', dest='midon0', help='midpoint type for skinny diamonds')
	parser.add_option('-N', type='float', default='0.0', dest='midon1', help='midpoint type for fat diamonds')
	parser.add_option('-S', type='int', default='5', dest='symmetry', help='degrees of symmetry [5]')
	parser.add_option('-n', type='int', default='30', dest='maxmax', help='number of lines to use [30]')
	parser.add_option('-w', type='float', default='20.0', dest='window', help='clipping window')

	(options, args) = parser.parse_args()

	fillon = bool(options.fillon)
	zfill = bool(options.zfill)
	midon = [int(options.midon0), int(options.midon1)]
	symmetry = int(options.symmetry)
	maxmax = int(options.maxmax)

	plotter = PSQuasiPlotter()
	plotter.offsetx = 2.0					# lower left corner of picture */
	plotter.offsety = 2.0
	plotter.scale = float(options.scale)
	plotter.rotate = bool(options.rotate)
	plotter.window = float(options.window) / float(options.magnifya)
	plotter.fillon = fillon
	plotter.ps_header()
	quasi(zfill, midon, symmetry, maxmax, plotter)
	plotter.ps_footer()
	
	exit(0)


# Uncomment main() to use as command line tool compatible with the original
# quasi.c which produces postscript to stdout.
#main()
