tcnc
====

An Inkscape extension that generates G-code suitable for a four axis CNC machine
controlled by LinuxCNC v2.4+. Where one axis (A) rotates about the Z axis and is
kept tangent to movement along the X and Y axis.

It is currenty used to produce output for a painting machine based on
a modified Fletcher 6100 CNC mat cutter controlled by EMC2 (now LinuxCNC).

Like gcodetools, tcnc uses biarc approximation to convert cubic Bezier splines
into circular arc segments. Gcodetools also has a tangent knife mode but it has
bugs and the code is very difficult to modify due to it's fairly opaque coding style.

If you are wondering how to build an Inkscape extension then it might
be helpful to look at this source code. There are some reusable components
as well that make writing extensions a little easier.
I suggest taking a look at the Eggbot extension as well.
Documentation about Inkscape extensions is poor
and the best way to learn how to write one is to look at previous attempts.

Installing tcnc
---------------
### For Mac OS X and Linux:
Copy the contents of the tcnc/src folder into ~/.config/inkscape/extensions then
restart Inkscape.

### Windows
Copy the contents of tcnc/src into wherever the Inkscape extension folder is.

Notes
-----
The tcnc/lib folder just contains copies of the Python modules that ship with
Inkscape. This makes it a little easier to develop locally using Eclipse or similar
IDE. Don't copy these to your Inkscape extensions folder - they should already
be there.

