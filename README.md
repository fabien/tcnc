tcnc
====

An Inkscape extension that generates G-code suitable for a four axis CNC cutter (XYZ plus tangent A).

This extension is designed to be compatible with a modified Fletcher 6100 cnc mat cutter
controlled by EMC2 (now LinuxCNC).

This is pre-pre-alpha software and hasn't been tested much, so if you need something
that works I would suggest getting the gcodetools plugin for now. Like gcodetools,
tcnc uses biarc approximation to convert cubic Bezier splines into circular arc segments.

However, if you are wondering how to build an Inkscape extension then it might
be helpful to look at the tcnc source code. I suggest taking a look at the Eggbot
extension as well. Documentation about Inkscape extensions is very poor
and the best way to learn how to write one is to look at previous attempts.

Installing tcnc
---------------

I haven't tried it on Windows yet but the idea is to just copy the contents of tcnc/src into
wherever the Inkscape extension folder is.

### For Mac OS X and Linux:
Copy the contents of the tcnc/src folder into ~/.config/inkscape/extensions then
restart Inkscape.
