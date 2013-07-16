from lib import geom
p = (
    geom.P(2066.4118, 1986.3414),
    geom.P(1995.2052, 1674.3645),

    geom.P(1995.2052, 1674.3645),
    geom.P(2066.4118, 1362.3876),

    geom.P(2066.4118, 1362.3876),
    geom.P(2137.6185, 1674.3645),

    geom.P(2137.6185, 1674.3645),
    geom.P(2066.4118, 1986.3414),

    geom.P(2066.4118, 1986.3414),
    geom.P(1995.2052, 1674.3645),

    geom.P(1995.2052, 1674.3645),
    geom.P(2066.4118, 1362.3876),

    geom.P(2066.4118, 1362.3876),
    geom.P(2137.6185, 1674.3645),

    geom.P(2137.6185, 1674.3645),
    geom.P(2066.4118, 1986.3414),

    geom.P(2066.4118, 1986.3414),
    geom.P(2386.4118, 1986.3414),

    geom.P(2386.4118, 1986.3414),
    geom.P(2457.6185, 1674.3645),

    geom.P(2457.6185, 1674.3645),
    geom.P(2137.6185, 1674.3645),
)
L = (
    geom.Line(p[0], p[1]),
    geom.Line(p[2], p[3]),
    geom.Line(p[4], p[5]),
    geom.Line(p[6], p[7]),
    geom.Line(p[8], p[9]),
    geom.Line(p[10], p[11]),
    geom.Line(p[12], p[13]),
    geom.Line(p[14], p[15]),
    geom.Line(p[16], p[17]),
    geom.Line(p[18], p[19]),
    geom.Line(p[20], p[21]),
)

p0 = geom.P(0, 0)
p1 = geom.P(1, 0)
p2 = geom.P(1, 1)
p3 = geom.P(0, 1)
p4 = geom.P(-1, 1)
p5 = geom.P(-1, 0)
p6 = geom.P(-1, -1)
p7 = geom.P(0, -1)
p8 = geom.P(1, -1)