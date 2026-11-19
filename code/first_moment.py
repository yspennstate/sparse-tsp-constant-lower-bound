"""First-moment constants quoted in the paper (standard library only).

b_d = (d/e) (v_d d!)^(-1/d)            the first-moment lower bound for c(0+) in dimension d
a_d = e (d-1) (v_d d!)^(-1/(d-1))      the formula of the correction note
g_d = sqrt(d / (2 pi e))               the common limit shape of beta_d (Rhee) and b_d

Also checks the exact count  int_{r_1+..+r_m <= R} prod (d v_d r_i^(d-1)) dr = (v_d d!)^m R^(dm) / (dm)!
against a direct numerical integration for small m, d, and the 27/26 hand bound.

usage: python first_moment.py
"""
import sys
from math import pi, e, lgamma, log, exp, sqrt, factorial


def log_vd(d):
    return (d / 2) * log(pi) - lgamma(d / 2 + 1)


def b(d):
    return (d / e) * exp(-(log_vd(d) + lgamma(d + 1)) / d)


def a_note(d):
    return e * (d - 1) * exp(-(log_vd(d) + lgamma(d + 1)) / (d - 1))


def g(d):
    return sqrt(d / (2 * pi * e))


def simplex_check(d, m, R=1.0, steps=400):
    """Numerical check of the exact count for m = 2 by a midpoint rule in r_1."""
    assert m == 2
    vd = exp(log_vd(d))
    h = R / steps
    tot = 0.0
    for i in range(steps):
        r1 = (i + 0.5) * h
        # inner integral of d v_d r2^(d-1) over [0, R - r1] is v_d (R - r1)^d
        tot += d * vd * r1 ** (d - 1) * vd * (R - r1) ** d * h
    exact = (vd * factorial(d)) ** 2 * R ** (2 * d) / factorial(2 * d)
    return tot, exact


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    for d in (2, 3):
        t, x = simplex_check(d, 2)
        print("exact count, d=%d m=2: midpoint %.6f  formula %.6f" % (d, t, x))
    print("b_2 = sqrt(2/pi)/e = %.6f   (check %.6f)" % (b(2), sqrt(2 / pi) / e))
    print("b_3 = %.6f" % b(3))
    print("a_2 = e/(2 pi) = %.6f   (check %.6f)" % (a_note(2), e / (2 * pi)))
    print("27/26 hand bound: %.6f" % (b(2) * (27 / 26) ** (1 / 6)))
    for d in (2, 3, 10, 100, 1000, 10 ** 4, 10 ** 6):
        print("d=%-8d b_d/g_d = %.6f   a_d/g_d = %.6f   (e^2 = %.6f)" % (d, b(d) / g(d), a_note(d) / g(d), e * e))


if __name__ == "__main__":
    main()
