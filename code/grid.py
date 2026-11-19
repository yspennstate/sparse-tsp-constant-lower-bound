"""Exact rational inputs for the transfer certificate (standard library only).

Writes grid.txt: the radius-bin weights w_j and the cosine upper bounds C_t used by cert.cpp and check.py.

Radii.  Under the density exp(-|u|)/(2 pi) on R^2 the radius has density r exp(-r), so
P(R in [x, y]) = S(x) - S(y) with S(x) = (1 + x) exp(-x).  Bin j is [j/s, (j+1)/s], j = 0..n-1, and
w_j = floor(Q * (S_lo(j/s) - S_hi((j+1)/s))) <= Q * P(R in bin j).  The wildcard weight is Q - sum w_j.
exp(x) is enclosed by its Taylor polynomial of degree N and the tail bound x^(N+1)/(N+1)! / (1 - x/(N+2)).

Angles.  Bin t is [2 pi t/A, 2 pi (t+1)/A].  The largest value of cos on it is cos(2 pi j/A) with
j = min(t, A-1-t), and 0 <= 2 pi j/A < pi.  Cos decreases on [0, pi], so cos(2 pi j/A) <= cos(2 pi_lo j/A),
where pi_lo < pi comes from Machin's formula.  The cosine series is cut after the term of degree 2K and the
absolute value of the next term is added (the terms decrease from degree 2 on when |x| <= pi).
C_t = ceil(D * upper bound) >= D * max cos on bin t.

Usage: python grid.py [n s A]   (defaults 96 12 64; Q = 2^32, D = 2^30)
"""
import sys
from fractions import Fraction as F
from math import factorial

Q = 2 ** 32
D = 2 ** 30


def exp_bounds(x, N=80):
    """Lower and upper rational bounds for exp(x), 0 <= x <= 8."""
    assert 0 <= x <= 8
    term, s = F(1), F(1)
    for k in range(1, N + 1):
        term = term * x / k
        s += term
    nxt = term * x / (N + 1)
    return s, s + nxt / (1 - F(x) / (N + 2))


def S_bounds(x):
    lo, hi = exp_bounds(x)
    return (1 + x) / hi, (1 + x) / lo


def atan_bounds(y, K=60):
    """Alternating series with decreasing terms: an even number of terms is a lower bound, odd an upper bound."""
    def partial(m):
        return sum(F((-1) ** k) * y ** (2 * k + 1) / (2 * k + 1) for k in range(m))
    return partial(2 * K), partial(2 * K + 1)


def pi_bounds():
    a_lo, a_hi = atan_bounds(F(1, 5))
    b_lo, b_hi = atan_bounds(F(1, 239))
    lo, hi = 16 * a_lo - 4 * b_hi, 16 * a_hi - 4 * b_lo
    # round outward to 40 decimals so later powers stay small; still lower and upper bounds
    S = 10 ** 40
    return F((lo * S).numerator // (lo * S).denominator, S), F(-((-hi * S).numerator // (hi * S).denominator), S)


def cos_upper(x, K=40):
    assert 0 <= x <= 4
    s = sum(F((-1) ** k) * x ** (2 * k) / factorial(2 * k) for k in range(K + 1))
    return s + x ** (2 * K + 2) / factorial(2 * K + 2)


def ceil_frac(fr):
    return -((-fr.numerator) // fr.denominator)


def main():
    n, s, A = (int(a) for a in sys.argv[1:4]) if len(sys.argv) >= 4 else (96, 12, 64)
    pi_lo, pi_hi = pi_bounds()
    assert pi_lo < pi_hi and pi_hi - pi_lo < F(3, 10 ** 40)
    w = []
    for j in range(n):
        lo_j, _ = S_bounds(F(j, s))
        _, hi_j1 = S_bounds(F(j + 1, s))
        m = lo_j - hi_j1
        w.append(max(0, (Q * m).numerator // (Q * m).denominator) if m > 0 else 0)
    wstar = Q - sum(w)
    assert wstar > 0
    C = []
    for t in range(A):
        j = min(t, A - 1 - t)
        x = 2 * pi_lo * j / A
        assert 0 <= x < pi_lo
        C.append(ceil_frac(D * cos_upper(x)))
    assert all(-D <= c <= D for c in C)
    with open("grid.txt", "w") as f:
        f.write("%d %d %d %d %d\n" % (n, s, A, Q, D))
        f.write(" ".join(map(str, w)) + "\n")
        f.write("%d\n" % wstar)
        f.write(" ".join(map(str, C)) + "\n")
    print("n=%d s=%d A=%d  sum w=%d  wildcard=%d (%.6g of Q)  min w=%d" % (n, s, A, sum(w), wstar, wstar / Q, min(w)))
    print("C[0..3] =", C[:4], " C[31..32] =", C[31:33])


if __name__ == "__main__":
    main()
