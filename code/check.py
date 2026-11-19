"""Independent check of the certificate written by cert.cpp (Python 3 with numpy).

1. Inputs.  Re-verifies grid.txt with enclosures that share nothing with grid.py: exp(-x) by its alternating
   series, pi by Euler's formula pi = 4 arctan(1/2) + 4 arctan(1/3), and cos x = 1 - 2 sin^2(x/2) with a lower
   bound for sin.  Checks 0 <= w_j <= Q P(R in bin j), wildcard = Q - sum w_j > 0, and C_t >= D max cos on bin t.
2. Image.  Recomputes Mv from v without the threshold structure used in cert.cpp: for every ordinary cell the
   swap test is decided in floating point, every comparison closer than 1e-3 to a tie is redone exactly in
   integers in the form  L^2 > X + Y  and  (L^2 + Y - X)^2 > 4 L^2 Y,  and the sums over (a, t) are matrix-vector
   products on 20-bit halves of v (exact in float64).  The result must equal the saved Mv entry by entry.
3. Inequality.  10^12 (Mv)_i <= N Q A v_i for every state i, with Python integers.
4. Constant.  c(0+) >= sqrt(2/pi)/e * rho^(-1/2) > 0.32417 via 2 > 0.32417^2 * pi_hi * e_hi^2 * rho.

usage: python check.py grid.txt prefix N
"""
import sys
from fractions import Fraction as F
from math import factorial
import numpy as np


def read_grid(path):
    with open(path) as f:
        n, s, A, Q, D = map(int, f.readline().split())
        w = list(map(int, f.readline().split()))
        wstar = int(f.readline())
        C = list(map(int, f.readline().split()))
    assert len(w) == n and len(C) == A
    return n, s, A, Q, D, w, wstar, C


def expneg_bounds(x, N=160):
    """exp(-x) for 0 <= x <= 8 by the alternating series; terms decrease for k > x, |remainder| <= next term."""
    terms = [F((-1) ** k) * x ** k / factorial(k) for k in range(N + 2)]
    s = sum(terms[: N + 1])
    r = abs(terms[N + 1])
    return s - r, s + r


def atan_bounds(y, K=80):
    def partial(m):
        return sum(F((-1) ** k) * y ** (2 * k + 1) / (2 * k + 1) for k in range(m))
    return partial(2 * K), partial(2 * K + 1)


def pi_bounds_euler():
    a = atan_bounds(F(1, 2))
    b = atan_bounds(F(1, 3))
    lo, hi = 4 * (a[0] + b[0]), 4 * (a[1] + b[1])
    S = 10 ** 40   # round outward to 40 decimals
    return F((lo * S).numerator // (lo * S).denominator, S), F(-((-hi * S).numerator // (hi * S).denominator), S)


def sin_lower(y, K=30):
    """0 <= y <= 2: partial sum ending in a negative term is a lower bound (terms decrease from the first)."""
    return sum(F((-1) ** k) * y ** (2 * k + 1) / factorial(2 * k + 1) for k in range(2 * K))


def check_inputs(n, s, A, Q, D, w, wstar, C):
    assert all(x >= 0 for x in w) and wstar > 0 and sum(w) + wstar == Q
    slack = None
    for j in range(n):
        x0, x1 = F(j, s), F(j + 1, s)
        lo0, _ = expneg_bounds(x0)
        _, hi1 = expneg_bounds(x1)
        mass_lo = (1 + x0) * lo0 - (1 + x1) * hi1
        assert w[j] <= Q * mass_lo, j
        sl = Q * mass_lo - w[j]
        slack = sl if slack is None or sl < slack else slack
    pi_lo, pi_hi = pi_bounds_euler()
    assert F(314159, 100000) < pi_lo < pi_hi < F(314160, 100000)
    for t in range(A):
        j = min(t, A - 1 - t)
        y = pi_lo * j / A             # half of the lower bound for the angle 2 pi j / A
        sl = sin_lower(y)
        assert sl >= 0
        cos_hi = 1 - 2 * sl * sl
        assert C[t] >= D * cos_hi, t
    print("inputs ok: %d radius masses (least slack %.3g in units of 1/Q), %d cosine bounds" % (n, float(slack), A))
    return pi_hi


def exact_forb(a, c, z, Ct, Cu, D):
    X = D * (a * a + z * z) + 2 * a * z * Ct
    Y = D * (c * c + z * z) + 2 * c * z * Cu
    L2 = D * (a + c) ** 2
    if not L2 > X + Y:
        return False
    return (L2 + Y - X) ** 2 > 4 * L2 * Y


def forbidden_at(z, n, A, D, Carr):
    """Boolean array F[a, c, t, u] (ordinary a, c) for the swap test at middle radius z."""
    if z == 0:
        # X = D a^2, Y = D c^2, L^2 = D (a + c)^2: then T = 2 D a c and T^2 = 4 X Y, so no cell is forbidden
        return np.zeros((n, n, A, A), dtype=bool), 0
    a = np.arange(n, dtype=np.float64)
    X = D * (a[:, None] ** 2 + z * z) + 2.0 * a[:, None] * z * Carr[None, :]   # (a, t)
    sx = np.sqrt(X)
    L = np.sqrt(float(D)) * (a[:, None] + a[None, :])                            # (a, c)
    g = sx[:, None, :, None] + sx[None, :, None, :] - L[:, :, None, None]
    Fz = g < 0
    amb = np.argwhere(np.abs(g) < 1e-3)
    for (ia, ic, it, iu) in amb:
        Fz[ia, ic, it, iu] = exact_forb(int(ia), int(ic), z, int(Carr[it]), int(Carr[iu]), D)
    return Fz, len(amb)


def main():
    gpath, prefix, N = sys.argv[1], sys.argv[2], int(sys.argv[3])
    n, s, A, Q, D, w, wstar, C = read_grid(gpath)
    pi_hi = check_inputs(n, s, A, Q, D, w, wstar, C)
    N1 = n + 1
    v = np.loadtxt(prefix + ".v", dtype=np.int64).reshape(N1, N1, A)
    Mv_saved = [int(x) for x in open(prefix + ".Mv")]
    assert v.min() >= 1 and v.max() < 2 ** 41
    W = w + [wstar]
    Carr = np.array(C, dtype=np.float64)
    Mv = [0] * (N1 * N1 * A)
    namb = 0
    F_next, k = forbidden_at(0, n, A, D, Carr)
    namb += k
    for b in range(N1):
        R = int(v[:, b, :].sum())
        if b < n:
            F_lo = F_next
            F_next, k = forbidden_at(b + 1, n, A, D, Carr)
            namb += k
            cell = (F_lo & F_next).transpose(0, 2, 1, 3).reshape(n * A, n * A).astype(np.float64)  # (a t) x (c u)
            vb = v[:n, b, :].reshape(n * A)
            lo = (vb & ((1 << 20) - 1)).astype(np.float64)
            hi = (vb >> 20).astype(np.float64)
            f_lo = (lo @ cell).astype(np.int64)
            f_hi = (hi @ cell).astype(np.int64)
            forb = (f_hi << 20) + f_lo                      # forbidden mass, indexed by (c, u)
            forb = forb.reshape(n, A)
        for c in range(N1):
            for u in range(A):
                Sbcu = R if (b == n or c == n) else R - int(forb[c, u])
                Mv[(b * N1 + c) * A + u] = W[c] * Sbcu
        if b % 16 == 0:
            print("b=%d done" % b, flush=True)
    bad = sum(1 for i in range(len(Mv)) if Mv[i] != Mv_saved[i])
    print("ambiguous float comparisons redone exactly: %d" % namb)
    print("image mismatches against the saved Mv: %d of %d" % (bad, len(Mv)))
    assert bad == 0
    vf = v.reshape(-1)
    QA = Q * A
    for i in range(len(Mv)):
        assert 10 ** 12 * Mv[i] <= N * QA * int(vf[i]), i
    print("all %d inequalities 10^12 (Mv)_i <= %d Q A v_i hold" % (len(Mv), N))
    s_e = sum(F(1, factorial(k)) for k in range(40))
    e_hi = s_e + F(2, factorial(40))          # tail sum_{k>=40} 1/k! < 2/40!
    rho = F(N, 10 ** 12)
    lhs = F(32417, 100000) ** 2 * pi_hi * e_hi ** 2 * rho
    assert 2 > lhs
    print("2 > 0.32417^2 pi e^2 rho holds; so c(0+) > 0.32417  (slack %.4g)" % float(2 - lhs))


if __name__ == "__main__":
    main()
