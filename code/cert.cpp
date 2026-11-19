// Transfer certificate for the overlapping exchange event (Section 5 of the paper).
//
// State (a, b, t): radius bins of two consecutive increments (n = wildcard) and the angle bin of their
// relative angle.  M acts by
//     (Mv)(b, c, u) = w_c * sum_{a, t} H(a, b, c, t, u) v(a, b, t),
// where H = 0 exactly when a, b, c are ordinary bins and the swap test below succeeds at both endpoints
// z = b, b+1 of the middle radius; otherwise H = 1.
//
// Swap test (radii in units of 1/s, cosines scaled by D = 2^30, L = 2^15 (a + c) = sqrt(D) (a + c)):
//     X = D(a^2 + z^2) + 2 a z C_t,  Y = D(c^2 + z^2) + 2 c z C_u,  T = L^2 - X - Y,
//     forbidden at z  <=>  T > 0 and T^2 > 4 X Y   (<=> sqrt X + sqrt Y < L).
//
// The program finds a positive integer vector v by power iteration in floating point (any positive v is
// admissible), computes Mv exactly in unsigned 128-bit arithmetic, and reports the least numerator N with
//     10^12 (Mv)_i <= N * Q * A * v_i  for every i,
// so that T v <= rho v with T = M/(QA) and rho = N/10^12.
//
// usage: cert grid.txt out_prefix iterations
// writes out_prefix.v, out_prefix.Mv (one integer per line, state order (a*(n+1)+b)*A+t) and prints N.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

typedef __int128 i128;
typedef unsigned __int128 u128;

static int n, s, A;
static long long Q, D;
static std::vector<long long> W;  // size n+1, W[n] = wildcard weight
static std::vector<long long> C;  // size A

static bool forb_z(long long a, long long c, long long z, long long Ct, long long Cu) {
    i128 X = (i128)D * (a * a + z * z) + (i128)2 * a * z * Ct;
    i128 Y = (i128)D * (c * c + z * z) + (i128)2 * c * z * Cu;
    if (X < 0 || Y < 0) { fprintf(stderr, "negative X or Y\n"); exit(2); }
    i128 L = (i128)32768 * (a + c);
    i128 T = L * L - X - Y;
    if (T <= 0) return false;
    return T * T > (i128)4 * X * Y;
}

static std::string u128s(u128 x) {
    if (x == 0) return "0";
    std::string r;
    while (x > 0) { r.push_back('0' + (int)(x % 10)); x /= 10; }
    std::reverse(r.begin(), r.end());
    return r;
}

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: cert grid.txt out_prefix iterations\n"); return 1; }
    FILE* f = fopen(argv[1], "r");
    if (!f) { perror("grid"); return 1; }
    if (fscanf(f, "%d %d %d %lld %lld", &n, &s, &A, &Q, &D) != 5) return 1;
    if (D != (1LL << 30)) { fprintf(stderr, "D must be 2^30\n"); return 1; }
    W.assign(n + 1, 0);
    long long sum = 0;
    for (int j = 0; j <= n; j++) { if (fscanf(f, "%lld", &W[j]) != 1) return 1; sum += W[j]; }
    if (sum != Q) { fprintf(stderr, "weights do not sum to Q\n"); return 1; }
    C.assign(A, 0);
    for (int t = 0; t < A; t++) if (fscanf(f, "%lld", &C[t]) != 1) return 1;
    fclose(f);
    const int iters = atoi(argv[3]);
    const int N1 = n + 1;
    const size_t NS = (size_t)N1 * N1 * A;

    // fold index of each angle bin; C is non-increasing in the fold index
    std::vector<int> fold(A);
    for (int t = 0; t < A; t++) fold[t] = std::min(t, A - 1 - t);
    const int H2 = A / 2;
    std::vector<long long> Cf(H2);
    for (int j = 0; j < H2; j++) Cf[j] = C[j];
    for (int j = 1; j < H2; j++) if (Cf[j] > Cf[j - 1]) { fprintf(stderr, "C not monotone\n"); return 1; }

    // J[(b*n + c)*n*A + a*A + t] = least fold index j such that every u with fold(u) >= j is forbidden
    // for the ordinary cell (a, b, c, t, u); H2 means no u is forbidden.  Monotone because X, Y grow with C.
    std::vector<uint8_t> J((size_t)n * n * n * A);
    long long nforb = 0;
    for (int b = 0; b < n; b++)
        for (int c = 0; c < n; c++)
            for (int a = 0; a < n; a++)
                for (int t = 0; t < A; t++) {
                    int j = H2;
                    while (j > 0 && forb_z(a, c, b, C[t], Cf[j - 1]) && forb_z(a, c, b + 1, C[t], Cf[j - 1])) j--;
                    J[(((size_t)b * n + c) * n + a) * A + t] = (uint8_t)j;
                    nforb += 2 * (H2 - j);
                }
    fprintf(stderr, "forbidden ordinary cells: %lld\n", nforb);

    auto idx = [&](int a, int b, int t) { return ((size_t)a * N1 + b) * A + t; };

    // generic application of M with accumulator type Acc
    auto applyM = [&](auto const& v, auto& out, auto zero) {
        typedef decltype(zero) Acc;
        std::vector<Acc> bucket(H2 + 1);
        for (int b = 0; b <= n; b++) {
            Acc R = zero;
            for (int a = 0; a <= n; a++) for (int t = 0; t < A; t++) R += (Acc)v[idx(a, b, t)];
            for (int c = 0; c <= n; c++) {
                if (b == n || c == n) {
                    for (int u = 0; u < A; u++) out[idx(b, c, u)] = (Acc)W[c] * R;
                    continue;
                }
                std::fill(bucket.begin(), bucket.end(), zero);
                for (int a = 0; a < n; a++) {
                    const uint8_t* Jr = &J[(((size_t)b * n + c) * n + a) * A];
                    for (int t = 0; t < A; t++) bucket[Jr[t]] += (Acc)v[idx(a, b, t)];
                }
                // forbidden mass for fold index j: sum of buckets with threshold <= j
                std::vector<Acc> pref(H2);
                Acc run = zero;
                for (int j = 0; j < H2; j++) { run += bucket[j]; pref[j] = run; }
                for (int u = 0; u < A; u++) out[idx(b, c, u)] = (Acc)W[c] * (R - pref[fold[u]]);
            }
        }
    };

    std::vector<double> x(NS, 1.0), y(NS);
    double lam = 0;
    for (int it = 0; it < iters; it++) {
        applyM(x, y, 0.0);
        double mx = 0;
        for (size_t i = 0; i < NS; i++) mx = std::max(mx, y[i]);
        double mr = 0;
        for (size_t i = 0; i < NS; i++) mr = std::max(mr, y[i] / x[i]);
        lam = mr / ((double)Q * A);
        for (size_t i = 0; i < NS; i++) x[i] = y[i] / mx;
        if (it % 10 == 0 || it == iters - 1) fprintf(stderr, "iter %d  max ratio %.15f\n", it, lam);
    }
    // integer vector: scale 2^40, round down, add one
    const double S40 = 1099511627776.0;
    std::vector<u128> v(NS), Mv(NS);
    for (size_t i = 0; i < NS; i++) v[i] = (u128)(unsigned long long)std::floor(S40 * x[i]) + 1;
    applyM(v, Mv, (u128)0);
    const u128 DEN = (u128)1000000000000ULL;
    u128 N = 0;
    size_t argmax = 0;
    for (size_t i = 0; i < NS; i++) {
        u128 lhs = DEN * Mv[i];
        u128 rhs1 = (u128)Q * A * v[i];
        u128 q = (lhs + rhs1 - 1) / rhs1;
        if (q > N) { N = q; argmax = i; }
    }
    FILE* fv = fopen((std::string(argv[2]) + ".v").c_str(), "w");
    FILE* fm = fopen((std::string(argv[2]) + ".Mv").c_str(), "w");
    unsigned long long vmin = ~0ULL, vmax = 0;
    for (size_t i = 0; i < NS; i++) {
        fprintf(fv, "%s\n", u128s(v[i]).c_str());
        fprintf(fm, "%s\n", u128s(Mv[i]).c_str());
        vmin = std::min(vmin, (unsigned long long)v[i]);
        vmax = std::max(vmax, (unsigned long long)v[i]);
    }
    fclose(fv); fclose(fm);
    printf("states %zu  v in [%llu, %llu]  floating ratio %.15f\n", NS, vmin, vmax, lam);
    printf("N %s  (rho = N/10^12, attained at state %zu)\n", u128s(N).c_str(), argmax);
    return 0;
}
