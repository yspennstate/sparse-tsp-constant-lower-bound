# Lower bounds for the sparse travelling salesman constant of Aldous and Krikun

Paper and programs. Let T_m be the length of the shortest path from the origin through m distinct points of a Poisson
process of unit intensity in R^d. Aldous and Krikun proved that T_m/m converges in probability to a constant c(0+).
The paper shows that the first moment method gives c(0+) >= (d/e)(v_d Gamma(d+1))^(-1/d), that no first moment count gives
more, and that c(0+)/sqrt(d) -> (2 pi e)^(-1/2). It also proves c(0+) > 0.32417 in the plane with a transfer
certificate checked in exact integer arithmetic.

- `main.tex`, `main.pdf`: the paper.
- `code/grid.py`: writes `grid.txt`, the integer radius weights and cosine bounds (standard library only).
- `code/cert.cpp`: finds the certificate vector v, computes Mv exactly and prints the least numerator N.
- `code/check.py`: independent check of grid.txt, v, Mv, the inequalities and the final constant (needs numpy).
- `code/first_moment.py`: the first moment constants of Section 2.
- `code/SHA256SUMS`: hashes of grid.txt, cert.cpp, check.py, grid.py and the two vectors.

The vectors `final.v` (7.1 MB) and `final.Mv` (13.9 MB), one integer per line in the state order
(a*(n+1)+b)*A+t, are attached to the release v1.0.

## Reproducing the certificate

```
cd code
python grid.py                       # rewrites grid.txt; compare with SHA256SUMS
clang++ -O2 -std=c++17 cert.cpp -o cert
./cert grid.txt final 400            # about a minute on a Mac mini; prints N 819869595611
python check.py grid.txt final 819869595611
```

`check.py` also takes about a minute. To check the released vectors without running `cert`,
download `final.v` and `final.Mv` from the release into `code/`, verify them with `sha256sum -c SHA256SUMS`, and run
the last command.

Two controls show that the check can fail: `python check.py grid.txt final 819869595610` stops with an assertion at
state 590292, and a copy of `final.v` with one entry increased by one (with `final.Mv` unchanged) gives 756 image
mismatches.
