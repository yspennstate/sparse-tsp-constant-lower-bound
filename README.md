# Lower bounds for the sparse travelling salesman constant of Aldous and Krikun

Paper and programs. Let T_m be the length of the shortest path from the origin through m distinct points of a Poisson
process of unit intensity in R^d. Aldous and Krikun proved that T_m/m converges in probability to a constant c(0+).
The exponential tilt of Aldous and Krikun, with its integral evaluated exactly, gives
c(0+) >= b_d = (d/e)(v_d Gamma(d+1))^(-1/d), which is 0.29352... in the plane, a value that also follows from the proof of a theorem of Carlsson. The paper shows that the expected number
of paths behind this bound tends to infinity for every larger threshold, and that c(0+)/sqrt(d) -> (2 pi e)^(-1/2).
Counting only paths that cannot be shortened by exchanging two consecutive points, it proves c(0+) > 0.32417 in the
plane, with a transfer operator and a positive vector checked in exact integer arithmetic.

- `main.tex`, `main.pdf`: the paper.
- `code/grid.py`: writes `grid.txt`, the integer length weights and cosine bounds (standard library only).
- `code/cert.cpp`: finds the positive vector v, computes Mv exactly and prints the least numerator N.
- `code/check.py`: independent check of grid.txt, v, Mv, the inequalities and the final constant (needs numpy).
- `code/first_moment.py`: the first moment constants of Section 2.
- `code/SHA256SUMS`: hashes of grid.txt, cert.cpp, check.py, grid.py and the two vectors.
- `LICENSE`: CC BY 4.0 for the paper, MIT for the programs.

The vectors `final.v` (7.1 MB) and `final.Mv` (13.9 MB), one integer per line in the state order
(a*(n+1)+b)*A+t, are in `code/`.

## Reproducing the vector v and the check

```
cd code
python grid.py                       # rewrites grid.txt; compare with SHA256SUMS
clang++ -O2 -std=c++17 cert.cpp -o cert
./cert grid.txt final 400            # prints N 819869595611
python check.py grid.txt final 819869595611
```

`check.py` takes about four minutes on a laptop; `outputs/check.txt` is its output, and `outputs/check_neg.txt` the control below. To check the saved vectors without running `cert`,
verify them in `code/` with `sha256sum -c SHA256SUMS` and run the last command.

Controls: `python check.py grid.txt final 819869595610` stops with an assertion at
state 590292, and a copy of `final.v` with one entry increased by one (with `final.Mv` unchanged) gives 756 image
mismatches.
