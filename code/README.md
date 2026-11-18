# Verification for cubic dimension dependence in coordinate hit-and-run

`verify.py` checks the rational identities in the proof for dimensions 2 through 100, checks 798 boundary fibers, and samples 30,000 independent uniform points in each of dimensions 2, 3, 4, 8, 16, 32, 64 and 128. It requires Python 3 and NumPy. It writes `verification.json` in the working directory. The saved run uses NumPy 2.2.3.

Run from this directory:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 verify.py
python3 exact_chords.py
```

The sampler first draws the slice radius with density `n*r**(n-1)`, then samples the corresponding box. Samples are independent; no Markov-chain burn-in is used. Seeds are `410000+n`. Floating-point estimates concern the Rayleigh quotient of the linear function `t`, which bounds the spectral gap from above. They do not estimate the spectral gap itself. The theorem is an all-dimensional written proof, independent of these finite checks.

`exact_chords.py` integrates coordinate chord lengths with rational arithmetic and checks the actual Rayleigh quotient against its upper envelope, including two-dimensional, equal-scale and short-axis cases. It uses only the Python standard library and writes `exact_chords.json`. Its polygon-integration functions are adapted from the separate verification of the construction.

In the coordinates `(y,t,z)`, the observable `t` has eigenvalue `1-1/n`. The corresponding control entry in `verification.json` records its separation from 1; it does not identify the gap of the entire chain.

The manuscript is `../main.tex`. Compile it from its directory with `pdflatex -interaction=nonstopmode -halt-on-error main.tex`, twice. The repository is https://github.com/yspennstate/coordinate-hit-and-run-cubic-sharpness.
