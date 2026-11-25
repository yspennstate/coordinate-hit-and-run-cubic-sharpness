# Cubic dimension dependence in coordinate hit-and-run is sharp

Yitzchak Shmalo, 25 November 2026.

Kook and Vempala (arXiv:2608.16878, Theorem 1.1) proved that for every convex body K in R^n containing a unit ball, the spectral gap of coordinate hit-and-run is at least c / (n^3 C_PI), where C_PI is the Poincaré constant of the uniform distribution on K. The resulting mixing bound of O(n^3 C_PI log(M/eps)) steps from a start at chi-square divergence M is, in their words, "conjectured to be asymptotically tight". The paper shows that the factor n^3 cannot be lowered.

For every n >= 2 it gives an explicit centrally symmetric convex polytope K_n containing the closed unit ball, a double cone over a box rotated through 45 degrees in two coordinates, with

- 4n^3 / (3(n+2)) <= C_PI(Unif(K_n)) <= 16(n^2+n-1) / pi^2,
- gap(K_n) <= 12 / (n^3 (n^2-1)),
- gap(K_n) C_PI(Unif(K_n)) <= 320 / (pi^2 n^3).

So C_PI = Theta(n^2), the gap is Theta(n^-5), and their product is Theta(n^-3).

## Files

- `main.tex`, `main.pdf`: the paper.
- `CREDITS.md`: which statements are known and which are new, and whose work they rest on.
- `code/verify.py`: checks the rational identities of the proof for n = 2, ..., 100, checks 798 boundary fibres, and estimates the Rayleigh quotient of the linear test function by independent uniform sampling. Needs NumPy. Output: `code/verification.json`.
- `code/exact_chords.py`: computes the coordinate chord second moments in exact rational arithmetic and checks the Rayleigh quotient against its upper bound. Standard library only. Output: `code/exact_chords.json`.

## License

The paper is licensed CC BY 4.0 and the code MIT.
