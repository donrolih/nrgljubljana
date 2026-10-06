# `abroaden`

Adaptive broadening of NRG spectral data, after S.-S. B. Lee and A. Weichselbaum, Phys. Rev. B 94, 235127 (2016).

`nrg` with `adaptive_broadening=true` writes, next to each spectral function, a file `<name>.abin`: the discrete
spectral weights binned in the frequency ω and in their broadening width σ = |d ln|ω|/dz|, which measures how far a
weight moves when the discretization shift z changes. `abroaden` turns these files into a continuous spectral
function:

1. it reads `i/<name>` for i = 1, …, Nz and averages the weights cell by cell (Eq. 15 of the paper);
2. it broadens each cell with a log-Gaussian (Eq. 17a) of width alpha·σ/lnΛ;
3. it convolves the result over all frequencies with −f′(ω) at the effective temperature γ (Eq. 16 and 21).

The log-Gaussian comes first and the −f′ step second, as the paper requires; the second step replaces the crossover
to a Gaussian at low frequencies that `broaden` uses.

## Usage

```text
abroaden [options] <name> <Nz> <alpha> <T> <Lambda>
```

- `alpha` is the width that a weight with σ = lnΛ gets, i.e. the `alpha` of conventional broadening: in the
  fixed-point regime, where σ = lnΛ, adaptive and conventional broadening coincide. The paper's α corresponds to
  `alpha`·Nz/lnΛ.
- `T` is the temperature of the calculation, `Lambda` the discretization parameter.

| Option | Meaning |
| --- | --- |
| `-x <ratio>` | γ = ratio·T for the −f′ convolution; default 0.1, the value the paper finds typical |
| `-S <sigma_min>` | floor on σ; default none |
| `-X` | choose γ automatically (Sec. V.2); see below |
| `-e <threshold>` | threshold on ε_rms for `-X`; default 0.002 |
| `-o` | read `<name>` directly, for Nz = 1 |
| `-I` | records carry an imaginary part (`reim=true` in `nrg`); the real part is broadened |
| `-m`, `-M`, `-r`, `-L` | the output frequency mesh, as in `broaden` |
| `-j <n>` | worker count (default: `OMP_NUM_THREADS`, or 1) |
| `-v`, `-V`, `-h` | verbose, version, help |

The result is written to `spec.dat` (ω, A(ω)).

## Automatic γ

With `-X`, γ is chosen as in Sec. V.2 of the paper. For γ/T on 41 logarithmically spaced points in [0.01, 1],
the −f′ convolution of the log-broadened curve is evaluated at 11 equally spaced points on [−γΛ², γΛ²], a quadratic is
fitted through them by least squares, and ε_rms (Eq. 24), the standard error of the fit divided by the mean, measures
how smooth the curve is around ω = 0. Small γ leaves the discreteness of the data visible (large ε_rms); large γ
overbroadens real structure (ε_rms rises again).

The chosen γ is the smallest one at which ε_rms falls below the threshold (`-e`, default 0.002), interpolated in
(ln γ, ln ε_rms) between the two scan points around the crossing. If ε_rms is below the threshold already at
γ/T = 0.01, that value is taken; if it never falls below, the minimum of ε_rms is taken, with a warning. The scan
is written to `gamma_scan.dat` (γ/T, ε_rms). The paper's rule for an irregular window between under- and
overbroadening (the geometric mean of the largest underbroadening and the smallest overbroadening γ) is a judgement
on the shape of the scan and is left to the user, who can set γ with `-x`. `-X` requires a mesh reaching below
0.01·T (`-m`) and a positive spectrum near ω = 0.

## Input format

`.abin`: headerless native doubles, records (ω, σ, Re w) or, with `reim=true`, (ω, σ, Re w, Im w). The runs of
different z must bin on the same ω grid, which `nrg` ensures by requiring explicit `emin` and `emax` with
`adaptive_broadening`.
