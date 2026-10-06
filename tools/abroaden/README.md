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
| `-o` | read `<name>` directly, for Nz = 1 |
| `-I` | records carry an imaginary part (`reim=true` in `nrg`); the real part is broadened |
| `-m`, `-M`, `-r`, `-L` | the output frequency mesh, as in `broaden` |
| `-j <n>` | worker count (default: `OMP_NUM_THREADS`, or 1) |
| `-v`, `-V`, `-h` | verbose, version, help |

The result is written to `spec.dat` (ω, A(ω)).

## Input format

`.abin`: headerless native doubles, records (ω, σ, Re w) or, with `reim=true`, (ω, σ, Re w, Im w). The runs of
different z must bin on the same ω grid, which `nrg` ensures by requiring explicit `emin` and `emax` with
`adaptive_broadening`.
