# Tests of abroaden

`abroaden` (tools/abroaden) broadens the `.abin` files that `nrg` writes with `adaptive_broadening`. The run scripts
generate their input with `make_input.pl` (24 delta peaks of total weight 1, both signs, from 1e-5 to 0.9) and store
no references.

| Test | What it checks |
| --- | --- |
| `abroaden1_conventional` | With σ = lnΛ for every weight, `abroaden` equals `broaden` with the same alpha and `-x`, to 1e-10 relative. |
| `abroaden2_zavg` | Two runs with weights 0.5 and 1.5 times a spectrum average to the spectrum (Eq. 15). |
| `abroaden3_weight` | A spectrum of weight 1 with σ varying from peak to peak broadens to weight 1 within 1e-3. |
| `abroaden4_errors` | Invalid arguments and a missing file fail. |
| `abroaden5_gamma` | `-X`: 41 points in `gamma_scan.dat`, the chosen γ/T in [0.01, 1], and `-x` at the chosen value reproduces the spectrum. |
| `abroaden6_wide` | A few cells with σ = 8 lnΛ next to the usual ones: all cells together equal the sum of the narrow and the wide cells alone (the grouping by kernel width keeps every kernel's own cutoff). |
