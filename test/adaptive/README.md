# Tests of adaptive broadening

Adaptive broadening (Lee and Weichselbaum, PRB 94, 235127 (2016)) gives each spectral weight its own width from how
its energy moves with the discretization shift `z`. `nrg` computes dE/dz of every state with `adaptive_broadening`,
from the chain at `z` and at a shifted `z` in the `D` block of `data`.

Each test has its own `run` script and stores no reference: it compares runs with each other.

| Test | What it checks |
| --- | --- |
| `adaptive_hf` | Hellmann–Feynman: dE/dz of every state at every shell against central differences of the energies of runs at z ± dz, with no truncation. |

`adaptive_hf` builds `data` like `test/mixchain+nrgrun/mixchain_qs_vs_nrgchain`: the flat-band chain from `nrgchain`
at z = 0.5 and 0.5 ± 1e-4, the SIAM/QS template of `test/mixchain+nrgrun/templates/QS` and its `instantiate`, plus a
`D` block with a neighbouring chain. All three runs dump `dEdz.nrg` (`adaptive_broadening_dumpenergies`), which holds
the energies and dE/dz at full precision, and `compare_hf.pl` compares the middle run with the differences.

The energies are in units of the shell, which depends on z, so the expected value is
`-ln(Lambda) (E_i - E_GS) + [(E_i - E_GS)(z+dz) - (E_i - E_GS)(z-dz)] / (2 dz)`. The tolerance, 1e-3 relative, covers
the one-sided derivative of the chain in the `D` block (an error of order dz) and the central difference (dz^2).

Near an avoided crossing inside [z-dz, z+dz] the difference quotient of a single level means nothing, while
Hellmann–Feynman still holds; such pairs were seen with gaps of 2e-5 at z and 2e-4 at z ± dz. States of a subspace
whose neighbouring levels come closer than 1e-3 in any of the three runs are therefore compared through their sum,
which is smooth through the crossing (the trace of dH/dz over the group). On the first run this formed 544 groups among the
8774 states; the largest deviation was 3.7e-4.
