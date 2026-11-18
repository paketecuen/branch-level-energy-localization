# Branch-Level Energy Localization in Three-Phase Loads

MATLAB code and data of the test cases of the manuscript

> F. G. Montoya, F. de Leon, F. M. Arrabal-Campos and A. Alcayde,
> "Branch-Level Energy Localization in Three-Phase Loads: Resolving
> Indeterminacy in Time-Domain," submitted to *IEEE Transactions on Power
> Delivery*, 2026. Preprint: [arXiv:2606.07076](https://arxiv.org/abs/2606.07076).

F. G. Montoya, F. M. Arrabal-Campos and A. Alcayde are with the Universidad de
Almería, Spain. F. de Leon is with New York University, USA. Corresponding
author: Francisco G. Montoya.

The manuscript decomposes the instantaneous terminal power of a three-phase
load into Joule dissipation and stored-energy rates, branch by branch, for a
selected lumped equivalent. This repository contains the scripts that generate
the waveforms, the branch terms, the balance residuals and the figures of the
single-phase illustration and of the six three-phase test cases. All the cases
are synthetic. No measured data are used.

## Contents

| Path | Content |
|---|---|
| `matlab/` | One script per test case and `run_all_cases.m`. |
| `data/waveforms/` | Waveforms and branch terms of each case, as MAT-files written by the scripts. |
| `data/figures/` | Figures written by the scripts. |
| `data/run_log.txt` | Console output of `run_all_cases`: identified parameters and balance residuals. |
| `data/fourwire_reference/` | Curves of the two four-wire figures of the manuscript, as CSV files. |

## Test cases

| Section | Script | Figure written | Content |
|---|---|---|---|
| II-G | `paradox_czarnecki_CL` | `fig_p2_czarnecki_CL_illustration.pdf` | Single-phase parallel LC branch with zero Budeanu reactive power. |
| IV-A | `paradox_balanced_3ph_RL` | `fig_p2_balanced_3ph_rl_energy.pdf` | Balanced series RL load: constant three-phase power and per-phase exchange. |
| IV-A | `paradox_balanced_RL_compare` | `fig_p2_balanced_RL_compare.pdf` | Branch terms and instantaneous terms of IEEE Std. 1459 for phase a. |
| IV-B | `paradox_resistive_unbalanced` | `fig_p2_resistive_unbalanced_pq.pdf` | Single resistor between two phases, with CPC active currents and p-q coordinates. |
| IV-C | `paradox_nonlinear_switched` | `fig_p2_nonlinear_switched.pdf` | Delta-connected TRIAC-switched resistors. |
| IV-D | `paradox_topology_indeterminacy` | `fig_p2_threewire_indeterminacy.pdf` | Delta-parallel and wye-series equivalents of the same three-wire terminal data. |
| IV-E | `paradox_fluctuating_phase` | `fig_p2_fluctuating_phase.pdf` | Three-phase fluctuating-phase load. |
| IV-F | `fourwire_nonlinear` | `fig_fourwire_powers.pdf`, `fig_fourwire_duality.pdf` | Unbalanced four-wire nonlinear load and series-parallel map of Corollary 2. |

## Running the scripts

The scripts were run with MATLAB R2024a. No toolbox is required.

```matlab
cd matlab
run_all_cases
```

Each script can also be run on its own. The outputs are written to `results/`,
which is not tracked. The files in `data/waveforms/`, `data/figures/` and
`data/run_log.txt` are the outputs of one run of `run_all_cases`.

The relative balance residuals between the terminal power and the sum of the
branch terms reported in `data/run_log.txt` are between 1e-16 and 6e-16.

## Four-wire case (Section IV-F)

The load is connected between each line and the neutral, and the neutral
conductor is measured. The supply is
`v_x = 120*sqrt(2)*sin(w t - theta_x) + 0.12*sqrt(2)*sin(3 w t)` at 50 Hz, with
`theta = (0, 2*pi/3, -2*pi/3)`.

| Phase | Branch | Parameters |
|---|---|---|
| a | Hysteretic inductor in parallel with a capacitor | Capacitor of 500 µF. Piecewise-linear loop with slopes of 10 mH and 1 mH and knees at flux linkages of ±0.15 Wb and ±0.40 Wb. |
| b | Linear parallel RLC | 2 Ω, 10 mH and 2500 µF. |
| c | Series RL with a freewheeling diode, connected through a transistor | 0.75 Ω and 0.4 mH. The transistor conducts from 15.33 ms to 22.33 ms in each period of 20 ms. |

`fourwire_nonlinear` computes the waveforms from this specification. The main
results are the following.

- Phase a contains no resistor. The loss is the area of the hysteresis loop,
  123.75 J per cycle (6187.5 W), and it appears as the cycle average of the
  power of the inductive branch.
- In phase c, the inductor keeps 2.15 J at turn-off. This energy is dissipated
  in the resistor through the diode and does not flow at the terminals.
- The neutral current is 197 A rms.
- The series parameters defined by Corollary 2 for phase b reproduce the three
  power terms of the parallel branch within 6e-12 W. The running integral of
  the series inductive power equals the energy stored in the parallel
  inductor within the error of the numerical integration.

The two four-wire figures of the manuscript were obtained in 2021 from the
parameters identified on the same circuit. Their curves are stored in
`data/fourwire_reference/`, and the script reports the difference with respect
to them. The relative rms difference of the phase currents is 0.6 % for phase
a and below 0.04 % for phases b and c. The curve of the current of phase a in
the manuscript has a mean value of 1.05 A, which the symmetric loop does not
reproduce.

The series elastance `S_s` of Corollary 2 depends on the integration constant
of the primitive of the current. The script uses the zero-mean primitive,
whereas the corresponding panel of the manuscript was computed with a
primitive that includes a large constant. The values of `S_s` are therefore
different, while `R_s`, `L_s` and the three power terms coincide.

### Reference files

| File | Columns |
|---|---|
| `powers_phase_a.csv`, `powers_phase_b.csv`, `powers_phase_c.csv` | Time in ms, voltage, current, active power, power of the inductive branch, power of the capacitive branch, their sum, and terminal power. |
| `duality_phase_b.csv` | Time in ms, series parameters, and the three power terms of the parallel and series representations. |

## Scope

The cases are simulations with known parameters. Measurement noise,
synchronization errors and the reconstruction of derivatives and primitives
from sampled records are not considered here. Code comments are in English.

## License

- Code in `matlab/`: MIT License (`LICENSE`).
- Data in `data/`: Creative Commons Attribution 4.0 International
  (`LICENSE-DATA.md`).

## Citation

```bibtex
@misc{montoya2026branchlevel,
  author        = {Montoya, Francisco G. and de Leon, Francisco and Arrabal-Campos, Francisco M. and Alcayde, Alfredo},
  title         = {Branch-Level Energy Localization in Three-Phase Loads: Resolving Indeterminacy in Time-Domain},
  year          = {2026},
  eprint        = {2606.07076},
  archivePrefix = {arXiv},
  primaryClass  = {eess.SP}
}
```
