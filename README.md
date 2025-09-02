# GRH Certificate Verifier (finite, rigorous mode)

This repository contains a reference Python implementation of the Part III “03” verifier algorithm for GRH-style certificates. It follows the paper’s phase-by-phase flow (Phases 0–6) and exposes the ledger/routine primitives described in the sketch: ArchBlock, RamBlock, PrimeBlockBL, PrimeTailBound (heat), RS_Positivity_LB, SANorm, RStarBudget, Q_of_X_eps, Eps_from_q_X, CheckA5Gap.

The core implementation lives in `grh_verifier.py` and is designed to be plug-compatible with external infrastructure that supplies interval-valued contributions for the archimedean, ramified, and prime-local terms.

## Usage

- Requirements: Python 3.9+ (no external dependencies).
- Quick demo: runs a small synthetic certificate
  - `python3 grh_verifier.py`
- Programmatic use:

```python
from decimal import Decimal
from grh_verifier import Interval, TestFunction, PrimeLocalTerm, RamifiedTerm, Certificate, verify_certificate

Phi = TestFunction(family="heat", a=Decimal("0.8"))
cert = Certificate(
    m=2, K="Q", Qpi=Decimal("1000"), t_star=Decimal("1.0"), Phi_test=Phi,
    precision_bits=200, test_family="heat",
    prime_block_mode="rs_lower_bound",  # or "evaluate"
    proof_mode="ann_odd_only",
    # Geometric side placeholders (replace with rigorous intervals)
    arch_value=Interval.point("0.05"),
    ram_terms=[RamifiedTerm(j=1, value=Interval.point("0.01"))], j_max=1,
    # Amplifier / cutoff
    q_modulus=97, band_limit_X=Decimal("6.0"),
    # A5 gap + invariance
    eta_a5_gap=Interval.point("0.40"),
    epsilon_invariance=Interval(Decimal("0"), Decimal("0.20")),
    alpha_inv=Decimal("0.5"), kappa_inv=Decimal("1.0"),
)
report = verify_certificate(cert)
print(report.result)
print(report.details)
```

## Mapping to Phases (0–6)

- Phase 0 — Numerical Rigor: `set_precision_bits` sets Decimal precision from `precision_bits`.
- Phase 1 — Initialization & Admissibility: `TestFunction.check_admissible`; `SANorm` returns an interval bound for `||Phi||_{SA(a)}`.
- Phase 2 — Geometric Side: `ArchBlock`, `RamBlock`, and prime block handling
  - Heat: truncated sum over provided `unram_local` + `PrimeTailBound_heat`.
  - BL: `PrimeBlockBL` includes terms with `2k log p <= X`.
  - RS-positivity mode: `RS_Positivity_LB` gives prime lower bound `[0,0]`.
- Phase 3 — Infrastructure (R* Budget): `RStarBudget` computes `B_{R*}` and verifies `|Ram| <= B_{R*}` via interval containment.
- Phase 4 — Amplified Positivity: checks `Prime - (Arch+Ram)` (evaluate mode) or `0 - UB(Arch+Ram)` (RS mode) is not strictly negative on intervals.
- Phase 5 — ν-Annihilation: computes `sqrt(1-eta)` via `CheckA5Gap` and combines with `epsilon_invariance`; optional tightening from `q >= Q(X,eps)` via `Eps_from_q_X`.
- Phase 6 — Proof Mode: either runs a finite `test_net` (returning `GRH_Verified` if all pass) or returns `AnnOdd_Verified` for the single supplied test.

## API Overview

- `Interval`: closed intervals with outward rounding. Key semantics:
  - Comparisons to zero use the upper endpoint (e.g., `lt0()` means `hi < 0`).
  - Arithmetic is outward-rounded for `+,-,*,/`; `exp/log/sqrt` use conservative float-based padding.
- `TestFunction`: admissible tests; families: `"heat"` (requires `a > 1/2`) and `"BL"`.
- `PrimeLocalTerm`, `RamifiedTerm`: carry interval contributions per `(p^k)` and ramified index `j`.
- `Certificate`: all inputs and knobs (precision, test family, prime mode, proof mode, constants, cutoffs).
- `verify_certificate(cert) -> VerifierReport`: runs Phases 0–6 and returns a result with details.

## Plug Points

- Archimedean: set `cert.arch_value` with a rigorous interval or swap `ArchBlock` to call your gamma-factor routine.
- Ramified: populate `cert.ram_terms` (and `j_max`) or make `RamBlock` call your producer.
- Prime block:
  - Heat: populate `cert.unram_local` and `cert.prime_tail_bound`; optional `p_max/k_max` cutoffs.
  - BL: populate `cert.unram_local` and set `cert.band_limit_X` (cutoff in `2k log p <= X`).
- Annihilation / invariance: set `eta_a5_gap`, `epsilon_invariance`, `alpha_inv`, `kappa_inv`, and optionally `q_modulus`.
- R* constants: set `Ainf`, `beta_m`, `C_Rstar` per your ledger.

## Rigor Notes

- `exp/log/sqrt` currently use float evaluation with symmetric padding for simplicity. For certification use, replace with arb/MPFR-backed implementations or rigorous Decimal routines with directed rounding.
- Basic interval arithmetic uses `decimal` with `ROUND_FLOOR/ROUND_CEILING` to ensure outward rounding.
- Domain checks guard against division by intervals containing 0 and invalid domains for `log`/`sqrt`.

## Repository Layout

- `grh_verifier.py`: main implementation and a tiny synthetic demo in `__main__`.
- `ef_generators.py`: rigorous EF-side generators (Φ spec, Arch/prime/ramified helpers).
- `rigor_quad.py`: validated (and padded-fallback) quadrature utilities.
- `rigor_backend.py`: backend facade for rigorous interval math.
- `zeta_gl1_cert.py`: GL(1) Riemann zeta example (RS mode and BL evaluate mode).
- `grh_cli.py`: CLI runner for modes and parameters.

## Zeta GL(1) quick start

- One-shot demo:
  - `python3 zeta_gl1_cert.py` (RS, BL evaluate, proof_of_GRH, heat evaluate)
- CLI (recommended):
  - `./grh_cli.py --mode max`  # runs proof_of_GRH (RS) + BL + heat
  - `./grh_cli.py --mode all --X 8.0 --a 0.85 --tau 2.5`
  - `./grh_cli.py --mode bl --X 7.0`  (single mode)

Modes
- `rs`: RS-positivity baseline (optionally `--proof` to run a small test net)
- `bl`: BL evaluate mode (generates (p,k) terms with `2k log p <= X`)
- `heat`: Heat evaluate mode (generates (p,k) terms with `2k log p <= X` + Gaussian tail bound)
- `all`: Runs rs + bl + heat
- `max` (default): Runs proof_of_GRH (rs) + bl + heat
