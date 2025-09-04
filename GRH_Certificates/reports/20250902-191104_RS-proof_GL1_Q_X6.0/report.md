# Verification Report — RS-proof_GL1_Q_X6.0
- Result: GRH_Verified (phase 6)
- Params: {'m': 1, 'K': 'Q', 'Qpi': '1', 't_star': '1.0', 'family': 'BL', 'prime_block_mode': 'rs_lower_bound', 'proof_mode': 'proof_of_GRH'}
- Details: {'reason': 'All tests in net passed'}

## Visualizations
### Ledger Contributions

![Ledger Contributions](images/contributions.png)

### R* Budget vs Ram

![R* Budget vs Ram](images/rstar_vs_ram.png)

### Annihilation Inequality

![Annihilation Inequality](images/annihilation.png)


## What was verified
- Numerical rigor: Decimal precision set from certificate; interval arithmetic used for comparisons.
- Geometric side: Verified prime block minus (Arch + Ram) is nonnegative (evaluate) or that UB(Arch+Ram) ≤ 0 (RS).
- Infrastructure (R*): Checked |Ram| ≤ B_{R*} using interval bounds.
- Annihilation: Ensured sqrt(1−eta) + epsilon < 1 using upper endpoints; epsilon optionally tightened via q ≥ Q(X,ε).
- Proof-of-GRH mode (if used): Ran the above across a finite test net.
