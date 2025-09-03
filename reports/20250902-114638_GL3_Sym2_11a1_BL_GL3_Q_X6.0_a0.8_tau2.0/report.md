# Verification Report — GL3_Sym2_11a1_BL_GL3_Q_X6.0_a0.8_tau2.0
- Result: GRH_Falsified (phase 4)
- Params: {'m': 3, 'K': 'Q', 'Qpi': '11', 't_star': '1.0', 'family': 'BL', 'prime_block_mode': 'evaluate', 'proof_mode': 'ann_odd_only'}
- Details: {'reason': 'Amplified positivity fails', 'GeomAmp': ('-0.7277045912410989654404049137365', '-0.7277045912410989654404049137365'), 'Prime': ('-0.7277045912410989654404049137365', '-0.7277045912410989654404049137365'), 'Arch+Ram': ('0', '0')}

## Visualizations
### Ledger Contributions

![Ledger Contributions](images/contributions.png)

### R* Budget vs Ram

![R* Budget vs Ram](images/rstar_vs_ram.png)

### Annihilation Inequality

![Annihilation Inequality](images/annihilation.png)

### Prime-power Contributions

![Prime-power Contributions](images/prime_terms.png)


## What was verified
- Numerical rigor: Decimal precision set from certificate; interval arithmetic used for comparisons.
- Geometric side: Verified prime block minus (Arch + Ram) is nonnegative (evaluate) or that UB(Arch+Ram) ≤ 0 (RS).
- Infrastructure (R*): Checked |Ram| ≤ B_{R*} using interval bounds.
- Annihilation: Ensured sqrt(1−eta) + epsilon < 1 using upper endpoints; epsilon optionally tightened via q ≥ Q(X,ε).
- Proof-of-GRH mode (if used): Ran the above across a finite test net.
