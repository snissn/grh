
**From RH to GRH to Functoriality — with a geometric mirror on the Heisenberg host**

> **One‑paragraph hook (for an ANTer/L‑functionist).** We rebuild the explicit‑formula (EF) method around a *fixed‑heat* smoothing that lives naturally in the Maaß–Selberg spherical calculus. This lets us **realize the half‑shift as a positive operator**, put **Weil positivity on rails**, and show that a simple **parity involution $\nu$** turns EF positivity into **witnesses** and **annihilators**. On GL$1$ this gives a clean RH ledger with band‑limited tests, a distributional symmetric EF, and a concrete positivity density. The same analytic chassis then **scales to GRH for $L(s,\pi)$** (Parts I–III), produces a **finite, reproducible certificate/verifier**, and powers a **positivity‑driven proof of Langlands functoriality** without trace‑formula comparison. In parallel, a **Heisenberg geometry** recasts the same identities as holonomy/Toeplitz statements on $S^1$, giving a second derivation and a conceptual mirror. Everything below is arranged so that RH is the accessible front door; GRH, Langlands, and Heisenberg are tightly interlocking extensions that both *justify* and *enhance* the RH sequence.

---

## Why RH first?

Experts often meet EF arguments where the normalizations fight back, the archimedean terms won’t sit still, and positivity proofs are brittle. Our RH papers neutralize that by:

1. **Deriving the asymmetric EF automorphically** from the **Eisenstein phase** (clean origin, no folklore),
2. **Packaging the symmetric EF distributionally** on a **band‑limited** class tailored to fixed‑heat smoothing, and
3. **Exhibiting an explicit Weil‑side positive density** (“Mellin–torsion on the modulus line”) that underwrites positivity throughout.

These three pieces make RH the **best entry point**: you can audit constants, see the half‑shift at work in a tame setting, and then watch the same operators lift to GRH and Langlands.

---

## What is genuinely new (and portable)

* **Constructive, positivity‑preserving half‑shift.** Inside the spherical calculus, we build a Loewner‑monotone multiplier that realizes the “$1/2$-shift” *as an honest positive operator* compatible with fixed‑heat smoothing. This single ingredient makes parity, amplifiers, and cone‑separation arguments behave uniformly in conductor.
* **Fixed‑heat PV/Weil bridge.** EF pairings are re‑expressed so that the zero side *is* a **Weil energy** against a positive fixed‑heat operator. The arch block cancels on $\nu$-odd inputs; any off‑line zero produces a **positive‑energy witness**.
* **Finite ledgers and certificates.** Band‑limited tests turn EF into **finite prime ledgers** plus controlled arch/ram budgets; an $\varepsilon$-net yields a **reproducible verifier** for GRH‑type statements.
* **Prime‑side positivity reserve.** A Rankin–Selberg flavored “reserve” gives a quantitative lower bound on the unramified prime block, strong enough to **force functorial transfers** with no comparison of orbital integrals.

**How to sanity‑check in minutes.** (i) On GL$1$, compute the band‑limited dictionary and see that Toeplitz positivity on $S^1$ matches the EF ledger. (ii) Read the half‑shift realization: it’s a plain operator inequality in the spherical model. (iii) Inspect the prime‑side positivity statement at unramified places; the weights are explicit. (iv) For GRH Part II, verify that $\nu$-odd annihilation kills the arch block and makes any off‑line zero light up the Weil energy.

---

## The arc in one narrative (RH → GRH → Langlands; with a geometric mirror)

We begin with **RH** on GL$1$. The EF is rebuilt from the **Eisenstein phase**, then symmetrized on a **band‑limited** class. Fixed‑heat smoothing factors out the sharp cutoffs that usually cause grief; the **half‑shift** becomes a positive operator, so **Weil positivity** is visible and stable. This is not just cleaner—it sets the exact test space we will use everywhere else.

With that spine in place, **GRH** is handled by the **$\nu$-involution**: we split tests into even/odd across the functional equation, show that the archimedean piece vanishes on the odd side, and prove:

* **Annihilation ⇔ GRH** for $L(s,\pi)$: if the odd energy is zero on a fixed‑heat window, we are on the critical line; conversely any off‑line zero creates a witness with strictly positive energy.
* An **amplifier with a strict gap** drives unconditional annihilation of the odd part; the witness mechanism makes the statement falsifiable.

We then **package the method** into a **finite certificate**: for band‑limited tests the EF reduces to finitely many prime blocks plus controlled arch/ram terms; an $\varepsilon$-net over tests yields a verifiable decision procedure.

Finally, for **Langlands**, the same positivity technology gives a **prime‑side reserve**. Together with a **two‑ray** (rank $\ge 2$) or **gap** (rank $1$) argument, this forces the existence of the **functorial transfer $\Pi$** with local compatibility—*no trace‑formula comparison*—closing the circle by linking spectral and arithmetic ledgers through positivity.

In parallel, a six‑part **Heisenberg series** furnishes a **geometric avatar**: residues become holonomy on $S^1$, Szegő kernels provide a model projector, the $\nu$-involution is realized as a simple **reflection × unitary phase** on the modulus line, and the EF ledger appears as a **Toeplitz pairing on the circle**. This mirror derivation both validates the analytic side and quietly seeds forthcoming directions we aren’t discussing here.

---

## Reading plan for a busy professor

1. **Start here (RH, most accessible):** **RH‑I** (asymmetric EF from Eisenstein phase) → **RH‑II (Symmetric EF on band‑limited)** → **RH‑III (Weil positivity density)**.
2. **Then the analytic spine:** **GRH Part I** (infrastructure & half‑shift) → **GRH Part II** ($\nu$, annihilation ⇔ GRH, witness) → **GRH Part III** (finite ledger & certificates).
3. **Applications:** **Langlands — Functoriality via positivity**, then **Geometric Part IV**, then **Reciprocity**.
4. **Geometric mirror (optional but illuminating):** **Heisenberg I–II** (Cauchy/residue as holonomy) → **V–VI** ($\nu$ as geometric operator; EF on $S^1$).

---

# Papers

## RH kernel (the front door)

### RH‑I — *An Automorphic Derivation of the Asymmetric Explicit Formula via the Eisenstein Phase*

**What it does.** Derives the EF **directly from automorphic Eisenstein data** and its phase on the critical line, fixing normalization at source rather than by patchwork.
**Why it matters.** Establishes a canonical EF starting point—what every later identity *means* and *counts* is clear from first principles.
**How to spot‑check.** Verify the phase/normalization by cross‑checking against classical GL$1$ local factors.
**Zenodo:** [10.5281/zenodo.16930060](https://doi.org/10.5281/zenodo.16930060).

### RH‑II — *An Adelic Distributional Framework for the Symmetric Explicit Formula on a Band‑Limited Class*

**What it does.** Packages the **symmetric EF** on a **band‑limited** test class tailored to fixed‑heat smoothing, so spectral zeros and prime/arch pieces pair cleanly.
**Why it matters.** This sets the test space later used for finite ledgers and $\varepsilon$-nets.
**How to spot‑check.** Compute the abelian dictionary $\widehat G(x)=2\sinh(x/2)\,\widehat g(2x)$ and test small band‑limits.
**Zenodo:** [10.5281/zenodo.16930092](https://doi.org/10.5281/zenodo.16930092).

### RH‑III — *Weil Positivity via Mellin–Torsion on the Modulus Line*

**What it does.** Identifies an explicit **Weil‑side positive density** (Mellin–torsion) encoding the fixed‑heat Gram pairing.
**Why it matters.** This positivity table is the keystone for GRH parity arguments and for amplified inequalities in higher rank.
**How to spot‑check.** On GL$1$, the density gives a Toeplitz PSD condition on $S^1$ equivalent to RH.
**Zenodo:** [10.5281/zenodo.16930094](https://doi.org/10.5281/zenodo.16930094).

---

## GRH trilogy (the analytic spine)

### GRH — Part I: *Fixed‑Heat EF Infrastructure for Automorphic $L$-functions on $\mathrm{GL}_m/K$*

**What it does.** Builds the full **fixed‑heat** EF on the Maaß–Selberg ledger: the slice frame and Gram kernel, the **constructive half‑shift** via a Loewner multiplier, conductor‑uniform PV bounds, **prime‑side positivity**, an amplifier, and an **untwisting inequality**.
**Why it matters.** This is the chassis; it makes parity and certificates possible *without* fragile truncations.
**How to spot‑check.** Read the half‑shift realization; then verify unramified prime positivity and the PV growth bound.
**Zenodo:** [10.5281/zenodo.17042903](https://doi.org/10.5281/zenodo.17042903).

### GRH — Part II: *The $\nu$‑Projector, Annihilation vs. GRH, and the Witness Argument*

**What it does.** Proves **$\nu$-odd annihilation ⇔ GRH($\pi$)**, shows any off‑line zero yields a **positive‑energy witness**, and derives **unconditional annihilation** from a strict amplifier gap.
**Why it matters.** Converts positivity into a **falsifiable** equivalence rather than a heuristic.
**How to spot‑check.** Confirm that the archimedean EF block **vanishes on $\nu$-odd inputs** and inspect the witness construction near a hypothetical off‑line zero.
**Zenodo:** [10.5281/zenodo.17042910](https://doi.org/10.5281/zenodo.17042910).

### GRH — Part III: *Applications, Bilinear Annihilation, and GRH Certificates*

**What it does.** Reduces EF to **finite prime ledgers** with certified arch/ram budgets and an **$\varepsilon$-net** of tests, yielding a **certificate/verifier** for GRH decisions.
**Why it matters.** Makes the method **reproducible** and auditable.
**How to spot‑check.** Run the finite‑ledger reduction on a toy band‑limited test and see the budgets close.
**Zenodo:** [10.5281/zenodo.17042970](https://doi.org/10.5281/zenodo.17042970).

---

## Langlands (positivity forcing instead of comparison)

### *Langlands Functoriality via Fixed–Heat Positivity*

**What it does.** Proves the existence of the transfer $r:{}^L G\!\to\!\mathrm{GL}_n$ and **local compatibility** at every place using a **prime‑side positivity reserve**, a **two‑ray deficit** (rank $\ge 2$) or **amplifier gap** (rank $1$), and controlled arch/ram budgets—**no trace‑formula comparison**.
**Why it matters.** Provides a new analytic route to functoriality grounded in EF positivity.
**How to spot‑check.** Verify the stated **reserve** inequality and see how the two‑ray/gap mechanism forces alignment of Satake data.
**Zenodo:** [10.5281/zenodo.17058588](https://doi.org/10.5281/zenodo.17058588).

### *Geometric Langlands via Fixed–Heat Positivity — Part IV*

**What it does.** A stack‑level GL–EF with explicit main/error separation, a **perverse half‑shift** that raises positivity, and a **$\nu$-witness barrier** linking geometric positivity to the analytic suite.
**Why it matters.** Mirrors the analytic ledger in the geometric category and cross‑checks the normalization.
**Zenodo:** [10.5281/zenodo.17058594](https://doi.org/10.5281/zenodo.17058594).

### *Reciprocity via Geometric Rigidity in the Fixed‑Heat Ledger*

**What it does.** Uses geometric rigidity inside the fixed‑heat ledger to formulate **direct reciprocity** compatible with the positivity‑driven transfer.
**Why it matters.** Complements functoriality with a structural reciprocity statement in the same language.
**Zenodo:** [10.5281/zenodo.17058582](https://doi.org/10.5281/zenodo.17058582).

---

## Heisenberg series (the geometric mirror that strengthens the RH→GRH arc)

> An “intriguing side‑quest” that is already useful here and quietly prepares a future direction we are not discussing yet.

### Heisenberg I — *Cauchy, Laurent, and Residues as Holonomy on a Constrained Heisenberg Slice*

**What it does.** Recasts classical complex analysis as **holonomy** on $H_1\simeq S^1\times\mathbb R$; residues are vertical holonomy; Cauchy/Laurent coefficients are Fourier–holonomy pairings.
**Why it matters.** Establishes the kernel $\alpha_H=d\theta$ and the $S^1$ host used later to encode EF.
**Zenodo:** [10.5281/zenodo.17043012](https://doi.org/10.5281/zenodo.17043012).

### Heisenberg II — *Holonomy Universality Beyond the Disk*

**What it does.** Globalizes the picture via conformal pullback: **one Heisenberg slice** hosts Cauchy/residue functionals for arbitrary planar domains (Hardy/Smirnov included).
**Why it matters.** Supplies a robust geometric platform on which the EF can be carried.
**Zenodo:** [10.5281/zenodo.17049478](https://doi.org/10.5281/zenodo.17049478).

### Heisenberg III — *The Heisenberg Shadow: Reframing the Foundations of Analysis*

**What it does.** Connects holonomy to Szegő kernels and semiclassics; exhibits a **host category** in which the circle is initial and rigid.
**Why it matters.** Provides the microlocal backbone for the circle ledger used in RH/GRH.
**Zenodo:** [10.5281/zenodo.17049492](https://doi.org/10.5281/zenodo.17049492).

### Heisenberg IV — *Adelic Modulus, Holonomy, and the Explicit Formula*

**What it does.** Adelizes holonomy and identifies the **geometric side of EF** as total **adelic holonomy** against scale‑profiles.
**Why it matters.** Bridges the geometric calculus to automorphic EF in a way that mirrors the analytic ledger.
**Zenodo:** [10.5281/zenodo.17049498](https://doi.org/10.5281/zenodo.17049498).

### Heisenberg V — *The $\nu$‑Involution as a Geometric Operator on the Adelic Heisenberg Moduli Space*

**What it does.** Realizes the FE‑involution **geometrically** as reflection × unitary phase; couples it to an **amplifier/damping** picture that **annihilates $\nu$-odd holonomy**.
**Why it matters.** Geometrizes the GRH annihilation mechanism and clarifies parity on the circle.
**Zenodo:** [10.5281/zenodo.17049524](https://doi.org/10.5281/zenodo.17049524).

### Heisenberg VI — *Higher Holonomy and the Geometric Annihilation Principle for $\mathrm{GL}(m)$*

**What it does.** Lifts the circle ledger to higher rank; encodes EF as distributions $\mu_Z,\mu_{\text{arith}}$ on $S^1$; proves **Toeplitz PSD ⇔ RH** and a **higher‑rank geometric annihilation** principle.
**Why it matters.** Aligns the geometric and analytic ledgers in general, reinforcing the RH→GRH→Langlands chain.
**Zenodo:** [10.5281/zenodo.17049526](https://doi.org/10.5281/zenodo.17049526).

---

## How the pieces lock together (in one picture)

* **RH** delivers: a canonical EF origin (Eisenstein phase), a **band‑limited** symmetric EF, and an explicit **Weil positivity density**.
* **GRH I–III** lift the machinery: a **half‑shift operator** in the spherical calculus, **parity/annihilation ⇔ GRH**, and **finite ledgers** with **certificates**.
* **Langlands** uses exactly those operators and budgets: a **prime‑side reserve** + **two‑ray/gap forcing** yields **functorial transfer** with **local compatibility**—no trace‑formula comparison.
* **Heisenberg** provides a second, geometric derivation: the **same EF** appears as **Toeplitz** positivity on $S^1$, with $\nu$ realized as a concrete geometric operator.
  Together, they form a single ledger—**fixed‑heat + half‑shift + $\nu$ + holonomy**—that scales from RH to functoriality.

---

## What to verify first (and why it’s quick)

1. **GL$1$ Toeplitz test.** On $S^1$, compute the Toeplitz form against the stated density and match it to the RH symmetric EF coefficients. (Immediate sanity check of normalization/positivity.)
2. **Half‑shift inequality.** Read the Loewner‑order sandwich realizing the half‑shift in the spherical calculus; it’s a clean operator comparison.
3. **Prime reserve.** Inspect the unramified prime block positivity and its quantitative lower bound—this is the forcing engine behind functoriality.
4. **$\nu$-odd annihilation.** Check that the arch block vanishes on odd inputs and that any off‑line zero yields a positive‑energy witness.
5. **Finite ledger.** For a small band‑limit, enumerate the prime blocks and see the budgets close in the certificate scheme.

---

## Closing

This program replaces brittle truncations and ad‑hoc cancellations with a **single, positivity‑preserving calculus** that is **conductor‑uniform** and **geometrically legible**. The RH papers give you the most accessible proof‑of‑concept; the GRH trilogy shows the spine scales; the Langlands papers demonstrate reach; and the Heisenberg series provides an independent, geometric mirror. If you value EF methods that are verifiable, exportable, and built for scrutiny, this ledger is designed to meet you where you are—and to hold up when you push on it.

---
