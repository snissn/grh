# Adelic Geometric Pairing \(Q(g)\) — One‑Page Report

**Goal.** Evaluate the adelic geometric pairing
\[
Q(g)=\sum_v\langle D_v,W\rangle=\sum_p\langle D_p,W\rangle+\langle D_\infty,W\rangle
\]
for **band‑limited** tests \(g\) and verify the identity \(Q(g)=P(g)-A(g)\) with machine‑precision checks on the discrete part and tight tolerances at the Archimedean place.

---

## Normalization (Maaß–Selberg) and dictionary

- **Test class.** \(g\) is real, even, with cosine profile \(g_b\in C_c^\infty(\mathbb R)\) and
  \(g(\tau)=\dfrac{1}{\pi}\int_0^\infty g_b(x)\cos(x\tau)\,dx.\) :contentReference[oaicite:3]{index=3}
- **Weil test on the modulus line.** \(W(u)=\dfrac{1}{\pi}(1-e^{-u})\,g_b(2u)\) for \(u\ge 0\).  
  Finite places: \(D_p(u)=\sum_{k\ge1}(\log p)\,\delta(u-k\log p)\).  
  Infinite place: \(K_\infty(u)=-\dfrac{e^{-u}}{1-e^{-2u}}=-\dfrac{1}{2\sinh u}\). :contentReference[oaicite:4]{index=4}
- **Arithmetic block.** 
  \[
  P(g)=\frac{1}{\pi}\sum_{p}\sum_{k\ge1}(\log p)\bigl(1-p^{-k}\bigr)\,g_b\!\left(2k\log p\right).
  \]
  Band‑limiting makes this a finite sum. :contentReference[oaicite:5]{index=5}
- **Archimedean block.**
  \[
  A(g)=\frac{1}{2\pi}\int_{\mathbb R} g(\tau)\,\Re\{\psi(i\tau)-\psi(\tfrac12+i\tau)\}\,d\tau,\qquad
  \langle D_\infty,W\rangle=-A(g).
  \]
  The equality is verified numerically and matches the Eisenstein‑phase normalization. :contentReference[oaicite:6]{index=6} :contentReference[oaicite:7]{index=7}
- **Useful \(u\)-side cross‑check.** 
  \[
  A(g)=\frac{1}{2\pi}\int_0^X \frac{g_b(t)}{e^{t/2}+1}\,dt
  \]
  (with \(X=\) support radius), used as a high‑accuracy reference for the \(\tau\)-integral. :contentReference[oaicite:8]{index=8}

---

## What the notebook computes

1. **Discrete (finite) block** two ways:  
   (i) sampling \(W(k\log p)\) with weights \(\log p\);  
   (ii) the closed form \(\tfrac{1}{\pi}(\log p)(1-p^{-k})g_b(2k\log p)\).  
   The difference tests machine precision.

2. **Archimedean block** two ways:  
   (i) \(u\)-side quadrature \(\int_0^{X/2}K_\infty(u)W(u)\,du\) (near‑zero handled by \(u=t^2\) and a small‑\(u\) series for \(1/(2\sinh u)\));  
   (ii) \(\tau\)-side Simpson integral of \(A(g)\) with \(g(\tau)\) from cosine inversion on \([0,X]\). The integrand is even; we integrate over \([\tau_0,T]\) and double to avoid the digamma singularity at \(0\).

3. **Main identity** \(Q(g)=P(g)-A(g)\) using the \(u\)-side \(A(g)\) as reference, and reporting the \(\tau\)-side value as an independent check. The normalization aligns with the Eisenstein‑phase derivation of the explicit formula and its geometric/spectral blocks. :contentReference[oaicite:9]{index=9}

---

## Acceptance checks (and tolerances)

- **Silent support (no prime spikes):** if \(\operatorname{supp} g_b \subset (0,2\log 2)\), then \(P(g)=0\) and  
  \(Q(g)=\langle D_\infty,W\rangle=-A(g)\) numerically.
- **Discrete block exactness:** \(|P_{\text{closed}}-P_{\text{sample}}|<10^{-12}\) (machine precision).
- **Archimedean agreement:** \(|\langle D_\infty,W\rangle + A_u|<10^{-12}\) and \(|\langle D_\infty,W\rangle + A_\tau|<10^{-6}\).
- **Main identity:** \(|Q-(P-A_u)|<10^{-12}\).

The notebook executes all tests automatically and raises an error if any tolerance is violated.

---

## How to run

Open `adelic_pairing_verification.ipynb` and run all cells. The final cell prints **“All acceptance tests PASSED.”** and displays a summary table with:
\[
P_{\text{closed}},\,P_{\text{sample}},\,\langle D_\infty,W\rangle,\,A_u,\,A_\tau,\,Q,\,P-A_u
\]
and the small residuals used for validation.

---

### References / provenance

- **Eisenstein‑phase route & Maaß–Selberg normalization** (derivation of \(A(g),P(g),J(g)\) and the asymmetric explicit formula). :contentReference[oaicite:10]{index=10}  
- **Adelic distributional framework on the modulus line** (Weil test \(W\), local distributions \(D_v\), and the identity \(\langle D_\infty,W\rangle=-A(g)\)). :contentReference[oaicite:11]{index=11}  
- **Mellin–torsion / \(u\)-side cross‑checks and positivity lens** (closed \(u\)-side form for \(A(g)\)). :contentReference[oaicite:12]{index=12}

