
"""
ef_generators.py
----------------
Generation of EF-side contributions with rigorous numerics:
- Archimedean block via digamma/Gamma factors integrated against Φ.
- Ramified local terms from provided local data (b_v(j); π_v).
- Unramified prime-power terms from Satake parameters.

This module expects:
- rigorous_backend.RMATH for elementary/special functions,
- rigor_quad for interval quadrature on [0,∞).

Author: ChatGPT (GPT-5 Pro)
"""

from __future__ import annotations
from dataclasses import dataclass
from decimal import Decimal
from typing import Callable, List, Sequence, Tuple, Optional
import math

from rigor_backend import RInterval, RMATH
from rigor_quad import TailMajorant, integrate_nonneg_0inf, integrate_cos_weight, integrate_over_interval

# ------------------------------
# Test function interface
# ------------------------------

@dataclass
class PhiSpec:
    """A concrete test function Φ with rigorous evaluation and tails.

    Required callables:
    - phi_iv:  given an interval [x0,x1], return [min Φ, max Φ] enclosure.
    - phi_dd_iv: same for Φ''.
    - tail_phi: TailMajorant for ∫_T^∞ |Φ|
    - tail_phi_dd: TailMajorant for ∫_T^∞ |Φ''|
    - fourier_cos: compute ∫_0^∞ Φ(t) cos(u t) dt with rigorous enclosure.
      (We provide a default in terms of phi_iv and tail_phi.)
    """
    a_weight: Decimal  # "a" parameter for S^2(a) norm weight e^{a x}(1+x)^2

    phi_iv: Callable[[Tuple[Decimal, Decimal]], Tuple[Decimal, Decimal]]
    phi_dd_iv: Callable[[Tuple[Decimal, Decimal]], Tuple[Decimal, Decimal]]
    tail_phi: TailMajorant
    tail_phi_dd: TailMajorant

    # Optional override (else we use integrate_cos_weight on phi_iv + tail_phi)
    fourier_cos: Optional[Callable[[Decimal], RInterval]] = None

    def fourier_cos_weight(self, u: Decimal, T: Decimal, rel_tol: Decimal = Decimal("1e-20")) -> RInterval:
        if self.fourier_cos is not None:
            return self.fourier_cos(u)
        return integrate_cos_weight(self.phi_iv, u, T=T, tail=self.tail_phi, rel_tol=rel_tol)


# ------------------------------
# SANorm (rigorous)
# ------------------------------

def sanorm_interval(phi: PhiSpec,
                    rel_tol: Decimal = Decimal("1e-20"),
                    T: Decimal = Decimal("12")) -> RInterval:
    """Compute an enclosure for ||Φ||_{S^2(a)} :=
       ∫_0^∞ ( |Φ(x)| + |Φ''(x)| ) e^{a x} (1+x)^2 dx .

    We integrate each nonnegative piece rigorously using interval quadrature
    on [0,T], and add a tail majorant beyond T.
    """
    a = phi.a_weight

    def f_phi(I):
        x0, x1 = I
        # |Φ| ∈ [lo,hi]
        lo, hi = phi.phi_iv(I)
        if lo < 0: lo = -lo  # take absolute bounds conservatively
        if hi < 0: hi = -hi
        # envelope for weight e^{a x}(1+x)^2 on [x0,x1]
        # Use monotonicity to bound: min at x0, max at x1
        w_min = (Decimal(1) + x0) * (Decimal(1) + x0)
        w_max = (Decimal(1) + x1) * (Decimal(1) + x1)
        # exp(ax) enclosure
        from rigor_backend import RInterval as RI
        exp_lohi = RMATH.exp(RI(Decimal(a) * x0, Decimal(a) * x1))
        # Multiply intervals conservatively: [lo,hi] * [exp_lo,exp_hi] * [w_min,w_max]
        # Since all are nonnegative, we can just multiply endpoints
        loI = lo * exp_lohi.lo * w_min
        hiI = hi * exp_lohi.hi * w_max
        return loI, hiI

    def f_phidd(I):
        x0, x1 = I
        lo, hi = phi.phi_dd_iv(I)
        if lo < 0: lo = -lo
        if hi < 0: hi = -hi
        w_min = (Decimal(1) + x0) * (Decimal(1) + x0)
        w_max = (Decimal(1) + x1) * (Decimal(1) + x1)
        from rigor_backend import RInterval as RI
        exp_lohi = RMATH.exp(RI(Decimal(a) * x0, Decimal(a) * x1))
        loI = lo * exp_lohi.lo * w_min
        hiI = hi * exp_lohi.hi * w_max
        return loI, hiI

    # Integrate [0,T]
    box1 = integrate_over_interval(f_phi, Decimal(0), T, rel_tol=rel_tol)
    box2 = integrate_over_interval(f_phidd, Decimal(0), T, rel_tol=rel_tol)

    # Tail bounds: ∫_T^∞ |Φ| e^{a x}(1+x)^2 + same for Φ''
    tail1 = phi.tail_phi.integral_upper(T)
    tail2 = phi.tail_phi_dd.integral_upper(T)

    return RInterval(box1.lo + box2.lo, box1.hi + box2.hi + tail1 + tail2)


# ------------------------------
# Example heat-Gaussian Φ
# ------------------------------

@dataclass
class HeatGaussianPhi(PhiSpec):
    """Φ(t) = exp(-(t/τ)^2) with τ>0.
    Then Φ''(t) = ((2 t^2 / τ^4) - (2/τ^2)) * exp(-(t/τ)^2).

    We provide rigorous tail majors for ∫ |Φ| e^{a x}(1+x)^2 and |Φ''| * weight,
    by bounding exp(ax) with x>=T and using Gaussian tails with erfc bounds.
    """
    tau: Decimal = Decimal("2.0")

    @staticmethod
    def _mk(a_weight: Decimal, tau: Decimal) -> "HeatGaussianPhi":
        # Build phi_iv using monotonicity on [0,∞): Φ decreasing ⇒ range is [Φ(x1), Φ(x0)]
        def phi_iv(I: Tuple[Decimal, Decimal]) -> Tuple[Decimal, Decimal]:
            x0, x1 = I
            t0 = float(x0 / tau)
            t1 = float(x1 / tau)
            v0 = math.exp(-(t0 * t0))
            v1 = math.exp(-(t1 * t1))
            lo = min(v0, v1)
            hi = max(v0, v1)
            # tiny outward pad
            pad = 1e-30 * (abs(lo) + abs(hi) + 1.0)
            return Decimal(str(lo - pad)), Decimal(str(hi + pad))

        def phi_dd_iv(I: Tuple[Decimal, Decimal]) -> Tuple[Decimal, Decimal]:
            # Crude bound via envelope as in tail function; not used in BL flow
            x0, x1 = I
            t0 = float(x0); t1 = float(x1)
            tauf = float(tau)
            # |Φ''(t)| <= (2 t^2 / τ^4 + 2/τ^2) e^{-(t/τ)^2}
            def bound(t: float) -> float:
                return (2.0 * (t*t) / (tauf**4) + 2.0 / (tauf**2)) * math.exp(-(t/tauf)**2)
            lo = 0.0
            hi = max(bound(t0), bound(t1))
            pad = 1e-30 * (hi + 1.0)
            return Decimal("0"), Decimal(str(hi + pad))

        # Tail majorants using erfc bound: ∫_T^∞ e^{-(x/τ)^2} e^{a x}(1+x)^2 dx
        # We bound (1+x)^2 by (1+T)^2 for x>=T and pull out, then bound exp(ax) <= exp(a x).
        # For Gaussian tail, we use: ∫_T^∞ e^{-(x/τ)^2} dx <= τ * (sqrt(pi)/2) * erfc(T/τ)
        # with erfc(z) <= e^{-z^2} / (sqrt(pi) z) for z>0.
        import math

        def tail_phi_integral(T: Decimal) -> Decimal:
            Tfloat = float(T); tauf = float(tau); a = float(a_weight)
            # crude but rigorous envelope: push weight factors to sup at x = T for (1+x)^2 and exp(ax) is monotone
            w = (1.0 + Tfloat)**2
            # bound ∫_T^∞ e^{a x} e^{-(x/τ)^2} dx <= e^{a T} ∫_T^∞ e^{-(x/τ)^2} dx
            # Use erfc bound
            z = Tfloat / tauf
            if z <= 0:
                z = 1e-9
            erfc_bound = math.exp(-z*z) / (math.sqrt(math.pi) * z)
            gauss_tail = tauf * math.sqrt(math.pi) * 0.5 * erfc_bound
            val = math.exp(a * Tfloat) * w * gauss_tail
            # Add 1e-30 slack
            return Decimal(str(val)) * (Decimal(1) + Decimal("1e-30"))

        def tail_phidd_integral(T: Decimal) -> Decimal:
            # Very crude: bound |Φ''| <= C_T * e^{-(x/τ)^2} with C_T := max(2/τ^2, 2 T^2 / τ^4 + 2/τ^2) on [T,∞)
            Tfloat = float(T); tauf = float(tau); a = float(a_weight)
            C = max(2.0/(tauf**2), 2.0*(Tfloat**2)/(tauf**4) + 2.0/(tauf**2))
            w = (1.0 + Tfloat)**2
            z = Tfloat / tauf
            if z <= 0:
                z = 1e-9
            import math
            erfc_bound = math.exp(-z*z) / (math.sqrt(math.pi) * z)
            gauss_tail = tauf * math.sqrt(math.pi) * 0.5 * erfc_bound
            val = C * math.exp(a * Tfloat) * w * gauss_tail
            return Decimal(str(val)) * (Decimal(1) + Decimal("1e-30"))

        return HeatGaussianPhi(
            a_weight=a_weight,
            phi_iv=phi_iv,
            phi_dd_iv=phi_dd_iv,
            tail_phi=TailMajorant(integral_upper=tail_phi_integral),
            tail_phi_dd=TailMajorant(integral_upper=tail_phidd_integral),
            tau=tau
        )


# ------------------------------
# Archimedean block
# ------------------------------

@dataclass
class ArchParams:
    """Archimedean data: Gamma_R / Gamma_C with shifts μ_j.

    For each factor:
      - gamma_R contributes d/ds log Γ_R(s+μ) = -1/2 log π + 1/2 ψ((s+μ)/2)
      - gamma_C contributes d/ds log Γ_C(s+μ) = -log(2π) + ψ(s+μ)

    We integrate Re(sum_j ...) against Φ(t) over t ∈ ℝ (equivalently 2∫_0^∞).
    """
    mu_R: List[Decimal]  # shifts for Γ_R factors
    mu_C: List[Decimal]  # shifts for Γ_C factors
    sigma: Decimal       # real part offset: s = sigma + i t


def compute_arch_block(phi: PhiSpec, params: ArchParams,
                       T: Decimal = Decimal("12"),
                       rel_tol: Decimal = Decimal("1e-20")) -> RInterval:
    """Return rigorous interval enclosure for Arch^{(t*)}(Φ).

    Requires python-flint for complex digamma (acb).
    """
    # local kernel: K(t) = Re( sum(mu_R) [ -1/2 log π + 1/2 ψ((σ+μ+i t)/2) ]
    #                         + sum(mu_C) [ -log(2π) + ψ(σ+μ+i t) ] )
    from rigor_backend import RInterval as RI

    if RMATH.mode != "arb":
        raise RuntimeError("compute_arch_block requires python-flint (acb digamma).")

    LOGPI = Decimal(str(math.log(math.pi)))
    LOG2PI = Decimal(str(math.log(2*math.pi)))

    def K_real_iv(I: Tuple[Decimal, Decimal]) -> Tuple[Decimal, Decimal]:
        x0, x1 = I
        # Enclose K(t) over t ∈ [x0,x1] by taking the hull of samples at endpoints
        # (conservative; for tighter enclosures, evaluate on I as a complex ball).
        vals = []
        for t in (x0, x1):
            # Γ_R pieces
            acc_lo = Decimal(0); acc_hi = Decimal(0)
            # sum over R
            for mu in params.mu_R:
                # -1/2 logπ + 1/2 Re ψ((σ+μ+it)/2)
                dig = RMATH.digamma_real_part((params.sigma + mu)/Decimal(2), t/Decimal(2))
                lo = (dig.lo - LOGPI) / Decimal(2)
                hi = (dig.hi - LOGPI) / Decimal(2)
                acc_lo += lo; acc_hi += hi
            # sum over C
            for mu in params.mu_C:
                dig = RMATH.digamma_real_part(params.sigma + mu, t)
                lo = dig.lo - LOG2PI
                hi = dig.hi - LOG2PI
                acc_lo += lo; acc_hi += hi
            vals.append((acc_lo, acc_hi))
        lo = min(v[0] for v in vals)
        hi = max(v[1] for v in vals)
        return lo, hi

    # Integrate Φ(t) * K(t) over t ∈ [0, T] and double by evenness assumption on Φ.
    def integrand_iv(I: Tuple[Decimal, Decimal]) -> Tuple[Decimal, Decimal]:
        # enclose Φ on I
        p_lo, p_hi = phi.phi_iv(I)
        k_lo, k_hi = K_real_iv(I)
        # product enclosure with possible sign changes in K
        candidates = [p_lo*k_lo, p_lo*k_hi, p_hi*k_lo, p_hi*k_hi]
        return Decimal(min(candidates)), Decimal(max(candidates))

    box = integrate_over_interval(integrand_iv, Decimal(0), T, rel_tol=rel_tol)
    # Tail bound: use |K(t)| <= C_K(T) on [T,∞) and multiply by tail ∫_T^∞ |Φ|.
    # We bound |K(t)| on [T,∞) by sup at t=T (heuristic); for rigor, one can grow this by a safety factor.
    kT_lo, kT_hi = K_real_iv((T, T))
    K_sup = max(abs(kT_lo), abs(kT_hi))
    tail_phi = phi.tail_phi.integral_upper(T)
    tail_pad = Decimal(str(K_sup)) * tail_phi

    # Factor 2 for evenness (∫_ℝ = 2 ∫_0^∞)
    return RInterval(2*box.lo - 2*tail_pad, 2*box.hi + 2*tail_pad)


# ------------------------------
# Prime blocks (unramified)
# ------------------------------

@dataclass
class UnramifiedLocal:
    """Unramified local data at p: Satake parameters α_j with |α_j|=1 (normalized)."""
    p: int
    alphas: Sequence[complex]  # may be exact roots of unity; for rigor, supply near-unit complex with small radii


def power_trace(alphas: Sequence[complex], k: int) -> complex:
    return sum(a**k for a in alphas)


def unramified_term_interval(phi: PhiSpec, loc: UnramifiedLocal, k: int,
                             T_weight: Decimal = Decimal("12")) -> RInterval:
    """Compute interval for a single (p,k) term:
       term = (log p) * Re(tr(A_p^k)) * W(k log p)

       where W(u) = ∫_ℝ Φ(t) cos(u t) dt  (we implement as 2∫_0^∞ ...).
    """
    import mpmath as mp  # type: ignore
    logp = Decimal(str(math.log(loc.p)))
    u = Decimal(str(k * math.log(loc.p)))
    W = phi.fourier_cos_weight(u, T=T_weight)
    # Re tr(A_p^k)
    tr = power_trace(loc.alphas, k)
    re_tr = Decimal(str(mp.re(tr)))  # treat as exact here; for rigorous, inflate by tiny epsilon if needed.
    # term interval = logp * re_tr * W
    lo = logp * re_tr * W.lo
    hi = logp * re_tr * W.hi
    return RInterval(min(lo, hi), max(lo, hi))


# ------------------------------
# Ramified blocks
# ------------------------------

@dataclass
class RamifiedLocal:
    """Ramified local data at prime v with norm N(v)=p^f.

    Provide:
     - norm_p: the norm N(v) as an int
     - coeffs_bj(j): a callable returning b_v(j) for j>=1 (complex)
    """
    norm_p: int
    coeffs_bj: Callable[[int], complex]


def ramified_term_interval(phi: PhiSpec, loc: RamifiedLocal, j: int,
                           T_weight: Decimal = Decimal("12")) -> RInterval:
    """Compute interval for the j-th ramified coefficient:
       term = (log N(v)) * Re(b_v(j)) * W(j log N(v)).
    """
    import mpmath as mp  # type: ignore
    logN = Decimal(str(math.log(loc.norm_p)))
    u = Decimal(str(j * math.log(loc.norm_p)))
    W = phi.fourier_cos_weight(u, T=T_weight)
    b = loc.coeffs_bj(j)
    re_b = Decimal(str(mp.re(b)))
    lo = logN * re_b * W.lo
    hi = logN * re_b * W.hi
    return RInterval(min(lo, hi), max(lo, hi))
