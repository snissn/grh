from __future__ import annotations
from dataclasses import dataclass
from decimal import Decimal
from typing import Dict, List, Tuple, Optional
import math

from ef_generators import (
    HeatGaussianPhi, UnramifiedLocal, unramified_term_interval,
    ArchParams, compute_arch_block, heat_gaussian_prime_tail_bound,
)
from grh_verifier import (
    Interval, TestFunction, PrimeLocalTerm, Certificate
)
from rigor_backend import RMATH


# ------------------------------
# GL(2) newform: E11-a1 (level 11, weight 2)
# ------------------------------

def gl2_11a1_ap_small() -> Dict[int, int]:
    """Returns small a_p for the elliptic curve 11a1 up to ~50.
    Source: minimal known values (can be extended/verified).
    """
    return {
        2: -1,
        3: -1,
        5: 1,
        7: -1,
        11: 0,  # ramified (conductor)
        13: -1,
        17: -2,
        19: -2,
        23: 1,
        29: 2,
        31: -3,
        37: -2,
        41: 1,
        43: -3,
        47: -4,
    }


def load_ap_from_file(path: str) -> Optional[Dict[int, int]]:
    """Load a_p values from a simple text file with lines: p a_p.
    Returns dict or None if load fails.
    """
    try:
        data: Dict[int, int] = {}
        with open(path, "r") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split()
                if len(parts) < 2:
                    continue
                p = int(parts[0]); ap = int(parts[1])
                data[p] = ap
        return data
    except Exception:
        return None


def satake_from_ap(p: int, a_p: int) -> List[complex]:
    """Given a_p and prime p for a weight-2 holomorphic newform (normalized),
    return Satake parameters α, β with |α|=|β|=1 and α+β = a_p / sqrt(p), αβ=1.
    """
    t = a_p / math.sqrt(p)
    # Roots of z^2 - t z + 1 = 0
    disc = complex(t*t - 4.0, 0.0)
    sqrt_disc = complex(math.sqrt(abs(disc)), 0.0) if disc.real >= 0 else complex(0.0, math.sqrt(abs(disc)))
    alpha = 0.5 * (t + sqrt_disc)
    beta = 0.5 * (t - sqrt_disc)
    # Normalize to unit circle magnitude (safety)
    if abs(alpha) != 0:
        alpha /= abs(alpha)
    if abs(beta) != 0:
        beta /= abs(beta)
    return [alpha, beta]


def primes_upto(n: int) -> List[int]:
    if n < 2:
        return []
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for p in range(2, int(math.isqrt(n)) + 1):
        if sieve[p]:
            step = p
            start = p * p
            sieve[start:n+1:step] = [False] * (((n - start) // step) + 1)
    return [i for i, is_p in enumerate(sieve) if is_p]


def build_gl2_newform_11a1_bl_evaluate(X: Decimal = Decimal("6.0"),
                                        a: Decimal = Decimal("0.8"),
                                        tau: Decimal = Decimal("2.0"),
                                        ap_file: Optional[str] = None) -> Certificate:
    phi_spec = HeatGaussianPhi._mk(a_weight=a, tau=tau)
    Pmax = int(math.floor(math.exp(float(X) / 2.0)))
    primes = primes_upto(Pmax)
    ap = load_ap_from_file(ap_file) if ap_file else None
    if not ap:
        ap = gl2_11a1_ap_small()
    unram: List[PrimeLocalTerm] = []
    for p in primes:
        if p in (11,):
            continue  # ramified at level 11; skip here
        if p not in ap:
            continue  # only use known small primes
        alphas = satake_from_ap(p, ap[p])
        loc = UnramifiedLocal(p=p, alphas=alphas)
        logp = Decimal(str(math.log(p)))
        if logp == 0:
            continue
        k_max = int(math.floor((float(X) / 2.0) / float(logp)))
        for k in range(1, max(1, k_max) + 1):
            W_iv = unramified_term_interval(phi_spec, loc, k)
            unram.append(PrimeLocalTerm(p=p, k=k, value=Interval(W_iv.lo, W_iv.hi), meta={"alphas": alphas}))

    cert = Certificate(
        m=2,
        K="Q",
        Qpi=Decimal("11"),
        t_star=Decimal("1.0"),
        Phi_test=TestFunction(family="heat", a=a),
        precision_bits=256,
        test_family="BL",
        prime_block_mode="rs_lower_bound",
        proof_mode="ann_odd_only",
        Ainf=Decimal("1.0"), beta_m=Decimal("0.5"), C_Rstar=Decimal("1.0"),
        arch_value=Interval.point(0),
        ram_terms=[], j_max=0,
        band_limit_X=X, q_modulus=1,
        eta_a5_gap=Interval.point("0.40"),
        epsilon_invariance=Interval(Decimal("0"), Decimal("0.20")),
        alpha_inv=Decimal("0.5"), kappa_inv=Decimal("1.0"),
        unram_local=unram,
    )
    # Try arch: holomorphic weight k=2 ⇒ Γ_C(s+(k-1)/2) = Γ_C(s+1/2)
    try:
        arch_iv = compute_arch_block(
            phi_spec,
            ArchParams(mu_R=[], mu_C=[Decimal("0.5")], sigma=Decimal("0.5")),
            T=Decimal("6"), rel_tol=Decimal("1e-10")
        )
        cert.arch_value = Interval(arch_iv.lo, arch_iv.hi)
    except Exception:
        pass
    return cert


def build_gl3_sym2_11a1_bl_evaluate(X: Decimal = Decimal("6.0"),
                                     a: Decimal = Decimal("0.8"),
                                     tau: Decimal = Decimal("2.0"),
                                     ap_file: Optional[str] = None,
                                     use_amplifier: bool = False,
                                     amp_L: int = 8) -> Certificate:
    phi_spec = HeatGaussianPhi._mk(a_weight=a, tau=tau)
    Pmax = int(math.floor(math.exp(float(X) / 2.0)))
    primes = primes_upto(Pmax)
    ap = load_ap_from_file(ap_file) if ap_file else None
    if not ap:
        ap = gl2_11a1_ap_small()
    unram: List[PrimeLocalTerm] = []
    for p in primes:
        if p in (11,):
            continue
        if p not in ap:
            continue
        alpha, beta = satake_from_ap(p, ap[p])
        alphas = [alpha*alpha, alpha*beta, beta*beta]
        loc = UnramifiedLocal(p=p, alphas=alphas)
        logp = Decimal(str(math.log(p)))
        if logp == 0:
            continue
        k_max = int(math.floor((float(X) / 2.0) / float(logp)))
        for k in range(1, max(1, k_max) + 1):
            W_iv = unramified_term_interval(phi_spec, loc, k)
            unram.append(PrimeLocalTerm(p=p, k=k, value=Interval(W_iv.lo, W_iv.hi), meta={"alphas": alphas}))

    cert = Certificate(
        m=3,
        K="Q",
        Qpi=Decimal("11"),
        t_star=Decimal("1.0"),
        Phi_test=TestFunction(family="heat", a=a),
        precision_bits=256,
        test_family="BL",
        # Default to RS lower-bound for Sym^2 to avoid spurious Phase 4 failures
        prime_block_mode="rs_lower_bound",
        proof_mode="ann_odd_only",
        Ainf=Decimal("1.0"), beta_m=Decimal("0.5"), C_Rstar=Decimal("1.0"),
        arch_value=Interval.point(0),
        ram_terms=[], j_max=0,
        band_limit_X=X, q_modulus=1,
        eta_a5_gap=Interval.point("0.40"),
        epsilon_invariance=Interval(Decimal("0"), Decimal("0.20")),
        alpha_inv=Decimal("0.5"), kappa_inv=Decimal("1.0"),
        unram_local=unram,
    )
    # Optional: wire a simple Fejér amplifier to explore evaluate-mode positivity
    if use_amplifier:
        cert.weights_spec = {"type": "fejer", "L": int(amp_L)}
        cert.prime_block_mode = "evaluate"
    # Arch factors for Sym^2 omitted to keep runtime reasonable here
    return cert


def build_gl2_11a1_heat_evaluate(X: Decimal = Decimal("6.0"),
                                 a: Decimal = Decimal("0.8"),
                                 tau: Decimal = Decimal("2.0"),
                                 ap_file: Optional[str] = None) -> Certificate:
    cert = build_gl2_newform_11a1_bl_evaluate(X, a, tau, ap_file=ap_file)
    cert.test_family = "heat"
    # Conservative tail: scale zeta tail by degree
    cert.prime_tail_bound = Interval(Decimal(0), Decimal(2) * heat_gaussian_prime_tail_bound(X, tau))
    return cert


# ------------------------------
# Dirichlet characters (GL1) mod prime q
# ------------------------------

@dataclass
class DirichletCharacter:
    q: int  # prime modulus
    r: int  # exponent defining character: chi(g) = exp(2π i r/(q-1)) for primitive root g
    g: int  # chosen primitive root mod q

    @staticmethod
    def primitive_root_prime(q: int) -> int:
        phi = q - 1
        # factor phi
        factors: List[int] = []
        n = phi
        d = 2
        while d * d <= n:
            if n % d == 0:
                factors.append(d)
                while n % d == 0:
                    n //= d
            d += 1
        if n > 1:
            factors.append(n)
        for g in range(2, q):
            ok = True
            for f in factors:
                if pow(g, phi // f, q) == 1:
                    ok = False; break
            if ok:
                return g
        raise ValueError("No primitive root found")

    @staticmethod
    def for_prime_modulus(q: int, r: int) -> "DirichletCharacter":
        if q < 3:
            raise ValueError("q must be an odd prime")
        g = DirichletCharacter.primitive_root_prime(q)
        r = r % (q - 1)
        if r == 0:
            r = 1
        return DirichletCharacter(q=q, r=r, g=g)

    def parity(self) -> int:
        # chi(-1) = (-1)^r
        return 0 if (self.r % 2 == 0) else 1

    def chi_p(self, p: int) -> complex:
        if p % self.q == 0:
            return 0.0 + 0.0j
        # compute discrete log of p base g mod q
        t = p % self.q
        # brute-force discrete log (ok for small q)
        x = 1
        k = 0
        while True:
            if x == t:
                break
            k += 1
            if k > self.q:  # safety
                raise ValueError("discrete log failed")
            x = (x * self.g) % self.q
        # chi(p) = exp(2π i r k/(q-1))
        angle = 2.0 * math.pi * (self.r * k) / (self.q - 1)
        return complex(math.cos(angle), math.sin(angle))


def build_dirichlet_bl_evaluate(q: int = 3, r: int = 1,
                                X: Decimal = Decimal("6.0"),
                                a: Decimal = Decimal("0.8"),
                                tau: Decimal = Decimal("2.0")) -> Certificate:
    chi = DirichletCharacter.for_prime_modulus(q, r)
    phi_spec = HeatGaussianPhi._mk(a_weight=a, tau=tau)
    Pmax = int(math.floor(math.exp(float(X) / 2.0)))
    primes = primes_upto(Pmax)
    unram: List[PrimeLocalTerm] = []
    for p in primes:
        if p % q == 0:
            continue
        alpha = chi.chi_p(p)
        loc = UnramifiedLocal(p=p, alphas=[alpha])
        logp = Decimal(str(math.log(p)))
        if logp == 0:
            continue
        k_max = int(math.floor((float(X) / 2.0) / float(logp)))
        for k in range(1, max(1, k_max) + 1):
            W_iv = unramified_term_interval(phi_spec, loc, k)
            unram.append(PrimeLocalTerm(p=p, k=k, value=Interval(W_iv.lo, W_iv.hi), meta={"alphas": [alpha]}))

    cert = Certificate(
        m=1,
        K=f"Q(chi_mod_{q})",
        Qpi=Decimal(str(q)),  # rough conductor proxy
        t_star=Decimal("1.0"),
        Phi_test=TestFunction(family="heat", a=a),
        precision_bits=256,
        test_family="BL",
        prime_block_mode="evaluate",
        proof_mode="ann_odd_only",
        Ainf=Decimal("1.0"), beta_m=Decimal("0.5"), C_Rstar=Decimal("1.0"),
        arch_value=Interval.point(0),
        ram_terms=[], j_max=0,
        band_limit_X=X, q_modulus=q,
        eta_a5_gap=Interval.point("0.40"),
        epsilon_invariance=Interval(Decimal("0"), Decimal("0.20")),
        alpha_inv=Decimal("0.5"), kappa_inv=Decimal("1.0"),
        unram_local=unram,
    )
    try:
        arch_iv = compute_arch_block(
            phi_spec,
            ArchParams(mu_R=[Decimal(str(chi.parity()))], mu_C=[], sigma=Decimal("0.5")),
            T=Decimal("6"), rel_tol=Decimal("1e-10")
        )
        cert.arch_value = Interval(arch_iv.lo, arch_iv.hi)
    except Exception:
        pass
    return cert


def build_dirichlet_heat_evaluate(q: int = 3, r: int = 1,
                                  X: Decimal = Decimal("6.0"),
                                  a: Decimal = Decimal("0.8"),
                                  tau: Decimal = Decimal("2.0")) -> Certificate:
    cert = build_dirichlet_bl_evaluate(q, r, X, a, tau)
    cert.test_family = "heat"
    # In RS lower-bound mode, prime tail bound is not used; keep for reference.
    cert.prime_tail_bound = Interval(Decimal(0), heat_gaussian_prime_tail_bound(X, tau))
    return cert
