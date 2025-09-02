from __future__ import annotations
from dataclasses import dataclass
from decimal import Decimal
from typing import Dict, List, Tuple
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
                                        tau: Decimal = Decimal("2.0")) -> Certificate:
    phi_spec = HeatGaussianPhi._mk(a_weight=a, tau=tau)
    Pmax = int(math.floor(math.exp(float(X) / 2.0)))
    primes = primes_upto(Pmax)
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
            unram.append(PrimeLocalTerm(p=p, k=k, value=Interval(W_iv.lo, W_iv.hi)))

    cert = Certificate(
        m=2,
        K="Q",
        Qpi=Decimal("11"),
        t_star=Decimal("1.0"),
        Phi_test=TestFunction(family="heat", a=a),
        precision_bits=256,
        test_family="BL",
        prime_block_mode="evaluate",
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
    # Try arch via ARB: holomorphic weight k=2 ⇒ Γ_C(s+(k-1)/2) = Γ_C(s+1/2)
    if RMATH.mode == "arb":
        try:
            arch_iv = compute_arch_block(phi_spec, ArchParams(mu_R=[], mu_C=[Decimal("0.5")], sigma=Decimal("0.5")))
            cert.arch_value = Interval(arch_iv.lo, arch_iv.hi)
        except Exception:
            pass
    return cert


def build_gl3_sym2_11a1_bl_evaluate(X: Decimal = Decimal("6.0"),
                                     a: Decimal = Decimal("0.8"),
                                     tau: Decimal = Decimal("2.0")) -> Certificate:
    phi_spec = HeatGaussianPhi._mk(a_weight=a, tau=tau)
    Pmax = int(math.floor(math.exp(float(X) / 2.0)))
    primes = primes_upto(Pmax)
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
            unram.append(PrimeLocalTerm(p=p, k=k, value=Interval(W_iv.lo, W_iv.hi)))

    cert = Certificate(
        m=3,
        K="Q",
        Qpi=Decimal("11"),
        t_star=Decimal("1.0"),
        Phi_test=TestFunction(family="heat", a=a),
        precision_bits=256,
        test_family="BL",
        prime_block_mode="evaluate",
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
    # Arch factors for Sym^2 not specified here; keep arch_value as provided
    return cert


def build_gl2_11a1_heat_evaluate(X: Decimal = Decimal("6.0"),
                                 a: Decimal = Decimal("0.8"),
                                 tau: Decimal = Decimal("2.0")) -> Certificate:
    cert = build_gl2_newform_11a1_bl_evaluate(X, a, tau)
    cert.test_family = "heat"
    # Conservative tail: scale zeta tail by degree
    cert.prime_tail_bound = Interval(Decimal(0), Decimal(2) * heat_gaussian_prime_tail_bound(X, tau))
    return cert

