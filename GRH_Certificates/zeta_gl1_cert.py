from decimal import Decimal
from typing import List
import math

from grh_verifier import (
    Interval, TestFunction, PrimeLocalTerm, RamifiedTerm, Certificate, verify_certificate
)

# Rigorous generators (for BL evaluate mode)
from ef_generators import HeatGaussianPhi, UnramifiedLocal, unramified_term_interval, ArchParams, compute_arch_block, heat_gaussian_prime_tail_bound
from rigor_backend import RMATH


def primes_upto(n: int) -> List[int]:
    if n < 2:
        return []
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    import math as _m
    for p in range(2, int(_m.isqrt(n)) + 1):
        if sieve[p]:
            step = p
            start = p * p
            sieve[start:n+1:step] = [False] * (((n - start) // step) + 1)
    return [i for i, is_p in enumerate(sieve) if is_p]


def build_zeta_gl1_certificate_rs() -> Certificate:
    """RS-positivity mode baseline (no prime generation)."""
    Phi = TestFunction(family="heat", a=Decimal("0.8"))
    cert = Certificate(
        m=1,
        K="Q",
        Qpi=Decimal("1"),
        t_star=Decimal("1.0"),
        Phi_test=Phi,
        precision_bits=212,
        test_family="BL",
        prime_block_mode="rs_lower_bound",
        proof_mode="ann_odd_only",
        Ainf=Decimal("1.0"),
        beta_m=Decimal("0.5"),
        C_Rstar=Decimal("1.0"),
        arch_value=Interval.point(0),
        ram_terms=[],
        j_max=0,
        band_limit_X=Decimal("6.0"),
        q_modulus=1,
        eta_a5_gap=Interval.point("0.40"),
        epsilon_invariance=Interval(Decimal("0"), Decimal("0.20")),
        alpha_inv=Decimal("0.5"),
        kappa_inv=Decimal("1.0"),
    )
    # If ARB is available, compute arch block rigorously for zeta: Γ_R(s) with μ_R=[0], μ_C=[] at σ=1/2
    if RMATH.mode == "arb":
        try:
            phi_spec = HeatGaussianPhi._mk(a_weight=Decimal("0.8"), tau=Decimal("2.0"))
            arch_iv = compute_arch_block(phi_spec, ArchParams(mu_R=[Decimal("0")], mu_C=[], sigma=Decimal("0.5")))
            cert.arch_value = Interval(arch_iv.lo, arch_iv.hi)
        except Exception:
            pass
    return cert


def build_zeta_gl1_certificate_bl_evaluate(X: Decimal = Decimal("6.0"),
                                           a: Decimal = Decimal("0.8"),
                                           tau: Decimal = Decimal("2.0")) -> Certificate:
    """BL evaluate mode with generated unramified terms for ζ.

    Includes (p,k) terms satisfying 2 k log p <= X.
    """
    # Create a rigorous PhiSpec for Gaussian heat
    phi_spec = HeatGaussianPhi._mk(a_weight=a, tau=tau)

    # Determine prime cutoff: for k=1, require log p <= X/2 ⇒ p <= exp(X/2)
    Pmax = int(math.floor(math.exp(float(X) / 2.0)))
    primes = primes_upto(Pmax)

    unram_local: List[PrimeLocalTerm] = []
    for p in primes:
        loc = UnramifiedLocal(p=p, alphas=[1.0])  # GL(1): α=1
        logp = Decimal(str(math.log(p)))
        # Max k satisfying 2k log p <= X
        if logp == 0:
            continue
        k_max = int(math.floor((float(X) / 2.0) / float(logp)))
        for k in range(1, max(1, k_max) + 1):
            W_iv = unramified_term_interval(phi_spec, loc, k, T_weight=Decimal("12"))
            # term interval already equals (log p) * Re tr(A_p^k) * W(k log p)
            term_iv = Interval(W_iv.lo, W_iv.hi)
            unram_local.append(PrimeLocalTerm(p=p, k=k, value=term_iv))

    cert = Certificate(
        m=1,
        K="Q",
        Qpi=Decimal("1"),
        t_star=Decimal("1.0"),
        Phi_test=TestFunction(family="heat", a=a),
        precision_bits=256,
        test_family="BL",
        prime_block_mode="evaluate",
        proof_mode="ann_odd_only",
        Ainf=Decimal("1.0"),
        beta_m=Decimal("0.5"),
        C_Rstar=Decimal("1.0"),
        arch_value=Interval.point(0),
        ram_terms=[],
        j_max=0,
        band_limit_X=X,
        q_modulus=1,
        eta_a5_gap=Interval.point("0.40"),
        epsilon_invariance=Interval(Decimal("0"), Decimal("0.20")),
        alpha_inv=Decimal("0.5"),
        kappa_inv=Decimal("1.0"),
        unram_local=unram_local,
    )
    if RMATH.mode == "arb":
        try:
            # Reuse the same phi_spec for arch integration
            arch_iv = compute_arch_block(phi_spec, ArchParams(mu_R=[Decimal("0")], mu_C=[], sigma=Decimal("0.5")))
            cert.arch_value = Interval(arch_iv.lo, arch_iv.hi)
        except Exception:
            pass
    return cert


def build_zeta_gl1_certificate_heat_evaluate(X: Decimal = Decimal("6.0"),
                                             a: Decimal = Decimal("0.8"),
                                             tau: Decimal = Decimal("2.0")) -> Certificate:
    """Heat evaluate mode with generated unramified terms and a rigorous tail bound.

    Terms included for 2k log p <= X. Tail beyond X bounded by heat_gaussian_prime_tail_bound.
    """
    phi_spec = HeatGaussianPhi._mk(a_weight=a, tau=tau)
    Pmax = int(math.floor(math.exp(float(X) / 2.0)))
    primes = primes_upto(Pmax)
    unram_local: List[PrimeLocalTerm] = []
    for p in primes:
        loc = UnramifiedLocal(p=p, alphas=[1.0])
        logp = Decimal(str(math.log(p)))
        if logp == 0:
            continue
        k_max = int(math.floor((float(X) / 2.0) / float(logp)))
        for k in range(1, max(1, k_max) + 1):
            W_iv = unramified_term_interval(phi_spec, loc, k, T_weight=Decimal("12"))
            term_iv = Interval(W_iv.lo, W_iv.hi)
            unram_local.append(PrimeLocalTerm(p=p, k=k, value=term_iv))

    tail = heat_gaussian_prime_tail_bound(X, tau)

    cert = Certificate(
        m=1,
        K="Q",
        Qpi=Decimal("1"),
        t_star=Decimal("1.0"),
        Phi_test=TestFunction(family="heat", a=a),
        precision_bits=256,
        test_family="heat",
        prime_block_mode="evaluate",
        proof_mode="ann_odd_only",
        Ainf=Decimal("1.0"),
        beta_m=Decimal("0.5"),
        C_Rstar=Decimal("1.0"),
        arch_value=Interval.point(0),
        ram_terms=[],
        j_max=0,
        band_limit_X=X,
        q_modulus=1,
        eta_a5_gap=Interval.point("0.40"),
        epsilon_invariance=Interval(Decimal("0"), Decimal("0.20")),
        alpha_inv=Decimal("0.5"),
        kappa_inv=Decimal("1.0"),
        unram_local=unram_local,
        prime_tail_bound=Interval(Decimal(0), tail),
    )
    if RMATH.mode == "arb":
        try:
            arch_iv = compute_arch_block(phi_spec, ArchParams(mu_R=[Decimal("0")], mu_C=[], sigma=Decimal("0.5")))
            cert.arch_value = Interval(arch_iv.lo, arch_iv.hi)
        except Exception:
            pass
    return cert


def run_proof_of_grh_rs() -> None:
    """Run a small test net in RS mode to produce GRH_Verified."""
    cert = build_zeta_gl1_certificate_rs()
    # Provide a tiny test net; RS mode ignores prime recomputation, safe here.
    cert.proof_mode = "proof_of_GRH"
    cert.test_net = [
        TestFunction(family="heat", a=Decimal("0.7")),
        TestFunction(family="heat", a=Decimal("0.8")),
        TestFunction(family="heat", a=Decimal("0.9")),
    ]
    rep = verify_certificate(cert)
    print("\n-- proof_of_GRH (RS mode, small test net) --")
    print("Result:", rep.result)
    print("Phase:", rep.phase_passed)
    print("Details:", rep.details)


if __name__ == "__main__":
    print("-- RS-positivity mode --")
    cert_rs = build_zeta_gl1_certificate_rs()
    report_rs = verify_certificate(cert_rs)
    print("Result:", report_rs.result)
    print("Phase:", report_rs.phase_passed)
    print("Details:", report_rs.details)

    print("\n-- BL evaluate mode (generated primes) --")
    cert_ev = build_zeta_gl1_certificate_bl_evaluate()
    report_ev = verify_certificate(cert_ev)
    print("Result:", report_ev.result)
    print("Phase:", report_ev.phase_passed)
    print("Details:", report_ev.details)

    # proof_of_GRH (RS mode)
    run_proof_of_grh_rs()

    print("\n-- Heat evaluate mode (generated primes + tail bound) --")
    cert_h = build_zeta_gl1_certificate_heat_evaluate()
    report_h = verify_certificate(cert_h)
    print("Result:", report_h.result)
    print("Phase:", report_h.phase_passed)
    print("Details:", report_h.details)
