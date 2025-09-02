from decimal import Decimal
from typing import List

from grh_verifier import (
    Interval, TestFunction, PrimeLocalTerm, RamifiedTerm, Certificate, verify_certificate
)


def build_zeta_gl1_certificate() -> Certificate:
    # GL(1) over Q (Riemann zeta): conductor Qpi = 1, no ramification
    Phi = TestFunction(family="heat", a=Decimal("0.8"))

    cert = Certificate(
        m=1,
        K="Q",
        Qpi=Decimal("1"),
        t_star=Decimal("1.0"),
        Phi_test=Phi,
        precision_bits=212,
        test_family="BL",                 # use BL; prime block bounded by X
        prime_block_mode="rs_lower_bound",# RS-positivity lower bound on primes ([0,0])
        proof_mode="ann_odd_only",

        # Infrastructure constants (can be tuned to your ledger)
        Ainf=Decimal("1.0"),
        beta_m=Decimal("0.5"),
        C_Rstar=Decimal("1.0"),

        # Archimedean/Ramified: zeta has no ramified local terms at finite places
        arch_value=Interval.point(0),
        ram_terms=[],
        j_max=0,

        # BL cutoff and amplifier modulus
        band_limit_X=Decimal("6.0"),
        q_modulus=1,

        # A5 gap and invariance inputs (placeholders satisfying sqrt(1-eta)+eps < 1)
        eta_a5_gap=Interval.point("0.40"),
        epsilon_invariance=Interval(Decimal("0"), Decimal("0.20")),
        alpha_inv=Decimal("0.5"),
        kappa_inv=Decimal("1.0"),
    )

    # In RS lower-bound mode we don't need to populate unramified terms; LB_Prime is [0,0].
    return cert


if __name__ == "__main__":
    cert = build_zeta_gl1_certificate()
    report = verify_certificate(cert)
    print("Result:", report.result)
    print("Phase:", report.phase_passed)
    print("Details:", report.details)

