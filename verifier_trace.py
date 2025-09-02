from __future__ import annotations
from dataclasses import asdict
from decimal import Decimal
from typing import Dict, Any, List, Optional

from grh_verifier import (
    Interval, Certificate, TestFunction,
    ArchBlock, RamBlock, PrimeBlockBL, PrimeTailBound_heat,
    RS_Positivity_LB, SANorm, RStarBudget,
    Q_of_X_eps, Eps_from_q_X, CheckA5Gap,
    set_precision_bits, verify_certificate,
)


def compute_trace(cert: Certificate) -> Dict[str, Any]:
    """Recompute key quantities to visualize and document the verification.

    Returns a dict with intervals and numeric summaries, plus the VerifierReport.
    """
    # Set precision for consistency
    set_precision_bits(cert.precision_bits)

    trace: Dict[str, Any] = {
        "params": {
            "m": cert.m, "K": cert.K, "Qpi": str(cert.Qpi), "t_star": str(cert.t_star),
            "family": cert.test_family, "prime_block_mode": cert.prime_block_mode,
            "proof_mode": cert.proof_mode,
        }
    }

    # Phase 1: SA norm
    M = SANorm(cert.Phi_test, a_weight=Decimal(cert.Phi_test.a or 1),
               t_min=Decimal(0), t_zero=cert.t_star, precision_bits=cert.precision_bits)
    trace["SA_norm"] = (str(M.lo), str(M.hi))

    # Phase 2: arch, ram, prime
    arch_iv = ArchBlock(cert.t_star, cert.Phi_test, cert.arch_value)
    ram_iv = RamBlock(cert.t_star, cert.Phi_test, cert.ram_terms, cert.j_max)
    trace["arch_iv"] = (str(arch_iv.lo), str(arch_iv.hi))
    trace["ram_iv"] = (str(ram_iv.lo), str(ram_iv.hi))

    if cert.prime_block_mode == "evaluate":
        if cert.test_family == "BL":
            prime_iv = PrimeBlockBL(cert.t_star, cert.Phi_test, cert.unram_local, cert.band_limit_X or Decimal(0))
            tail_iv = Interval.point(0)
        else:
            prime_iv = Interval.point(0)
            for term in cert.unram_local:
                if (cert.p_max is None or term.p <= (cert.p_max)) and (cert.k_max is None or term.k <= (cert.k_max)):
                    prime_iv = prime_iv + term.value
            tail_iv = PrimeTailBound_heat(cert.prime_tail_bound)
            prime_iv = prime_iv + tail_iv
        LB_Prime = prime_iv
    else:
        LB_Prime = RS_Positivity_LB(cert.Phi_test)
        tail_iv = Interval.point(0)

    trace["prime_iv"] = (str(LB_Prime.lo), str(LB_Prime.hi))
    trace["prime_tail"] = (str(tail_iv.lo), str(tail_iv.hi))

    arch_plus_ram = arch_iv + ram_iv
    trace["arch_plus_ram"] = (str(arch_plus_ram.lo), str(arch_plus_ram.hi))

    if cert.prime_block_mode == "evaluate":
        geom = LB_Prime - arch_plus_ram
        trace["geom_amp"] = (str(geom.lo), str(geom.hi))
    else:
        # RS: check 0 - UB(Arch+Ram)
        geom = Interval.point(0) - Interval.point(arch_plus_ram.hi)
        trace["geom_lb_rs"] = (str(geom.lo), str(geom.hi))

    # Phase 3: R* budget
    B_Rstar = RStarBudget(M, cert.Ainf, cert.Qpi, cert.beta_m, cert.C_Rstar, cert.t_star)
    trace["B_Rstar"] = (str(B_Rstar.lo), str(B_Rstar.hi))

    # Phase 5: annihilation
    sqrt1m_eta = None
    eps_iv = None
    if cert.eta_a5_gap is not None:
        sqrt1m_eta = CheckA5Gap(cert.eta_a5_gap)
        trace["sqrt_1m_eta"] = (str(sqrt1m_eta.lo), str(sqrt1m_eta.hi))
    if cert.epsilon_invariance is not None:
        eps_iv = cert.epsilon_invariance
        if cert.q_modulus is not None and cert.band_limit_X is not None and cert.alpha_inv is not None and cert.kappa_inv is not None:
            eps_cap = Eps_from_q_X(Decimal(cert.q_modulus), Decimal(cert.band_limit_X), Decimal(cert.alpha_inv), Decimal(cert.kappa_inv))
            eps_iv = Interval(Decimal(0), min(eps_iv.hi, eps_cap.hi))
        trace["epsilon_iv"] = (str(eps_iv.lo), str(eps_iv.hi))
    if sqrt1m_eta is not None and eps_iv is not None:
        annih_sum_hi = sqrt1m_eta.hi + eps_iv.hi
        trace["annih_sum_hi"] = str(annih_sum_hi)

    # Attach unramified term summaries if present
    if cert.unram_local:
        terms: List[Dict[str, Any]] = []
        for t in cert.unram_local:
            terms.append({"p": t.p, "k": t.k, "lo": str(t.value.lo), "hi": str(t.value.hi)})
        trace["unram_terms"] = terms

    # Now run the real verifier
    report = verify_certificate(cert)
    trace["report"] = {"result": report.result, "phase_passed": report.phase_passed, "details": report.details}

    return trace

