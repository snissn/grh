#!/usr/bin/env python3
# eisenstein_phase_validate.py  — public-facing report
#
# What we validate (paper → plain language):
#   Theorem 7.4 (H-form):  2 * sum_rho H(g; rho) = P(g) + J(g)
#   Under RH and summing only positive ordinates gamma>0:  4*S = P + J,  S = sum_{gamma>0} H(g; 1/2+i gamma)
#
# This program prints:
#   • WHAT is supposed to match (4*S vs P+J) and WHY (RH form of Thm 7.4)
#   • HOW each number was computed (finite sum vs numerical integral vs truncated zero-sum)
#   • EXACT vs APPROX status and approximation knobs (primes used, zeros used, last increment)
#   • A crisp VERDICT with residuals and convergence status
#
# References: Definitions of g, \hat g_b, Ĝ_b, and functionals P, J, H, A match the manuscript. (See Theorem 7.4, Prop. 7.1/7.3, Cor. 7.6.)

import argparse
import json
import sys
import time
from typing import Callable, List, Tuple, Optional, Dict

import mpmath as mp


# ---------- Pretty printing & verdict logic ----------

def nstr(x, dps=18):
    try:
        return mp.nstr(x, dps)
    except Exception:
        return str(x)

def verdict_from_residual(abs_resid, rel_resid, abs_tol, rel_tol, converged_flag):
    if abs_resid <= abs_tol or rel_resid <= rel_tol:
        return "PASS", "Residual is within tolerance."
    if not converged_flag:
        return "NEEDS_MORE_ZEROS", "Residual above tolerance and zero-sum not yet converged (increase --zeros-cap and/or precision)."
    return "FAIL", "Zero-sum converged but residual is above tolerance (increase precision or review parameters)."


# ---------- Test class ----------

def smooth_bump(u: mp.mpf) -> mp.mpf:
    # C^∞ bump on (-1,1)
    u = mp.mpf(u)
    if abs(u) < 1:
        return mp.e**(-1/(1-u*u))
    return mp.mpf('0')

class TestTransform:
    """
    Compactly supported cosine transform \hat g_b and derived transform \hat G_b:
      \hat g_b(x) = scale * φ(x/X) on [-X, X] (even)
      \hat G_b(x) = 2*sinh(|x|/2) * \hat g_b(2|x|)
    """
    def __init__(self, X: float, profile: Callable[[mp.mpf], mp.mpf] = smooth_bump, scale: float = 1.0):
        self.X = mp.mpf(X)
        self.profile = profile
        self.scale = mp.mpf(scale)

    def gb(self, x: mp.mpf) -> mp.mpf:
        x = mp.mpf(x)
        if abs(x) >= self.X:
            return mp.mpf('0')
        return self.scale * self.profile(x / self.X)

    def Gb(self, x: mp.mpf) -> mp.mpf:
        xx = abs(mp.mpf(x))
        return 2 * mp.sinh(xx / 2) * self.gb(2 * xx)


# ---------- P(g): finite prime–power sum (exact modulo FP) ----------

def sieve_primes_upto(n: int) -> List[int]:
    n = int(n)
    if n < 2: return []
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    r = int(n**0.5)
    for p in range(2, r + 1):
        if sieve[p]:
            step, start = p, p*p
            cnt = ((n - start) // step) + 1 if start <= n else 0
            if cnt > 0:
                sieve[start:n+1:step] = b"\x00" * cnt
    return [i for i in range(2, n + 1) if sieve[i]]

def compute_P_with_meta(tf: TestTransform) -> Tuple[mp.mpf, Dict[str, int]]:
    """
    P(g) = (1/pi) * sum_{p} sum_{k>=1} (log p) p^{-k/2} * \hat G_b(k log p).
    Finite by support: only terms with 2*k*log p < X survive.
    Returns (P, meta={#primes, #prime_powers, p_max}).
    """
    X = tf.X
    p_max = int(mp.floor(mp.e ** (X / 2)))  # k=1 bound
    primes = sieve_primes_upto(p_max)

    total = mp.mpf('0')
    count_pp = 0
    for p in primes:
        lp = mp.log(p)
        kmax = int(mp.floor(X / (2 * lp)))
        if kmax <= 0: continue
        for k in range(1, kmax + 1):
            total += (lp) * p ** (-mp.mpf(k) / 2) * tf.Gb(k * lp)
            count_pp += 1

    meta = {"num_primes": len(primes), "num_prime_powers": count_pp, "p_max": p_max}
    return total / mp.pi, meta


# ---------- J(g): numerical integral with a self-check error estimate ----------

def quad_with_split_estimate(func, a, b) -> Tuple[mp.mpf, mp.mpf]:
    """
    Compute integral on [a,b] two ways: single interval vs split midpoint.
    Return (value_single, |value_single - value_split|) as a crude error proxy.
    """
    val1 = mp.quad(func, [a, b])
    mid = (a + b) / 2
    val2 = mp.quad(func, [a, mid, b])
    return val1, abs(val1 - val2)

def compute_J_with_meta(tf: TestTransform) -> Tuple[mp.mpf, Dict[str, mp.mpf]]:
    """
    J(g) = (1/pi) ∫_0^{X/2} e^{-x} \hat g_b(2x) dx.
    Returns (J, meta={abs_error_est}).
    """
    upper = tf.X / 2
    f = lambda x: mp.e**(-x) * tf.gb(2 * x)
    val, err = quad_with_split_estimate(f, mp.mpf('0'), upper)
    return val / mp.pi, {"abs_error_est": err / mp.pi, "interval_upper": upper}


# ---------- Zero weights and RH sum with live residuals ----------

def H_weight(tf: TestTransform, beta: mp.mpf, gamma: mp.mpf) -> mp.mpf:
    upper = tf.X / 2
    integrand = lambda x: mp.e**(-beta * x) * tf.gb(2 * x) * mp.cos(gamma * x)
    return mp.quad(integrand, [0, upper]) / mp.pi

def sum_H_assuming_RH(tf: TestTransform,
                      n_zeros_cap: int,
                      tol: float,
                      batch: int,
                      progress: bool,
                      progress_every: int,
                      # for on-the-fly residual display:
                      target_value: Optional[mp.mpf] = None,
                      target_coeff: mp.mpf = mp.mpf('4')) -> Tuple[mp.mpf, Dict[str, mp.mpf]]:
    """
    S = sum_{gamma>0} H(g; 1/2 + i*gamma), adaptive in batches.
    If target_value is provided (e.g., P+J) prints current residual target_value - target_coeff * S_sofar.
    """
    S = mp.mpf('0')
    used = 0
    last_inc = mp.mpf('inf')
    converged = False

    def partial_sum(start: int, count: int) -> mp.mpf:
        subtot = mp.mpf('0')
        for n in range(start, start + count):
            rho = mp.zetazero(n)  # 0.5 + i gamma_n
            gamma = mp.im(rho)
            subtot += H_weight(tf, mp.mpf('0.5'), gamma)
        return subtot

    batch_index = 0
    while used < n_zeros_cap:
        to_take = min(batch, n_zeros_cap - used)
        inc = partial_sum(used + 1, to_take)
        S += inc
        used += to_take
        batch_index += 1
        last_inc = abs(inc)
        if progress and (batch_index % max(1, progress_every) == 0):
            if target_value is None:
                print(f"[zeros] used={used:5d}  |batch_inc|={nstr(last_inc)}  |S|={nstr(abs(S))}", flush=True)
            else:
                resid_now = target_value - target_coeff * S
                print(f"[zeros] used={used:5d}  |batch_inc|={nstr(last_inc)}  |S|={nstr(S)}  |residual so far|={nstr(abs(resid_now))}", flush=True)
        if last_inc < tol:
            converged = True
            break

    diags = {
        "zeros_used": used,
        "batch_size": batch,
        "last_batch_increment_abs": last_inc,
        "converged": converged,
        "tol": mp.mpf(tol),
    }
    return S, diags


# ---------- Archimedean A(g) (used in routes mode) ----------

_g_cache: Dict[Tuple[float, float, float], mp.mpf] = {}

def g_of_tau(tf: TestTransform, tau: mp.mpf) -> mp.mpf:
    key = (float(tf.X), float(tf.scale), float(tau))
    if key in _g_cache:
        return _g_cache[key]
    integrand = lambda x: tf.gb(x) * mp.cos(tau * x)
    val = (1/mp.pi) * mp.quad(integrand, [0, tf.X])
    _g_cache[key] = val
    return val

def compute_A_with_meta(tf: TestTransform, tau_max: float = 12.0, use_infinite_domain: bool = False) -> Tuple[mp.mpf, Dict[str, mp.mpf]]:
    K = lambda t: mp.re(mp.digamma(1j*t) - mp.digamma(mp.mpf('0.5') + 1j*t))
    integrand = lambda t: g_of_tau(tf, t) * K(t)
    if use_infinite_domain:
        val, err = quad_with_split_estimate(integrand, mp.mpf('0'), mp.inf)  # mp.quad supports [0, inf]
    else:
        val, err = quad_with_split_estimate(integrand, mp.mpf('0'), mp.mpf(tau_max))
    return (1/mp.pi) * val, {"abs_error_est": err / mp.pi, "tau_max": mp.mpf(tau_max), "infinite": bool(use_infinite_domain)}


# ---------- Public-facing printing ----------

def print_header_validate():
    print("\n" + "="*78)
    print("EISENSTEIN-PHASE EXPLICIT FORMULA — VALIDATION REPORT (RH specialization)")
    print("="*78)
    print("What should match (under RH, positive ordinates γ>0):  4 * S  ≈  P(g) + J(g)")
    print("S is the truncated sum over zeros S = Σ_{γ>0} H(g; 1/2 + i γ).")
    print("This is equivalent to Theorem 7.4 in the paper (Cor. 7.6 phrased for γ>0).")
    print("-"*78)

def print_inputs(X, mp_dps, zeros_cap, rh_tol, batch, gb_scale):
    print(f"Inputs:")
    print(f"  • Test support radius X           = {X}")
    print(f"  • Precision (mpmath, digits)      = {mp_dps}")
    print(f"  • Zero-sum cap / tol / batch      = {zeros_cap} / {rh_tol} / {batch}")
    print(f"  • Transform scale (for sensitivity)= {gb_scale}")
    print("-"*78)

def print_number(name, value, how, details, time_s=None):
    tstr = "" if time_s is None else f"  [time {time_s:.2f}s]"
    print(f"{name:<6} = {nstr(value)}{tstr}")
    print(f"   └─ how: {how}")
    if details:
        for k, v in details.items():
            print(f"      {k}: {v}")
    print()

def print_target_block(P, J, S, resid, rel, diags):
    print("-"*78)
    print("MAIN CHECK (public-facing):  does  4*S  match  P(g)+J(g) ?")
    print(f"  P(g) + J(g)     = {nstr(P + J)}   (target)")
    print(f"  4 * S_estimate  = {nstr(4*S)}   (from first {diags['zeros_used']} zeros)")
    print(f"  |difference|    = {nstr(abs(resid))}")
    print(f"  relative diff   = {nstr(rel)}")
    print()
    print("Interpretation:")
    print("  • If the last zero-sum batch increment is below the tolerance, the S estimate is")
    print("    numerically converged at this tolerance.")
    print("  • Otherwise, increase --zeros-cap and/or precision until the residual stabilizes.")
    print("-"*78)

def print_verdict(verdict, reason):
    print(f"VERDICT: {verdict} — {reason}")
    print("="*78 + "\n")

def print_timings(times_dict):
    print("Timings (seconds):")
    for k, v in times_dict.items():
        print(f"  {k:>8}: {v:.2f}")
    print("-"*78)

def print_next_steps(diags, resid, rel):
    print("Next steps (practical):")
    if not diags["converged"]:
        print("  • The zero-sum has NOT converged to the requested tolerance.")
        print("    Increase --zeros-cap and/or reduce --rh-tol (e.g., 1e-8 → 3e-9).")
        print("  • Larger batch sizes (--batch 150–250) can speed convergence per progress step.")
    else:
        print("  • The zero-sum is converged at the requested tolerance.")
        if rel > mp.mpf('1e-8'):
            print("  • Residual still large? Raise precision (--mp-dps) to reduce quadrature/FP error.")
    print()


# ---------- Validation mode ----------

def validate_corollary_RH_public(tf: TestTransform,
                                 mp_dps: int,
                                 n_zeros_cap: int,
                                 rh_tol: float,
                                 batch: int,
                                 resid_tol: float,
                                 abs_resid_tol: float,
                                 progress: bool,
                                 progress_every: int):

    print_header_validate()
    print_inputs(float(tf.X), mp_dps, n_zeros_cap, rh_tol, batch, float(tf.scale))

    # P(g)
    t0 = time.time()
    P, Pmeta = compute_P_with_meta(tf)
    tP = time.time()

    # J(g)
    J, Jmeta = compute_J_with_meta(tf)
    tJ = time.time()

    # S (with live residuals against P+J)
    target_PplusJ = P + J
    S, diags = sum_H_assuming_RH(tf,
                                 n_zeros_cap=n_zeros_cap, tol=rh_tol, batch=batch,
                                 progress=progress, progress_every=progress_every,
                                 target_value=target_PplusJ, target_coeff=mp.mpf('4'))
    tS = time.time()

    # Print components with provenance
    print_number("P(g)", P,
                 how="FINITE prime–power sum (exact modulo floating point).",
                 details={
                     "prime_upper_bound_p_max": Pmeta["p_max"],
                     "number_of_primes_used": Pmeta["num_primes"],
                     "number_of_prime_power_terms": Pmeta["num_prime_powers"],
                 },
                 time_s=tP - t0)

    print_number("J(g)", J,
                 how="NUMERICAL integral on [0, X/2] using mp.quad.",
                 details={
                     "interval_upper=X/2": nstr(Jmeta["interval_upper"]),
                     "abs_error_estimate": nstr(Jmeta["abs_error_est"]),
                 },
                 time_s=tJ - tP)

    print_number("S", S,
                 how="TRUNCATED sum over positive ζ-zeros: S = Σ_{γ>0} H(g; 1/2 + i γ).",
                 details={
                     "zeros_used": diags["zeros_used"],
                     "last_batch_increment_abs": nstr(diags["last_batch_increment_abs"]),
                     "converged_to_rh_tol": bool(diags["converged"]),
                     "rh_tol_requested": nstr(diags["tol"]),
                     "batch_size": diags["batch_size"],
                 },
                 time_s=tS - tJ)

    # Main residual
    resid = (P + J) - 4*S
    rel = abs(resid) / max(mp.mpf('1e-30'), abs(P + J))
    print_target_block(P, J, S, resid, rel, diags)

    verdict, reason = verdict_from_residual(abs(resid), rel, abs_resid_tol, mp.mpf(resid_tol), bool(diags["converged"]))
    print_verdict(verdict, reason)

    print_timings({"P": tP - t0, "J": tJ - tP, "Zero-sum S": tS - tJ})
    print_next_steps(diags, resid, rel)

    return {
        "P": P, "J": J, "S_pos": S,
        "Residual_main": resid, "RelResidual_main": rel,
        "Zeros_used": diags["zeros_used"],
        "Last_batch_increment_abs": diags["last_batch_increment_abs"],
        "Converged": diags["converged"],
        "Verdict": verdict,
        "Times": {"P_s": tP - t0, "J_s": tJ - tP, "S_s": tS - tJ}
    }


# ---------- Routes mode (optional, unchanged targets; prints that 4S=P+J implies I_A=I_B) ----------

def compute_I_route_A(tf: TestTransform, tau_max: float, arch_infty: bool):
    A, Ameta = compute_A_with_meta(tf, tau_max=tau_max, use_infinite_domain=arch_infty)
    P, Pmeta = compute_P_with_meta(tf)
    return (A - P), {"A": A, "P": P, "Ameta": Ameta, "Pmeta": Pmeta}

def compute_I_route_B_RH(tf: TestTransform, tau_max: float, arch_infty: bool,
                         n_zeros_cap: int, rh_tol: float, batch: int,
                         progress: bool, progress_every: int):
    A, Ameta = compute_A_with_meta(tf, tau_max=tau_max, use_infinite_domain=arch_infty)
    J, Jmeta = compute_J_with_meta(tf)
    S, diags = sum_H_assuming_RH(tf, n_zeros_cap=n_zeros_cap, tol=rh_tol, batch=batch,
                                 progress=progress, progress_every=progress_every)
    IB = A + J - 4*S
    meta = {"A": A, "J": J, "Ameta": Ameta, "Jmeta": Jmeta}
    meta.update(diags)
    return IB, meta

def routes_public(tf: TestTransform, args):
    print("\n" + "="*78)
    print("EISENSTEIN-PHASE — ROUTE A vs. ROUTE B (under RH)")
    print("="*78)
    print("Targets:")
    print("  • Route A: I_A = A(g) - P(g)")
    print("  • Route B: I_B = A(g) + J(g) - 4*S  (S=Σ_{γ>0} H)")
    print("  • Equality I_B - I_A = (P+J) - 4S, same residual as main check.")
    print("-"*78)
    print_inputs(float(tf.X), mp.mp.dps, args.zeros_cap, args.rh_tol, args.batch, float(tf.scale))

    t0 = time.time()
    IA, left = compute_I_route_A(tf, args.tau_max, args.arch_infty)
    tA = time.time()
    IB, right = compute_I_route_B_RH(tf, args.tau_max, args.arch_infty,
                                     args.zeros_cap, args.rh_tol, args.batch,
                                     args.progress, args.progress_every)
    tB = time.time()

    P, A = left["P"], left["A"]
    J, S = right["J"], right.get("S_RH", None)
    deltaI = IB - IA
    resid = (P + J) - 4*S
    rel = abs(resid) / max(mp.mpf('1e-30'), abs(P + J))
    verdict, reason = verdict_from_residual(abs(resid), rel, args.abs_resid_tol, mp.mpf(args.resid_tol), bool(right["converged"]))

    print_number("A(g)", A, "NUMERICAL integral in τ (see meta).", details=left["Ameta"], time_s=tA - t0)
    print_number("P(g)", P, "FINITE prime–power sum (exact modulo FP).", details=left["Pmeta"])
    print_number("J(g)", J, "NUMERICAL integral on [0, X/2].", details=right["Jmeta"])
    print_number("S", right.get("S_RH"), "TRUNCATED sum over positive ζ-zeros.", {
        "zeros_used": right["zeros_used"],
        "last_batch_increment_abs": nstr(right["last_batch_increment_abs"]),
        "converged_to_rh_tol": bool(right["converged"]),
    }, time_s=tB - tA)

    print("-"*78)
    print(f"I_A = A - P           = {nstr(IA)}")
    print(f"I_B = A + J - 4*S     = {nstr(IB)}")
    print(f"ΔI   = I_B - I_A      = {nstr(deltaI)}  (should be 0 if 4S=P+J)")
    print_target_block(P, J, right.get("S_RH"), resid, rel, right)
    print_verdict(verdict, reason)
    print_timings({"IA": tA - t0, "IB": tB - tA})


# ---------- Smoke test (Appendix D(ii)) ----------

def smoke_public(args):
    Xcrit = 2 * mp.log(2)
    X = float(0.9 * Xcrit)  # strict inclusion
    tf = TestTransform(X=X)

    print("\n" + "="*78)
    print("ZERO-SUPPORT SMOKE TEST — primes silent: P(g) ≈ 0, check 4*S ≈ J(g)")
    print("="*78)
    print_inputs(X, args.mp_dps, args.zeros_cap, args.rh_tol, args.batch, 1.0)

    t0 = time.time()
    P, Pmeta = compute_P_with_meta(tf)
    tP = time.time()
    J, Jmeta = compute_J_with_meta(tf)
    tJ = time.time()
    S, diags = sum_H_assuming_RH(tf, n_zeros_cap=args.zeros_cap, tol=args.rh_tol, batch=args.batch,
                                 progress=args.progress, progress_every=args.progress_every)
    tS = time.time()

    print_number("P(g)", P, "FINITE prime–power sum (should be ~0 by design).", Pmeta, time_s=tP - t0)
    print_number("J(g)", J, "NUMERICAL integral on [0, X/2].", Jmeta, time_s=tJ - tP)
    print_number("S", S, "TRUNCATED sum over positive ζ-zeros.", {
        "zeros_used": diags["zeros_used"],
        "last_batch_increment_abs": nstr(diags["last_batch_increment_abs"]),
        "converged_to_rh_tol": bool(diags["converged"]),
    }, time_s=tS - tJ)

    resid = J - 4*S
    rel = abs(resid) / max(mp.mpf('1e-30'), abs(J))
    verdict, reason = verdict_from_residual(abs(resid), rel, args.abs_resid_tol, mp.mpf(args.resid_tol), bool(diags["converged"]))
    print_target_block(P=mp.mpf('0'), J=J, S=S, resid=resid, rel=rel, diags=diags)
    print_verdict(verdict, reason)
    print_timings({"P": tP - t0, "J": tJ - tP, "Zero-sum S": tS - tJ})


# ---------- CLI & presets ----------

def apply_preset(args):
    if not args.preset: return
    if args.preset == "quick":
        args.X = args.X or 6.0
        args.mp_dps = 40
        args.zeros_cap = 300
        args.rh_tol = 1e-7
        args.batch = 150
    elif args.preset == "standard":
        args.X = args.X or 6.0
        args.mp_dps = 50
        args.zeros_cap = 900
        args.rh_tol = 5e-9
        args.batch = 120
    elif args.preset == "deep":
        args.X = args.X or 8.0
        args.mp_dps = 60
        args.zeros_cap = 2000
        args.rh_tol = 1e-10
        args.batch = 100

def main():
    ap = argparse.ArgumentParser(
        description="Validates 4*S ?= P+J (RH specialization of Theorem 7.4) with public-facing output."
    )
    ap.add_argument("--mode", choices=["validate", "routes", "smoke"], default="validate",
                    help="validate: RH specialization 4*S ?= P+J; routes: compare I_A vs I_B; smoke: zero-support test.")
    ap.add_argument("--preset", choices=["quick", "standard", "deep"],
                    help="Preset for speed/accuracy knobs.")
    ap.add_argument("--X", type=float, default=6.0, help="Support radius for ĝ_b (supp in [-X, X]).")
    ap.add_argument("--mp-dps", type=int, default=40, help="mpmath precision (decimal digits).")
    ap.add_argument("--zeros-cap", dest="zeros_cap", type=int, default=600, help="Max positive-line zeros (adaptive stops earlier if converged).")
    ap.add_argument("--rh-tol", type=float, default=1e-8, help="Adaptive batch tolerance for zero-sum S.")
    ap.add_argument("--batch", type=int, default=100, help="Batch size for zero-sum adaptation.")
    ap.add_argument("--gb-scale", type=float, default=1.0, help="Scale factor for ĝ_b.")
    ap.add_argument("--tau-max", type=float, default=12.0, help="Upper limit for τ in A(g) on [0, τ_max] (routes mode).")
    ap.add_argument("--arch-infty", action="store_true", help="Integrate A(g) over [0, ∞) (routes mode; slower).")
    ap.add_argument("--resid-tol", type=float, default=1e-8, help="Relative residual tolerance for verdicts.")
    ap.add_argument("--abs-resid-tol", type=float, default=1e-12, help="Absolute residual tolerance for tiny targets.")
    ap.add_argument("--p-zero-tol", type=float, default=1e-20, help="Smoke test pass threshold for |P(g)|~0.")
    ap.add_argument("--json-out", type=str, default="", help="Optional path to write a JSON summary.")
    ap.add_argument("--progress", action="store_true", help="Print progress lines during the zero-sum.")
    ap.add_argument("--progress-every", type=int, default=1, help="Print progress after every N batches.")
    args = ap.parse_args()

    # Apply preset overrides
    apply_preset(args)

    mp.mp.dps = int(args.mp_dps)
    print(args)

    if args.mode == "smoke":
        smoke_public(args)
    elif args.mode == "routes":
        tf = TestTransform(X=args.X, scale=args.gb_scale)
        routes_public(tf, args)
    else:
        tf = TestTransform(X=args.X, scale=args.gb_scale)
        out = validate_corollary_RH_public(
            tf, mp_dps=args.mp_dps, n_zeros_cap=args.zeros_cap, rh_tol=args.rh_tol, batch=args.batch,
            resid_tol=args.resid_tol, abs_resid_tol=args.abs_resid_tol,
            progress=args.progress, progress_every=args.progress_every
        )
        if args.json_out:
            try:
                with open(args.json_out, "w", encoding="utf-8") as f:
                    json.dump({k: str(v) for k, v in out.items()}, f, indent=2)
                print(f"\nWrote JSON summary to: {args.json_out}", flush=True)
            except Exception as e:
                print(f"\n[warn] Failed to write JSON summary: {e}", flush=True)

if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:
        pass
    main()

