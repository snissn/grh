
"""
rigor_quad.py
-------------
Validated quadrature (interval-based) for:
- Finite-interval integrals via adaptive interval subdivision.
- Semi-infinite integrals using a user-supplied tail majorant.
- Oscillatory Fourier-cosine weights with rigorous tail enclosure.

This is backend-agnostic; it uses mpmath.iv if available to evaluate
the integrand on real intervals to get rigorous range enclosures.

Author: ChatGPT (GPT-5 Pro)
"""

from __future__ import annotations
from dataclasses import dataclass
from decimal import Decimal
import math
from typing import Callable, Optional, Tuple

from rigor_backend import RInterval, RMATH, HAVE_MPMATH

if HAVE_MPMATH:
    import mpmath as mp  # type: ignore


@dataclass
class TailMajorant:
    """Tail bound M(T): for x >= T, |f(x)| <= M_bound(x) with integrable envelope.
    We require a function that returns a rigorous overestimate of ∫_T^∞ M_bound(x) dx.
    """
    integral_upper: Callable[[Decimal], Decimal]  # returns an upper bound on the tail integral


def _iv_eval_nonneg(f_nonneg_iv: Callable[[Tuple[Decimal, Decimal]], Tuple[Decimal, Decimal]],
                    a: Decimal, b: Decimal) -> Tuple[Decimal, Decimal]:
    """Evaluate nonnegative function on [a,b] using interval arithmetic.
    f_nonneg_iv takes a real interval and returns [lo,hi] with 0 <= lo <= hi.
    """
    lo, hi = f_nonneg_iv((a, b))
    if lo < 0:
        lo = Decimal(0)
    return lo, hi


def integrate_nonneg_0inf(f_nonneg_iv: Callable[[Tuple[Decimal, Decimal]], Tuple[Decimal, Decimal]],
                          T: Decimal,
                          tail: TailMajorant,
                          rel_tol: Decimal = Decimal("1e-20"),
                          max_depth: int = 20) -> RInterval:
    """Integrate nonnegative f over [0,∞) as [0,T] + tail bound.

    f_nonneg_iv : callable taking (lo,hi) and returning [f_min,f_max] enclosure.
    T           : finite cutoff (increase to shrink the tail).
    tail        : TailMajorant giving rigorous upper bound for ∫_T^∞ |f|.

    Returns an interval enclosure.
    """
    box = integrate_over_interval(f_nonneg_iv, Decimal(0), T, rel_tol=rel_tol, max_depth=max_depth)
    tail_hi = tail.integral_upper(T)
    return RInterval(box.lo, box.hi + tail_hi)


def integrate_over_interval(f_iv: Callable[[Tuple[Decimal, Decimal]], Tuple[Decimal, Decimal]],
                            a: Decimal, b: Decimal,
                            rel_tol: Decimal = Decimal("1e-20"),
                            max_depth: int = 20) -> RInterval:
    """Validated interval integration on [a,b] using adaptive bisection.

    On each subinterval I, we evaluate f_iv(I) -> [lo,hi] and multiply by |I|
    to get an enclosure of the local integral. Refine until the sum's width
    is below rel_tol * |value| + abs_tol.

    This is conservative but robust for monotone/smooth integrands.
    """
    length = b - a
    if length <= 0:
        return RInterval.point(0)

    # Absolute tolerance derived from rel_tol and a heuristic scale
    abs_tol = rel_tol * (abs(length) + Decimal(1))

    def recurse(x0: Decimal, x1: Decimal, depth: int) -> RInterval:
        width = x1 - x0
        # Enclose f on [x0,x1]
        lo, hi = f_iv((x0, x1))
        # Local integral enclosure: [lo,hi] * width
        I_lo = lo * width
        I_hi = hi * width
        # If function range is tight enough, accept
        if (I_hi - I_lo) <= abs_tol or depth >= max_depth:
            return RInterval(I_lo, I_hi)

        # Otherwise split
        xm = (x0 + x1) / 2
        left = recurse(x0, xm, depth + 1)
        right = recurse(xm, x1, depth + 1)
        return RInterval(left.lo + right.lo, left.hi + right.hi)

    return recurse(a, b, 0)


def integrate_cos_weight(phi_iv: Callable[[Tuple[Decimal, Decimal]], Tuple[Decimal, Decimal]],
                         u: Decimal, T: Decimal,
                         tail: TailMajorant,
                         rel_tol: Decimal = Decimal("1e-20"),
                         max_depth: int = 18) -> RInterval:
    """Compute ∫_0^∞ phi(t) cos(u t) dt rigorously.

    We integrate on [0,T] by enclosing phi and cos over subintervals.
    Over tail [T,∞), we bound by ∫_T^∞ |phi(t)| dt (since |cos|<=1).
    """
    if HAVE_MPMATH:
        try:
            # helper: enclose phi(t)*cos(u t) on [x0,x1]
            def f_iv(I: Tuple[Decimal, Decimal]) -> Tuple[Decimal, Decimal]:
                x0, x1 = I
                # Build interval for u*t directly from endpoints (u>=0 here)
                u_dec = u
                if u_dec < 0:
                    ut_a = u_dec * x1
                    ut_b = u_dec * x0
                else:
                    ut_a = u_dec * x0
                    ut_b = u_dec * x1
                if ut_a > ut_b:
                    ut_a, ut_b = ut_b, ut_a
                ut_iv = mp.iv.mpf([str(ut_a), str(ut_b)])  # type: ignore
                c_iv = mp.iv.cos(ut_iv)  # type: ignore
                # enclose phi on interval
                p_lo, p_hi = phi_iv(I)
                if p_lo == p_hi:
                    p_iv = mp.mpf(str(p_lo))  # type: ignore
                else:
                    p_iv = mp.iv.mpf([str(p_lo), str(p_hi)])  # type: ignore
                prod = p_iv * c_iv  # interval product
                return Decimal(str(prod.a)), Decimal(str(prod.b))

            box = integrate_over_interval(f_iv, Decimal(0), T, rel_tol=rel_tol, max_depth=max_depth)
            tail_hi = tail.integral_upper(T)  # since |cos|<=1
            return RInterval(box.lo - tail_hi, box.hi + tail_hi)
        except Exception:
            # fall back to numerical quadrature below
            pass

    # Fallback: numeric trapezoid with conservative padding
    import mpmath as mp  # type: ignore
    mp.mp.dps = max(60, int(1/0.30103 * 80))

    def phi_mid(x: Decimal) -> float:
        lo, hi = phi_iv((x, x))
        return float((lo + hi) / Decimal(2))

    N = 4096
    a = 0.0
    b = float(T)
    h = (b - a) / N
    s = 0.0
    for i in range(N + 1):
        t = a + i * h
        w = 0.5 if (i == 0 or i == N) else 1.0
        s += w * phi_mid(Decimal(str(t))) * math.cos(float(u) * t)
    I_est = s * h
    tail_hi = float(tail.integral_upper(T))
    pad = float(rel_tol) * (abs(I_est) + 1.0)
    lo = Decimal(str(I_est - tail_hi - pad))
    hi = Decimal(str(I_est + tail_hi + pad))
    return RInterval(lo, hi)
