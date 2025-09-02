
"""
rigor_backend.py
----------------
Numeric backends for rigorous interval computations.

Priority of backends:
1) python-flint (ARB/ACB ball arithmetic) — if available.
2) mpmath interval arithmetic (iv.mpf) — portable, rigorous within mp.dps.
3) Decimal fallback — non-rigorous, last resort.

This module exposes a minimal facade used by the rest of the system.

Author: ChatGPT (GPT-5 Pro)
"""

from __future__ import annotations
from dataclasses import dataclass
from decimal import Decimal
from typing import Callable, Optional, Tuple, Union

# -----------------------------
# Attempt to import python-flint
# -----------------------------
HAVE_ARB = False
try:
    # python-flint typically installs as "flint"
    from flint import arb, acb, ctx  # type: ignore
    HAVE_ARB = True  # success
except Exception:  # pragma: no cover
    HAVE_ARB = False

# -----------------------------
# Attempt to import mpmath
# -----------------------------
HAVE_MPMATH = False
try:
    import mpmath as mp  # type: ignore
    HAVE_MPMATH = True
except Exception:  # pragma: no cover
    HAVE_MPMATH = False


# -----------------------------
# Utility
# -----------------------------
NumberLike = Union[int, float, Decimal]


@dataclass
class RInterval:
    """Real interval [lo, hi]."""
    lo: Decimal
    hi: Decimal

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError(f"Invalid RInterval: lo={self.lo} > hi={self.hi}")

    @staticmethod
    def point(x: NumberLike) -> "RInterval":
        d = Decimal(str(x)) if not isinstance(x, Decimal) else x
        return RInterval(d, d)

    def to_tuple(self) -> Tuple[Decimal, Decimal]:
        return (self.lo, self.hi)


# -----------------------------
# Facade for rigorous operations
# -----------------------------

class RigorousMath:
    """Facade that provides rigorous elementary functions and some special functions.
    Picks the best available backend.
    """

    def __init__(self):
        self.mode = "decimal"
        if HAVE_ARB:
            self.mode = "arb"
        elif HAVE_MPMATH:
            self.mode = "iv"

    # ---------------------
    # Precision control
    # ---------------------
    def set_digits(self, dps: int):
        """Set precision for the active backend (in decimal digits)."""
        if HAVE_ARB and self.mode == "arb":
            # python-flint uses bits; set via ctx.prec
            #  dps ~ bits * log10(2)  => bits ~ ceil(dps / log10(2))
            import math
            bits = max(64, int(math.ceil(dps / 0.3010299956639812)))
            ctx.prec = bits  # type: ignore
        elif HAVE_MPMATH and self.mode == "iv":
            import mpmath as mp  # type: ignore
            mp.mp.dps = max(50, dps)

    # ---------------------
    # Elementary functions
    # ---------------------
    def exp(self, iv: RInterval) -> RInterval:
        if self.mode == "arb":
            # Build an arb ball that covers [lo,hi] by mid+rad
            mid = (iv.lo + iv.hi) / Decimal(2)
            rad = (iv.hi - iv.lo) / Decimal(2)
            x = arb(str(mid))  # type: ignore
            if rad != 0:
                x = x.add_error(str(rad))  # type: ignore
            y = x.exp()  # type: ignore
            # Convert back to [lo,hi]
            lo = Decimal(str(y.lower()))  # type: ignore
            hi = Decimal(str(y.upper()))  # type: ignore
            return RInterval(lo, hi)
        elif self.mode == "iv":
            import mpmath as mp  # type: ignore
            xin = mp.iv.mpf([str(iv.lo), str(iv.hi)])  # type: ignore
            y = mp.iv.exp(xin)  # type: ignore
            return RInterval(Decimal(str(y.a)), Decimal(str(y.b)))
        else:
            # Decimal fallback with symmetric padding
            import math
            el = Decimal(str(math.exp(float(iv.lo))))
            eh = Decimal(str(math.exp(float(iv.hi))))
            pad = Decimal("1e-15") * (abs(el) + abs(eh) + Decimal(1))
            return RInterval(el - pad, eh + pad)

    def log(self, iv: RInterval) -> RInterval:
        if iv.lo <= 0:
            raise ValueError("log domain requires (0,∞)")
        if self.mode == "arb":
            mid = (iv.lo + iv.hi) / Decimal(2)
            rad = (iv.hi - iv.lo) / Decimal(2)
            x = arb(str(mid))  # type: ignore
            if rad != 0: x = x.add_error(str(rad))  # type: ignore
            y = x.log()  # type: ignore
            lo = Decimal(str(y.lower()))  # type: ignore
            hi = Decimal(str(y.upper()))  # type: ignore
            return RInterval(lo, hi)
        elif self.mode == "iv":
            import mpmath as mp  # type: ignore
            xin = mp.iv.mpf([str(iv.lo), str(iv.hi)])  # type: ignore
            y = mp.iv.log(xin)  # type: ignore
            return RInterval(Decimal(str(y.a)), Decimal(str(y.b)))
        else:
            import math
            ll = Decimal(str(math.log(float(iv.lo))))
            lh = Decimal(str(math.log(float(iv.hi))))
            pad = Decimal("1e-15") * (abs(ll) + abs(lh) + Decimal(1))
            return RInterval(ll - pad, lh + pad)

    def sqrt(self, iv: RInterval) -> RInterval:
        if iv.lo < 0:
            raise ValueError("sqrt domain requires [0,∞)")
        if self.mode == "arb":
            mid = (iv.lo + iv.hi) / Decimal(2)
            rad = (iv.hi - iv.lo) / Decimal(2)
            x = arb(str(mid))  # type: ignore
            if rad != 0: x = x.add_error(str(rad))  # type: ignore
            y = x.sqrt()  # type: ignore
            lo = Decimal(str(y.lower()))  # type: ignore
            hi = Decimal(str(y.upper()))  # type: ignore
            return RInterval(lo, hi)
        elif self.mode == "iv":
            import mpmath as mp  # type: ignore
            xin = mp.iv.mpf([str(iv.lo), str(iv.hi)])  # type: ignore
            y = mp.iv.sqrt(xin)  # type: ignore
            return RInterval(Decimal(str(y.a)), Decimal(str(y.b)))
        else:
            import math
            sl = Decimal(str(math.sqrt(float(iv.lo))))
            sh = Decimal(str(math.sqrt(float(iv.hi))))
            pad = Decimal("1e-15") * (abs(sl) + abs(sh) + Decimal(1))
            return RInterval(sl - pad, sh + pad)

    # ---------------------
    # Special functions (digamma)
    # ---------------------
    def digamma_real_part(self, sigma: Decimal, t: Decimal) -> RInterval:
        """Return an interval enclosing Re ψ(sigma + i t). Requires ACB (arb).
        If arb is unavailable, raises NotImplementedError.
        """
        if self.mode != "arb":
            raise NotImplementedError("digamma_real_part requires python-flint (arb/acb).")
        z = acb(str(sigma), str(t))  # type: ignore
        psi = z.digamma()            # type: ignore (acb.digamma)
        # acb.real returns arb
        re = psi.real()              # type: ignore
        lo = Decimal(str(re.lower()))  # type: ignore
        hi = Decimal(str(re.upper()))  # type: ignore
        return RInterval(lo, hi)


# Singleton instance
RMATH = RigorousMath()
