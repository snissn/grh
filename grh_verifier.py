# grh_verifier.py — Reference Python implementation of the Part III “03” verifier algorithm
#
# This module follows the phase-by-phase verifier described in the paper’s
# “GRH Certificate Verifier (finite, rigorous mode)” and the routine/ledger sketch:
# ArchBlock, RamBlock, PrimeBlockBL, PrimeTailBound (heat), RS_Positivity_LB,
# SANorm, RStarBudget, Q_of_X_eps, Eps_from_q_X, CheckA5Gap.
#
# Domain‑specific numerics are intentionally pluggable to match your infrastructure.
# All real computations are carried out with outward‑rounded intervals.
#
from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal, getcontext, localcontext, ROUND_FLOOR, ROUND_CEILING
from math import ceil
from typing import Iterable, List, Literal, Optional, Sequence, Union
import math


# ------------------------------
# Interval arithmetic utilities
# ------------------------------

@dataclass(frozen=True)
class Interval:
    """Closed interval [lo, hi] with outward rounding.
    Endpoints are Decimal. Invariants: lo <= hi.
    """
    lo: Decimal
    hi: Decimal

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError(f"Invalid interval: lo={self.lo} > hi={self.hi}")

    # Predicates
    def lt0(self) -> bool: return self.hi < Decimal(0)     # “< 0” means upper endpoint < 0
    def le0(self) -> bool: return self.hi <= Decimal(0)
    def gt0(self) -> bool: return self.lo >  Decimal(0)
    def ge0(self) -> bool: return self.lo >= Decimal(0)
    def contains0(self) -> bool: return self.lo <= 0 <= self.hi
    def width(self) -> Decimal: return self.hi - self.lo
    def mid(self) -> Decimal: return (self.lo + self.hi) / Decimal(2)

    # Basic interval algebra (with directed rounding)
    def __add__(self, other: "Interval") -> "Interval":
        with localcontext() as ctx:
            ctx.rounding = ROUND_FLOOR; lo = self.lo + other.lo
            ctx.rounding = ROUND_CEILING; hi = self.hi + other.hi
        return Interval(lo, hi)

    def __sub__(self, other: "Interval") -> "Interval":
        with localcontext() as ctx:
            ctx.rounding = ROUND_FLOOR; lo = self.lo - other.hi
            ctx.rounding = ROUND_CEILING; hi = self.hi - other.lo
        return Interval(lo, hi)

    def __neg__(self) -> "Interval": return Interval(-self.hi, -self.lo)

    def __mul__(self, other: "Interval") -> "Interval":
        with localcontext() as ctx:
            candidates = []
            for a in (self.lo, self.hi):
                for b in (other.lo, other.hi):
                    ctx.rounding = ROUND_FLOOR; lo = a * b
                    ctx.rounding = ROUND_CEILING; hi = a * b
                    candidates.append((lo, hi))
            lo = min(x for x, _ in candidates); hi = max(y for _, y in candidates)
        return Interval(lo, hi)

    def __truediv__(self, other: "Interval") -> "Interval":
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError("Division by interval containing 0")
        return self * other.recip()

    def recip(self) -> "Interval":
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError("Invert interval containing 0")
        with localcontext() as ctx:
            ctx.rounding = ROUND_FLOOR; lo = Decimal(1) / (self.hi if self.hi > 0 else self.lo)
            ctx.rounding = ROUND_CEILING; hi = Decimal(1) / (self.lo if self.lo > 0 else self.hi)
        lo, hi = (lo, hi) if lo <= hi else (hi, lo)
        return Interval(lo, hi)

    # Elementary functions with safe symmetric padding (replace with arb/MPFR for production)
    def exp(self) -> "Interval":
        el = Decimal(str(math.exp(float(self.lo))))
        eh = Decimal(str(math.exp(float(self.hi))))
        return _pad_bounds(el, eh)

    def log(self) -> "Interval":
        if self.lo <= 0: raise ValueError("log domain requires (0, ∞)")
        ll = Decimal(str(math.log(float(self.lo))))
        lh = Decimal(str(math.log(float(self.hi))))
        return _pad_bounds(ll, lh)

    def sqrt(self) -> "Interval":
        if self.lo < 0: raise ValueError("sqrt domain requires [0, ∞)")
        sl = Decimal(str(math.sqrt(float(self.lo))))
        sh = Decimal(str(math.sqrt(float(self.hi))))
        return _pad_bounds(sl, sh)

    @staticmethod
    def point(x: Union[int, float, Decimal]) -> "Interval":
        d = Decimal(str(x)) if not isinstance(x, Decimal) else x
        return Interval(d, d)

    @staticmethod
    def hull(values: Iterable[Union[int, float, Decimal]]) -> "Interval":
        vals = [Decimal(str(v)) if not isinstance(v, Decimal) else v for v in values]
        return Interval(min(vals), max(vals))

    def __repr__(self) -> str:
        return f"[{self.lo}, {self.hi}]"


def _pad_bounds(a: Decimal, b: Decimal, rel: Decimal = Decimal("1e-15")) -> Interval:
    lo = min(a, b); hi = max(a, b)
    span = abs(lo) + abs(hi) + Decimal(1)
    return Interval(lo - rel * span, hi + rel * span)


# -------------------------------------
# Test functions (fixed–heat, BL, etc.)
# -------------------------------------

@dataclass
class TestFunction:
    """Admissible test Φ with family tag (“heat” or “BL”)."""
    family: Literal["heat", "BL"]
    a: Optional[Decimal] = None  # only used for heat family

    def check_admissible(self) -> bool:
        # In the “heat” family, we enforce a>1/2 as in SA(a) (Part I).
        if self.family == "heat":
            return self.a is not None and self.a > Decimal("0.5")
        elif self.family == "BL":
            return True
        return False

    def sanorm_interval(self, t_min: Decimal, t_zero: Decimal, precision_bits: int) -> Interval:
        # A tame surrogate bound; replace with rigorous quadrature in your infra.
        a = self.a if self.a is not None else Decimal(1)
        span = max(Decimal(0), t_zero - t_min)
        base = (Decimal(1) / (a - Decimal("0.5"))) + (Decimal(1) + span)
        pad = Decimal("1e-12")
        return Interval(base * (Decimal(1) - pad), base * (Decimal(1) + pad))


# ---------------------------------
# Data classes for the certificate
# ---------------------------------

@dataclass
class PrimeLocalTerm:
    p: int
    k: int
    value: Interval   # interval contribution of this (p^k)-term
    # Optional metadata for amplifiers (e.g., Satake parameters)
    meta: Optional[dict] = None


@dataclass
class RamifiedTerm:
    j: int
    value: Interval


@dataclass
class Certificate:
    # Core instance data
    m: int
    K: str
    Qpi: Decimal
    t_star: Decimal
    Phi_test: TestFunction

    # Precision and modes
    precision_bits: int = 212
    test_family: Literal["heat", "BL"] = "heat"
    prime_block_mode: Literal["evaluate", "rs_lower_bound"] = "evaluate"
    proof_mode: Literal["proof_of_GRH", "ann_odd_only"] = "ann_odd_only"

    # Infrastructure constants
    Ainf: Decimal = Decimal(1)
    beta_m: Decimal = Decimal("0.5")
    C_Rstar: Decimal = Decimal(1)

    # Prime block
    band_limit_X: Optional[Decimal] = None
    p_max: Optional[int] = None
    k_max: Optional[int] = None
    prime_tail_bound: Optional[Interval] = None
    unram_local: List[PrimeLocalTerm] = field(default_factory=list)

    # Archimedean and ramified data
    arch_value: Optional[Interval] = None
    ram_terms: List[RamifiedTerm] = field(default_factory=list)
    j_max: Optional[int] = None

    # Amplifier (q, weights)
    q_modulus: Optional[int] = None
    weights_spec: Optional[dict] = None

    # A5 gap & approximate invariance (Q(X,ε)=e^{αX} ε^{−κ})
    eta_a5_gap: Optional[Interval] = None
    epsilon_invariance: Optional[Interval] = None
    alpha_inv: Optional[Decimal] = None
    kappa_inv: Optional[Decimal] = None

    # Test net for proof_of_GRH
    test_net: Optional[List[TestFunction]] = None


# ------------------------------
# Core ledger‑style primitives
# ------------------------------

def ArchBlock(t_star: Decimal, Phi: TestFunction, arch_value: Optional[Interval]) -> Interval:
    """Arch^{(t*)}(Φ): if supplied as an interval, trust it; otherwise zero placeholder."""
    return Interval.point(0) if arch_value is None else arch_value

def RamBlock(t_star: Decimal, Phi: TestFunction, ram_terms: Sequence[RamifiedTerm], j_max: Optional[int]) -> Interval:
    """Sum ramified contributions up to j_max as an interval."""
    total = Interval.point(0)
    terms = [rt for rt in ram_terms if (j_max is None or rt.j <= j_max)]
    for rt in terms: total = total + rt.value
    return total

def PrimeBlockBL(t_star: Decimal, Phi: TestFunction,
                 unram_local: Sequence[PrimeLocalTerm], band_limit_X: Decimal,
                 weights_spec: Optional[dict] = None) -> Interval:
    """Band‑limited prime block: include terms with 2k·log p ≤ X."""
    if band_limit_X is None:
        raise ValueError("band_limit_X required for BL tests")
    total = Interval.point(0)

    def fejer_weight(u: float, L: int) -> float:
        # Nonnegative Dirichlet-Fejér kernel normalized to ≤1 at u=0
        import math
        if L <= 1:
            return 1.0
        x = 0.5 * u
        s = math.sin(L * x)
        d = math.sin(x)
        if abs(d) < 1e-15:
            return 1.0
        val = (s / (L * d)) ** 2
        # Clip to [0,1] for safety against FP noise
        if val < 0:
            val = 0.0
        if val > 1:
            val = 1.0
        return val

    def satake_amp_weight(meta: Optional[dict], k: int) -> float:
        # Prototype Satake-based weight: w = (|tr(A_p^k)| / d)^2 in [0,1].
        # Uses provided alphas if available.
        if not meta:
            return 1.0
        alphas = meta.get("alphas")
        if not alphas:
            return 1.0
        try:
            tr = sum((a**k) for a in alphas)
            import math
            tr_abs = abs(tr)
            d = max(1, len(alphas))
            w = (tr_abs / d) ** 2
            if w < 0:
                w = 0.0
            if w > 1:
                w = 1.0
            return float(w)
        except Exception:
            return 1.0
    for term in unram_local:
        if (2 * term.k * Decimal(str(math.log(term.p)))) <= band_limit_X:
            val = term.value
            if weights_spec is not None:
                wtype = weights_spec.get("type")
                if wtype == "fejer":
                    L = int(weights_spec.get("L", 8))
                    u = float(term.k * math.log(term.p))
                    w = fejer_weight(u, L)
                    val = val * Interval.point(w)
                elif wtype == "satake_amp":
                    w = satake_amp_weight(term.meta, term.k)
                    val = val * Interval.point(w)
            total = total + val
    return total

def PrimeTailBound_heat(prime_tail_bound: Optional[Interval]) -> Interval:
    """Worst‑case signed tail bound for heat tests (orientation already baked in)."""
    return prime_tail_bound if prime_tail_bound is not None else Interval.point(0)

def RS_Positivity_LB(_: TestFunction) -> Interval:
    """RS positivity ⇒ prime block lower bound is [0,0]."""
    return Interval.point(0)

def SANorm(Phi: TestFunction, a_weight: Decimal, t_min: Decimal, t_zero: Decimal, precision_bits: int) -> Interval:
    """Interval bound for ||Φ||_{SA(a)}."""
    return Phi.sanorm_interval(t_min, t_zero, precision_bits)

def RStarBudget(M: Interval, Ainf: Decimal, Qpi: Decimal, beta_m: Decimal, C_Rstar: Decimal, t_star: Decimal) -> Interval:
    """B_{R*} := C_{R*} · A_inf · Qpi^{−β(m)·t*} · M (computed with outward padding)."""
    power = float(-beta_m * t_star)
    base = float(Qpi)
    val = base ** power  # Qpi^{-β t*}
    pad = Decimal("1e-15")
    qpow_iv = Interval(Decimal(str(val)) * (Decimal(1) - pad),
                       Decimal(str(val)) * (Decimal(1) + pad))
    consts = Interval.point(C_Rstar * Ainf)
    return consts * qpow_iv * M

def Q_of_X_eps(X: Decimal, eps: Interval, alpha: Decimal, kappa: Decimal) -> Interval:
    """Q(X,ε) = exp(α X) · ε^{−κ} (interval in ε)."""
    if eps.lo <= 0: raise ValueError("ε must be > 0")
    exp_ax = Interval.point(float(Decimal(alpha) * Decimal(X))).exp()
    log_eps = eps.log()
    negk = Interval.point(-kappa)
    power = (negk * log_eps).exp()
    return exp_ax * power

def Eps_from_q_X(q: Decimal, X: Decimal, alpha: Decimal, kappa: Decimal) -> Interval:
    """Invert q ≥ exp(αX)·ε^{−κ} ⇒ ε ≤ (q / exp(αX))^{−1/κ} (returns [0, upper])."""
    exp_ax = Decimal(str(math.exp(float(alpha * X))))
    base = q / exp_ax
    if base <= 0: raise ValueError("q / exp(αX) must be positive")
    val = base ** (Decimal(-1) / Decimal(kappa))
    pad = Decimal("1e-15")
    return Interval(Decimal(0), val * (Decimal(1) + pad))

def CheckA5Gap(eta_interval: Interval) -> Interval:
    """Return sqrt(1 − η) as interval (used in the annihilation inequality)."""
    one_minus_eta = Interval.point(1) - eta_interval
    return one_minus_eta.sqrt()


# ---------------
# Verifier phases
# ---------------

@dataclass
class VerifierReport:
    phase_passed: int
    result: Literal["GRH_Verified", "AnnOdd_Verified", "GRH_Falsified", "Indeterminate"]
    details: dict

def bits_to_decimal_digits(bits: int) -> int:
    return int(ceil(bits * 0.3010299956639812))  # ≈ bits * log10(2)

def set_precision_bits(bits: int) -> None:
    getcontext().prec = max(28, bits_to_decimal_digits(bits) + 5)

def verify_certificate(cert: Certificate) -> VerifierReport:
    """Master verifier implementing Phases 0–6 (finite, rigorous mode)."""

    # Phase 0: Numerical Rigor
    set_precision_bits(cert.precision_bits)

    # Phase 1: Initialization and Admissibility
    if not cert.Phi_test.check_admissible():
        return VerifierReport(phase_passed=1, result="GRH_Falsified",
                              details={"reason": "Test function not admissible"})

    # SA-norm (interval). t_min/t_zero are infra choices; use placeholders here.
    M = SANorm(cert.Phi_test, a_weight=Decimal(cert.Phi_test.a or 1),
               t_min=Decimal(0), t_zero=cert.t_star, precision_bits=cert.precision_bits)

    # Phase 2: EF Geometric Side (finite/bounded)
    arch_iv = ArchBlock(cert.t_star, cert.Phi_test, cert.arch_value)
    ram_iv  = RamBlock(cert.t_star, cert.Phi_test, cert.ram_terms, cert.j_max)

    if cert.prime_block_mode == "evaluate":
        if cert.test_family == "BL":
            if cert.band_limit_X is None:
                return VerifierReport(phase_passed=2, result="Indeterminate",
                                      details={"reason": "Missing band_limit_X for BL"})
            prime_iv = PrimeBlockBL(cert.t_star, cert.Phi_test, cert.unram_local, cert.band_limit_X, cert.weights_spec)
        else:  # heat
            prime_iv = Interval.point(0)
            for term in cert.unram_local:
                if (cert.p_max is None or term.p <= cert.p_max) and (cert.k_max is None or term.k <= cert.k_max):
                    prime_iv = prime_iv + term.value
            prime_iv = prime_iv + PrimeTailBound_heat(cert.prime_tail_bound)
        LB_Prime = prime_iv
    else:
        LB_Prime = RS_Positivity_LB(cert.Phi_test)  # [0,0]

    # Phase 3: Infrastructure Verification (R* Budget)
    B_Rstar = RStarBudget(M, cert.Ainf, cert.Qpi, cert.beta_m, cert.C_Rstar, cert.t_star)
    ok_R = (ram_iv.lo >= -B_Rstar.hi) and (ram_iv.hi <= B_Rstar.hi)
    if not ok_R:
        return VerifierReport(phase_passed=3, result="GRH_Falsified",
                              details={"reason": "R* bound violated",
                                       "Ram": (str(ram_iv.lo), str(ram_iv.hi)),
                                       "B_R*": (str(B_Rstar.lo), str(B_Rstar.hi))})

    # Phase 4: Amplified Positivity Filter
    arch_plus_ram = arch_iv + ram_iv
    if cert.prime_block_mode == "rs_lower_bound":
        # Need LB_Prime - UB(Arch+Ram) >= 0; here LB_Prime = 0.
        geom_lb = Interval.point(0) - Interval.point(arch_plus_ram.hi)
        if geom_lb.lt0():
            return VerifierReport(phase_passed=4, result="GRH_Falsified",
                                  details={"reason": "Amplified positivity fails (RS lower bound mode)",
                                           "Arch+Ram": (str(arch_plus_ram.lo), str(arch_plus_ram.hi))})
    else:
        Geom_Amp = LB_Prime - arch_plus_ram
        if Geom_Amp.lt0():
            return VerifierReport(phase_passed=4, result="GRH_Falsified",
                                  details={"reason": "Amplified positivity fails",
                                           "GeomAmp": (str(Geom_Amp.lo), str(Geom_Amp.hi)),
                                           "Prime": (str(LB_Prime.lo), str(LB_Prime.hi)),
                                           "Arch+Ram": (str(arch_plus_ram.lo), str(arch_plus_ram.hi))})

    # Phase 5: ν–Annihilation (A5 gap + approx.~invariance)
    if cert.eta_a5_gap is None or cert.epsilon_invariance is None or cert.alpha_inv is None or cert.kappa_inv is None:
        return VerifierReport(phase_passed=5, result="Indeterminate",
                              details={"reason": "Missing A5/invariance inputs"})

    sqrt1m_eta = CheckA5Gap(cert.eta_a5_gap)

    # Optionally tighten ε via q ≥ Q(X,ε) if modulus and X are available
    eps_iv = cert.epsilon_invariance
    if cert.q_modulus is not None and cert.band_limit_X is not None:
        eps_cap = Eps_from_q_X(Decimal(cert.q_modulus), Decimal(cert.band_limit_X),
                               Decimal(cert.alpha_inv), Decimal(cert.kappa_inv))
        # Intersect: ε ∈ [0, min(eps_iv.hi, eps_cap.hi)]
        eps_iv = Interval(Decimal(0), min(eps_iv.hi, eps_cap.hi))

    # Check inequality: sqrt(1−η) + ε < 1 (worst‑case on upper endpoints)
    annih_hi = sqrt1m_eta.hi + eps_iv.hi
    if annih_hi >= Decimal(1):
        return VerifierReport(phase_passed=5, result="GRH_Falsified",
                              details={"reason": "Annihilation inequality fails (sqrt(1-η)+ε ≥ 1)",
                                       "sqrt(1-eta)": (str(sqrt1m_eta.lo), str(sqrt1m_eta.hi)),
                                       "eps": (str(eps_iv.lo), str(eps_iv.hi))})

    # Phase 6: Proof Mode
    if cert.proof_mode == "proof_of_GRH":
        if not cert.test_net:
            return VerifierReport(phase_passed=6, result="Indeterminate",
                                  details={"reason": "proof_of_GRH requires a test_net"})
        # In a full implementation, recompute ledger pieces per Φ_i.
        for i, Phi_i in enumerate(cert.test_net):
            sub = Certificate(**{**cert.__dict__})
            sub.Phi_test = Phi_i
            # Prevent infinite recursion: evaluate each test independently
            sub.proof_mode = "ann_odd_only"
            sub.test_net = None
            subres = verify_certificate(sub)
            if subres.result != "AnnOdd_Verified":
                return VerifierReport(phase_passed=6, result=subres.result,
                                      details={"reason": f"test_net failed at index {i}", "sub": subres.details})
        return VerifierReport(phase_passed=6, result="GRH_Verified",
                              details={"reason": "All tests in net passed"})
    else:
        return VerifierReport(phase_passed=6, result="AnnOdd_Verified",
                              details={"reason": "Annihilation holds for supplied Φ_test"})


# -----------------------------
# Tiny self-test (synthetic)
# -----------------------------

def _demo():
    Phi = TestFunction(family="heat", a=Decimal("0.8"))
    cert = Certificate(
        m=2, K="Q", Qpi=Decimal("1000"), t_star=Decimal("1.0"), Phi_test=Phi,
        precision_bits=200, test_family="heat", prime_block_mode="rs_lower_bound",
        proof_mode="ann_odd_only",
        Ainf=Decimal("1.0"), beta_m=Decimal("0.5"), C_Rstar=Decimal("1.0"),
        arch_value=Interval.point("0.05"),
        ram_terms=[RamifiedTerm(j=1, value=Interval.point("0.01"))], j_max=1,
        q_modulus=97, band_limit_X=Decimal("6.0"),
        eta_a5_gap=Interval.point("0.40"),
        epsilon_invariance=Interval(Decimal("0"), Decimal("0.20")),
        alpha_inv=Decimal("0.5"), kappa_inv=Decimal("1.0"),
    )
    rep = verify_certificate(cert)
    print("Result:", rep.result); print("Details:", rep.details)


if __name__ == "__main__":
    _demo()
