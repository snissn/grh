#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
explicit_formula_ms_working.py  (oscillatory-safe)

Maaß–Selberg explicit formula check for gb (cosine-profile) on [a,b].
Profiles:
  rect   — hard rectangle gb=1_[a,b] (closed forms; converges slowly in zeros)
  smooth — C^∞ bump-edged rectangle (uses oscillatory-aware quadrature)

Identity (MS ledger):
  2 * sum_{rho} H(g; rho) = P(g) + 1/2 g(0) + 2 J(g) - A(g)

Key fixes in this version:
  • H_smooth_for_rho uses mp.quadosc with period=2π/γ for oscillatory stability.
  • zeros_from_mpmath coerces zetazero(n) → 1/2 + i*γ_n robustly.
"""

from __future__ import annotations
import argparse, csv
from typing import List, Optional
from mpmath import mp, zetazero

# ---------- utilities ----------
def ensure_mp(dps: int):
    mp.dps = int(dps)

def sieve_primes(nmax: int) -> List[int]:
    if nmax < 2: return []
    sieve = bytearray(b"\x01") * (nmax + 1)
    sieve[0:2] = b"\x00\x00"
    p = 2
    while p*p <= nmax:
        if sieve[p]:
            start = p*p
            step = p
            sieve[start:nmax+1:step] = b"\x00" * (((nmax - start)//step) + 1)
        p += 1
    return [i for i in range(2, nmax+1) if sieve[i]]

# ---------- profiles (rect & smooth) ----------
def g0_rect(a: float, b: float): return (mp.mpf(b)-mp.mpf(a))/mp.pi
def J_rect(a: float, b: float):  return (mp.e**(-a/2) - mp.e**(-b/2))/mp.pi
def A_rect(a: float, b: float):
    term = (mp.mpf(b)-mp.mpf(a))/2 - mp.log1p(mp.e**(b/2)) + mp.log1p(mp.e**(a/2))
    return term/mp.pi
def P_rect(a: float, b: float) -> mp.mpf:
    total = mp.mpf('0')
    pmax = int(mp.floor(mp.e**(b/2))) + 5
    for p in sieve_primes(pmax):
        lp = mp.log(p)
        kmin = max(1, int(mp.ceil(a/(2*lp))))
        kmax = int(mp.floor(b/(2*lp)))
        for k in range(kmin, kmax+1):
            total += lp * (1 - p**(-k))
    return total/mp.pi

# Smooth C^∞ edge ramp
def _bump01(t: mp.mpf) -> mp.mpf:
    if t <= 0: return mp.mpf('0')
    if t >= 1: return mp.mpf('1')
    s = mp.e**(-1/t)
    s2 = mp.e**(-1/(1-t))
    return s/(s+s2)

def gb_smooth(x: mp.mpf, a: float, b: float, eps: float) -> mp.mpf:
    if x < a or x > b: return mp.mpf('0')
    if eps <= 0: return mp.mpf('1')
    if x <= a+eps: return _bump01((x-a)/eps)
    if x >= b-eps: return _bump01((b-x)/eps)
    return mp.mpf('1')

def g0_smooth(a: float, b: float, eps: float) -> mp.mpf:
    return (1/mp.pi)*mp.quad(lambda x: gb_smooth(x,a,b,eps), [a,b])

def J_smooth(a: float, b: float, eps: float) -> mp.mpf:
    return (1/(2*mp.pi))*mp.quad(lambda x: mp.e**(-x/2)*gb_smooth(x,a,b,eps), [a,b])

def T_int_smooth(a: float, b: float, eps: float) -> mp.mpf:
    return (1/(2*mp.pi))*mp.quad(lambda x: gb_smooth(x,a,b,eps)/(1+mp.e**(-x/2)), [a,b])

def A_from_T(g0: mp.mpf, T_int: mp.mpf) -> mp.mpf:
    return mp.mpf('0.5')*g0 - T_int

def P_smooth(a: float, b: float, eps: float) -> mp.mpf:
    total = mp.mpf('0')
    pmax = int(mp.floor(mp.e**(b/2))) + 100
    for p in sieve_primes(pmax):
        lp = mp.log(p)
        kmin = max(1, int(mp.ceil(a/(2*lp))))
        kmax = int(mp.floor(b/(2*lp)))
        for k in range(kmin, kmax+1):
            x = 2*k*lp
            total += lp*(1 - p**(-k))*gb_smooth(x,a,b,eps)
    return total/mp.pi

# ---------- zero-side H ----------
def H_rect_for_rho(a: float, b: float, rho: complex) -> mp.mpf:
    beta = mp.mpf(mp.re(rho))
    gamma = mp.mpf(mp.im(rho))
    L = mp.mpf(a)/2; U = mp.mpf(b)/2
    den = beta*beta + gamma*gamma
    def F(x):  # CORRECT antiderivative
        return mp.e**(-beta*x)*(-beta*mp.cos(gamma*x)+gamma*mp.sin(gamma*x))/den
    return (F(U)-F(L))/mp.pi

def H_smooth_for_rho(a: float, b: float, eps: float, rho: complex) -> mp.mpf:
    beta = mp.mpf(mp.re(rho))
    gamma = mp.mpf(mp.im(rho))
    if gamma == 0:
        # fallback to non-oscillatory quad (rare)
        f = lambda x: mp.e**(-beta*x)*gb_smooth(2*x,a,b,eps)
        return (1/mp.pi)*mp.quad(f, [a/2,b/2])
    # Oscillatory-aware quadrature with known period
    period = 2*mp.pi/gamma
    f = lambda x: mp.e**(-beta*x)*gb_smooth(2*x,a,b,eps)*mp.cos(gamma*x)
    try:
        val = mp.quadosc(f, [a/2, b/2], period=period, maxn=100000)
    except Exception:
        # fallback: split into many subintervals of length ≈ period
        L = a/2; U = b/2
        nseg = max(50, int(mp.ceil((U-L)/(period/2))))
        xs = [L + i*(U-L)/nseg for i in range(nseg)] + [U]
        val = mp.mpf('0')
        for i in range(nseg):
            val += mp.quad(f, [xs[i], xs[i+1]])
    return (1/mp.pi)*val

# ---------- zeros ----------
def zeros_from_mpmath(N: int) -> List[complex]:
    zs = []
    for n in range(1, N+1):
        z = zetazero(n)
        # robust coercion: some mpmath versions return γ (real), some 1/2+iγ
        try:
            im = mp.im(z); re = mp.re(z)
        except Exception:
            im = mp.mpf(z); re = mp.mpf('0.5')
        else:
            if im == 0:  # real
                im = mp.mpf(re); re = mp.mpf('0.5')
        zs.append(re + 1j*im)
    return zs

def zeros_from_csv(csv_path: str, N: Optional[int]=None) -> List[complex]:
    gammas = []
    with open(csv_path,'r',newline='') as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            gammas.append(mp.mpf(row['gamma'].strip()))
    if N is not None:
        gammas = gammas[:N]
    return [mp.mpf('0.5') + 1j*g for g in gammas]

# ---------- partial sums ----------
def partial_zero_sum(H_vals: List[mp.mpf], gammas: List[mp.mpf], fejer: bool) -> List[mp.mpf]:
    vals = []
    S = mp.mpf('0')
    if not H_vals: return vals
    T = gammas[-1]
    for n, hv in enumerate(H_vals, start=1):
        w = (1 - gammas[n-1]/T) if fejer else 1
        S += w*hv
        vals.append(2*S)  # conjugate pair
    return vals

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", choices=["rect","smooth"], default="rect")
    ap.add_argument("--a", type=float, default=1.0)
    ap.add_argument("--b", type=float, default=2.0)
    ap.add_argument("--eps", type=float, default=0.15)
    ap.add_argument("--N", type=int, default=200)
    ap.add_argument("--dps", type=int, default=80)
    ap.add_argument("--use-mpmath", action="store_true")
    ap.add_argument("--zeros-csv", type=str, default=None)
    ap.add_argument("--fejer", action="store_true")
    ap.add_argument("--plot-file", type=str, default=None)
    args = ap.parse_args()

    ensure_mp(args.dps)
    a,b = float(args.a), float(args.b)
    if b<=a: raise SystemExit("Require b>a.")

    print("=== Maaß–Selberg explicit formula (MS ledger) ===")
    print(f"profile={args.profile}, a={a}, b={b}" + (f", eps={args.eps}" if args.profile=='smooth' else ""))
    print(f"precision (mp.dps) = {mp.dps}")

    # RHS blocks
    if args.profile=="rect":
        g0, J, A, P = g0_rect(a,b), J_rect(a,b), A_rect(a,b), P_rect(a,b)
    else:
        g0 = g0_smooth(a,b,args.eps)
        T_int = T_int_smooth(a,b,args.eps)
        J = J_smooth(a,b,args.eps)
        A = A_from_T(g0, T_int)
        P = P_smooth(a,b,args.eps)
    RHS = P + mp.mpf('0.5')*g0 + 2*J - A

    print(f"g(0)         = {g0}")
    print(f"J(g)         = {J}")
    print(f"A(g)         = {A}")
    print(f"P(g)         = {P}")
    print(f"RHS (P + 1/2 g0 + 2J - A) = {RHS}")

    # zeros
    if args.use_mpmath:
        zeros = zeros_from_mpmath(args.N)
    elif args.zeros_csv:
        zeros = zeros_from_csv(args.zeros_csv, N=args.N)
    else:
        print("No zeros source; use --use-mpmath or --zeros-csv.")
        return
    gammas = [mp.im(z) for z in zeros]
    print(f"Loaded {len(zeros)} zeros. First γ: {gammas[:5]} ... γ_N = {gammas[-1]}")

    # LHS
    if args.profile=="rect":
        H_vals = [H_rect_for_rho(a,b,z) for z in zeros]
    else:
        H_vals = [H_smooth_for_rho(a,b,args.eps,z) for z in zeros]
    LHS_vals = partial_zero_sum(H_vals, gammas, fejer=args.fejer)
    if not LHS_vals:
        print("No H values computed.")
        return
    LHS_N = LHS_vals[-1]
    print(f"LHS with N={args.N} zeros{' (Fejér)' if args.fejer else ''}: {LHS_N}")
    print(f"RHS - LHS(N) = {RHS - LHS_N}")

    # plot
    if args.plot_file:
        try:
            import matplotlib.pyplot as plt
            resids = [RHS - s for s in LHS_vals]
            ns = list(range(1,len(resids)+1))
            plt.figure()
            plt.plot(ns, resids)
            plt.xlabel("Number of zeros N")
            plt.ylabel("RHS - 2 Σ_{n≤N} H(g; ρ_n)")
            plt.title(f"Residual (profile={args.profile}, a={a}, b={b}" + (f", eps={args.eps}" if args.profile=='smooth' else "") + ")")
            plt.tight_layout()
            plt.savefig(args.plot_file, dpi=150)
            print(f"Saved plot to {args.plot_file}")
        except Exception as e:
            print(f"Plotting failed: {e}")

if __name__ == "__main__":
    main()


