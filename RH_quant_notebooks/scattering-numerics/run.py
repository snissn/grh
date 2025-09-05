#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Camera‑ready demo: Zeta zeros from phase — two independent quantizers

METHOD A (standard, “bound‑state ruler”):
    Q_RS(τ) = ( θ(2τ) + arg ζ(1/2 + i·2τ) ) / π
    Zeros are read off from half‑integers:  Q_RS(τ_n) = n + 1/2  ⇒  γ_n = 2 τ_n
    This is the classical Riemann–Siegel (RS) quantizer. It uses only analytic functions.

METHOD B (novel framing, “scattering ruler”):
    Q_raw(τ) = − ( φ(τ) − φ(τ0) ) / π,  where  φ(τ) = arg S_e(1/2 + iτ)
        (Eisenstein scattering phase from the GL₂ cusp channel)
    Learn a tiny, smooth baseline L(τ) by fitting the analytic mismatch
        D(τ) = S(τ) − Q_raw(τ),  where  S(τ) = Q_RS(τ)  (the analytic RS ruler)
    Define Q_adj(τ) = Q_raw(τ) + L(τ) and quantize:  Q_adj(τ_n) = n + 1/2  ⇒  γ_n = 2 τ_n

Guarantees:
* No zeros are used by either quantizer; they are only used AFTER the fact for evaluation.
* Baseline L(τ) is smooth (poly5 or blockwise cubic); it cannot inject steps/zeros.

Outputs:
* Console analytics (mean/RMS/max errors vs true zeros)
* CSV tables in the current directory:
    - rs_quantizer_table.csv
    - rslock_eisenstein_table.csv
* Plots (four figures): RS staircase & residuals, baseline D vs L, RS‑Lock staircase & residuals

Requirements:  pip install mpmath numpy pandas matplotlib
"""

import math
import numpy as np
import pandas as pd
import mpmath as mp
import matplotlib.pyplot as plt

# ------------------------ Numeric precision ------------------------
mp.mp.dps = 70  # working precision
rtpi = mp.sqrt(mp.pi)

# ------------------------ Helper utilities ------------------------
def unwrap(x):
    """Phase unwrapping (numpy.unwrap) with float array output."""
    return np.unwrap(np.asarray(x, dtype=float))

def true_zeros_up_to(Tmax, Kmax=1000):
    """
    Return true zeta zeros γ_n (imag parts) up to Tmax using mpmath.zetazero.
    Used ONLY for evaluation, never by the quantizers.
    """
    out = []
    k = 1
    while len(out) < Kmax:
        z = mp.zetazero(k)
        g = float(mp.im(z))
        if g <= Tmax + 1e-12:
            out.append(g)
            k += 1
        else:
            break
    return out

def make_targets(Q_values):
    """Compute the half‑integer target list covered by the array Q(τ)."""
    Qmin = float(np.min(Q_values))
    Qmax = float(np.max(Q_values))
    targets = []
    k = math.ceil(Qmin - 0.5)
    while k + 0.5 <= Qmax:
        targets.append(k + 0.5)
        k += 1
    return targets

def crossing_brackets(taus, Q_values, target):
    """Indices i where [taus[i], taus[i+1]] brackets a crossing Q=target."""
    Q_values = np.asarray(Q_values, dtype=float)
    d = Q_values - target
    idxs = np.where(np.diff(np.sign(d)) != 0)[0]
    return idxs.tolist()

def refine_bisection(f_eval, a, b, target, iters=28):
    """
    Bisection for f(x)=target on [a,b].
    If numeric issues break the bracket, fall back to midpoint.
    """
    fa = f_eval(a) - target
    fb = f_eval(b) - target
    if not (np.isfinite(fa) and np.isfinite(fb) and fa * fb <= 0):
        return a + (b - a) * 0.5
    aa, bb = a, b
    for _ in range(iters):
        m = 0.5 * (aa + bb)
        fm = f_eval(m) - target
        if fa * fm <= 0:
            bb, fb = m, fm
        else:
            aa, fa = m, fm
    return 0.5 * (aa + bb)

def stats_vs_true(pred_gamma, true_gamma):
    m = min(len(pred_gamma), len(true_gamma))
    if m == 0:
        return 0, [], [], [], None, None, None
    pg = pred_gamma[:m]
    tg = true_gamma[:m]
    err = [p - t for p, t in zip(pg, tg)]
    arr = np.array(err, dtype=float)
    return m, pg, tg, err, float(arr.mean()), float(np.sqrt(np.mean(arr**2))), float(np.max(np.abs(arr)))


# ------------------------ Classical RS blocks ------------------------
def theta_RS(T):
    """Riemann–Siegel theta: θ(T) = Im log Γ(1/4 + iT/2) − (T/2) log π."""
    z = mp.mpf('0.25') + 0.5j * mp.mpf(T)
    return float(mp.im(mp.loggamma(z)) - (mp.mpf(T) / 2.0) * mp.log(mp.pi))

def arg_zeta_half(T):
    """Principal arg ζ(1/2 + iT) as a continuous function (will be branch‑aligned in bisection)."""
    return float(mp.arg(mp.zeta(0.5 + 1j * mp.mpf(T))))


# ------------------------ Eisenstein scattering phase ------------------------
def Se_tau(tau):
    """
    Eisenstein scattering multiplier (up to a fixed scalar):
        S_e(1/2 + iτ) ∝ Γ(iτ)/Γ(1/2 + iτ) · ζ(2iτ)/ζ(1 + 2iτ)
    """
    t = mp.mpf(tau)
    num = mp.gamma(1j * t) * mp.zeta(2j * t)
    den = mp.gamma(mp.mpf('0.5') + 1j * t) * mp.zeta(mp.mpf('1.0') + 2j * t)
    return rtpi * num / den

def phi_tau(tau):
    """Eisenstein phase φ(τ) = arg S_e(1/2 + iτ) (principal argument)."""
    return float(mp.arg(Se_tau(tau)))


# ------------------------ Baseline fitters (smooth, zero‑free) ------------------------
def fit_baseline_poly5(taus, D):
    """Global degree‑5 polynomial fit to D(τ) = S(τ) − Q_raw(τ)."""
    coef = np.polyfit(taus, D, 5)
    def L_eval(x):
        return float(np.polyval(coef, x))
    L_fit = np.polyval(coef, taus)
    info = ('poly5', coef)
    return L_fit, L_eval, info

def fit_baseline_block3(taus, D, blocks=6):
    """
    Piecewise degree‑3 polynomials over equal τ‑blocks (C^0 continuous).
    This reduces endpoint bias and adapts to mild curvature.
    """
    n = len(taus)
    edges = np.linspace(0, n - 1, blocks + 1).astype(int)
    coefs = []
    L_fit = np.zeros_like(D, dtype=float)
    for b in range(blocks):
        i0, i1 = edges[b], edges[b + 1]
        if i1 - i0 < 4:
            i1 = min(i0 + 4, n - 1)
        coef = np.polyfit(taus[i0:i1], D[i0:i1], 3)
        coefs.append((i0, i1, coef))
        L_fit[i0:i1] = np.polyval(coef, taus[i0:i1])
        if b == blocks - 1 and i1 < n:
            L_fit[i1:] = np.polyval(coef, taus[i1:])
    def L_eval(x):
        t = float(x)
        for (i0, i1, coef) in coefs:
            lo, hi = taus[i0], taus[min(i1 - 1, len(taus) - 1)]
            if lo <= t <= hi:
                return float(np.polyval(coef, t))
        return float(np.polyval(coefs[-1][2], t))
    info = ('block3', coefs)
    return L_fit, L_eval, info


# ------------------------ Quantizers ------------------------
def quantize_RS(taus):
    """
    Method A: classical RS quantizer.
    Returns:
        pred_gamma_RS, Q_RS_array, F_RS_array (unwrapped radians), targets
    """
    T = 2.0 * np.asarray(taus, dtype=float)
    theta_vals = unwrap([theta_RS(t) for t in T])
    argz_vals  = unwrap([arg_zeta_half(t) for t in T])
    F_RS_arr = theta_vals + argz_vals      # radians
    Q_RS_arr = F_RS_arr / math.pi          # in π‑units

    targets = make_targets(Q_RS_arr)

    # Bisection with branch alignment anchored to left endpoint of each bracket
    pred_tau = []
    for target in targets:
        idxs = crossing_brackets(taus, Q_RS_arr, target)
        for i in idxs:
            a, b = float(taus[i]), float(taus[i + 1])
            F_anchor = float(F_RS_arr[i])  # radians at left anchor
            def Q_eval(x):
                Tx = 2.0 * float(x)
                F = theta_RS(Tx) + arg_zeta_half(Tx)  # principal
                k2 = round((F_anchor - F) / (2 * math.pi))
                F_aligned = F + 2 * math.pi * k2
                return F_aligned / math.pi
            pred_tau.append(refine_bisection(Q_eval, a, b, target, iters=26))
    pred_tau = sorted(set([round(t, 12) for t in pred_tau]))
    pred_gamma = [2.0 * t for t in pred_tau]
    return pred_gamma, Q_RS_arr, F_RS_arr, targets

def quantize_RSLock(taus, Q_RS_arr, baseline_method='poly5', blocks=6):
    """
    Method B: RS‑Lock Eisenstein quantizer (scattering phase + smooth baseline)
    Inputs:
        taus: grid
        Q_RS_arr: the analytic RS ruler S(τ) (already computed by quantize_RS)
        baseline_method: 'poly5' or 'block3'
    Returns:
        pred_gamma_adj, Q_adj_arr, (D_arr, L_fit_arr, L_eval, info), targets_adj
    """
    # Raw scattering staircase (unwrapped)
    phi_unw = unwrap([phi_tau(t) for t in taus])
    Q_raw = -(phi_unw - phi_unw[0]) / math.pi

    # Analytic mismatch D(τ) = S(τ) − Q_raw(τ)
    D = Q_RS_arr - Q_raw

    # Fit a tiny smooth baseline L(τ) to D(τ)
    if baseline_method == 'poly5':
        L_fit, L_eval, info = fit_baseline_poly5(taus, D)
    else:
        L_fit, L_eval, info = fit_baseline_block3(taus, D, blocks=blocks)

    # Adjusted staircase and targets
    Q_adj = Q_raw + L_fit
    targets_adj = make_targets(Q_adj)

    # Bisection with frozen branch relative to left anchor
    phi0 = float(phi_unw[0])
    pred_tau_adj = []
    for target in targets_adj:
        idxs = crossing_brackets(taus, Q_adj, target)
        for i in idxs:
            a, b = float(taus[i]), float(taus[i + 1])
            base = float(phi_unw[i])
            def Q_eval(x):
                ph = phi_tau(x)  # principal
                k2 = round((base - ph) / (2 * math.pi))
                ph_aligned = ph + 2 * math.pi * k2
                return -(ph_aligned - phi0) / math.pi + L_eval(float(x))
            pred_tau_adj.append(refine_bisection(Q_eval, a, b, target, iters=24))
    pred_tau_adj = sorted(set([round(t, 12) for t in pred_tau_adj]))
    pred_gamma_adj = [2.0 * t for t in pred_tau_adj]
    return pred_gamma_adj, Q_adj, (D, L_fit, L_eval, info), targets_adj


# ------------------------ Main script ------------------------
def main():
    # ---- Parameters (safe demo defaults) ----
    TAU_MIN, TAU_MAX, NUM = 0.05, 35.0, 2200   # ~ first ~17 zeros
    BASELINE_METHOD = 'poly5'                  # 'poly5' or 'block3'
    BLOCKS = 6                                 # only used if 'block3'
    SAVE_CSV = True
    PLOT = True
    EVAL_WITH_TRUE_ZEROS = True                # for residuals only

    # ---- Grid ----
    taus = np.linspace(TAU_MIN, TAU_MAX, NUM)
    # ---- True zeros (for evaluation only) ----
    true_g = true_zeros_up_to(2 * TAU_MAX, 1000) if EVAL_WITH_TRUE_ZEROS else []

    # ---- Method A: RS quantizer ----
    pred_gamma_RS, Q_RS_arr, F_RS_arr, targets_RS = quantize_RS(taus)
    m_RS, pg_RS, tg_RS, err_RS, mean_RS, rms_RS, max_RS = stats_vs_true(pred_gamma_RS, true_g)

    print("=== METHOD A: Classical explicit (RS quantizer) ===")
    print(f"Window τ ∈ [{TAU_MIN}, {TAU_MAX}] with {NUM} samples (γ ≤ {2*TAU_MAX:.1f})")
    print(f"Matched zeros: {m_RS}")
    print(f"  mean error = {mean_RS:.6e}" if mean_RS is not None else "  mean error = None")
    print(f"  rms  error = {rms_RS:.6e}"  if rms_RS  is not None else "  rms  error = None")
    print(f"  max  error = {max_RS:.6e}"  if max_RS  is not None else "  max  error = None")

    # ---- Method B: RS‑Lock Eisenstein ----
    pred_gamma_adj, Q_adj, (D, L_fit, L_eval, info), targets_adj = quantize_RSLock(
        taus, Q_RS_arr, baseline_method=BASELINE_METHOD, blocks=BLOCKS
    )

    m_adj, pg_adj, tg_adj, err_adj, mean_adj, rms_adj, max_adj = stats_vs_true(pred_gamma_adj, true_g)
    rms_baseline = float(np.sqrt(np.mean((D - L_fit) ** 2)))

    print("\n=== METHOD B: RS‑Lock Eisenstein (scattering) ===")
    print(f"Baseline fit method: {info[0]}")
    if info[0] == 'poly5':
        coef = info[1]
        print("  poly5 coefficients (highest degree first):")
        for j, c in enumerate(coef):
            print(f"    c{5 - j} = {c:+.6e}")
    else:
        print(f"  blockwise deg‑3, blocks={len(info[1])}")
    print(f"  baseline fit RMS on D = {rms_baseline:.3e}")
    print(f"Matched zeros: {m_adj}")
    print(f"  mean error = {mean_adj:.6e}" if mean_adj is not None else "  mean error = None")
    print(f"  rms  error = {rms_adj:.6e}"  if rms_adj  is not None else "  rms  error = None")
    print(f"  max  error = {max_adj:.6e}"  if max_adj  is not None else "  max  error = None")

    # ---- Save CSV tables ----
    if SAVE_CSV and EVAL_WITH_TRUE_ZEROS:
        df_RS = pd.DataFrame({
            'n': np.arange(1, m_RS + 1),
            'gamma_true': tg_RS,
            'gamma_pred_RS': pg_RS,
            'error_RS': err_RS
        })
        df_adj = pd.DataFrame({
            'n': np.arange(1, m_adj + 1),
            'gamma_true': tg_adj,
            'gamma_pred_RSLock': pg_adj,
            'error_RSLock': err_adj
        })
        df_RS.to_csv('rs_quantizer_table.csv', index=False)
        df_adj.to_csv('rslock_eisenstein_table.csv', index=False)
        print("\nSaved tables:")
        print("  RS quantizer table: rs_quantizer_table.csv")
        print("  RS‑Lock Eisenstein table: rslock_eisenstein_table.csv")

    # ---- Plots ----
    if PLOT:
        # A1: RS staircase with half‑integers and zeros
        plt.figure(figsize=(8, 4.6))
        plt.plot(taus, Q_RS_arr, label='Q_RS(τ)')
        qmin, qmax = float(np.min(Q_RS_arr)), float(np.max(Q_RS_arr))
        h = math.ceil(qmin - 0.5)
        while h + 0.5 <= qmax:
            y = h + 0.5
            plt.axhline(y=y, linestyle='--', linewidth=0.9)
            h += 1
        if EVAL_WITH_TRUE_ZEROS:
            for t in [g / 2.0 for g in tg_RS]:
                plt.axvline(x=t, linewidth=0.7)
        for t in [g / 2.0 for g in pg_RS]:
            plt.axvline(x=t, linestyle=':', linewidth=0.9)
        plt.xlabel('τ'); plt.ylabel('Q_RS(τ)')
        plt.title('Method A — Classical RS quantizer')
        plt.legend(); plt.tight_layout(); plt.show()

        # A2: RS residuals
        if EVAL_WITH_TRUE_ZEROS and m_RS > 0:
            plt.figure(figsize=(8, 4.2))
            plt.plot(np.arange(1, m_RS + 1), err_RS, marker='o', ms=3, lw=1)
            plt.xlabel('n'); plt.ylabel('γ_pred − γ_true')
            plt.title('Method A — Residuals')
            plt.tight_layout(); plt.show()

        # B1: Baseline D vs L
        plt.figure(figsize=(8, 4.6))
        plt.plot(taus, D, label='D(τ) = S(τ) − Q_raw(τ)')
        plt.plot(taus, L_fit, label='L_fit(τ) (smooth baseline)')
        plt.xlabel('τ'); plt.ylabel('phase (π units)')
        plt.title('Method B — Baseline learning for Eisenstein phase')
        plt.legend(); plt.tight_layout(); plt.show()

        # B2: RS‑Lock staircase with half‑integers and zeros
        plt.figure(figsize=(8, 4.6))
        plt.plot(taus, Q_adj, label='Q_RSLock(τ)')
        qmin2, qmax2 = float(np.min(Q_adj)), float(np.max(Q_adj))
        h = math.ceil(qmin2 - 0.5)
        while h + 0.5 <= qmax2:
            y = h + 0.5
            plt.axhline(y=y, linestyle='--', linewidth=0.9)
            h += 1
        if EVAL_WITH_TRUE_ZEROS:
            for t in [g / 2.0 for g in tg_adj]:
                plt.axvline(x=t, linewidth=0.7)
        for t in [g / 2.0 for g in pg_adj]:
            plt.axvline(x=t, linestyle=':', linewidth=0.9)
        plt.xlabel('τ'); plt.ylabel('Q_RSLock(τ)')
        plt.title('Method B — RS‑Lock Eisenstein quantizer')
        plt.legend(); plt.tight_layout(); plt.show()

        # B3: RS‑Lock residuals
        if EVAL_WITH_TRUE_ZEROS and m_adj > 0:
            plt.figure(figsize=(8, 4.2))
            plt.plot(np.arange(1, m_adj + 1), err_adj, marker='o', ms=3, lw=1)
            plt.xlabel('n'); plt.ylabel('γ_pred − γ_true')
            plt.title('Method B — Residuals')
            plt.tight_layout(); plt.show()


if __name__ == "__main__":
    main()

