
"""
Torsion–Mellin straightening of the prime-counting function.

Generates two figures used in the appendix:
  1) torsion_li_comparison.png
  2) torsion_derivative_check.png

Requirements: numpy, matplotlib
"""

import numpy as np
import math
import matplotlib.pyplot as plt

# --------------------------
# Primes and π(x)
# --------------------------
def primes_upto(n: int) -> np.ndarray:
    """Sieve of Eratosthenes up to n (inclusive). Returns sorted primes."""
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(math.isqrt(n)) + 1):
        if sieve[i]:
            sieve[i*i:n+1:i] = False
    return np.nonzero(sieve)[0]

# --------------------------
# Torsion kernel and derivatives
# K(u) = (1/4) * sech^2(u/2)
# --------------------------
def K_logistic(u: np.ndarray) -> np.ndarray:
    return 0.25 / (np.cosh(u/2.0) ** 2)

def K1_logistic(u: np.ndarray) -> np.ndarray:
    # derivative: K'(u) = - (1/4) sech^2(u/2) * tanh(u/2)
    return -0.25 * (1.0 / (np.cosh(u/2.0) ** 2)) * np.tanh(u/2.0)

def K2_logistic(u: np.ndarray) -> np.ndarray:
    # second derivative: K''(u) = (1/4) A B^2 - (1/8) A^2
    # where A = sech^2(u/2), B = tanh(u/2)
    A = 1.0 / (np.cosh(u/2.0) ** 2)
    B = np.tanh(u/2.0)
    return 0.25 * A * (B**2) - 0.125 * (A**2)

# --------------------------
# Convolution on a uniform u-grid via FFT with zero padding
# --------------------------
def convolve_uniform(u: np.ndarray, f: np.ndarray, kernel_fun, half_width: float) -> np.ndarray:
    """Continuous-style convolution (f * K)(u) with kernel sampled on the same grid step as u.
    u must be uniform; kernel support is [-half_width, half_width].
    """
    du = u[1] - u[0]
    n = len(u)
    # Build kernel on the same step
    t = np.arange(-half_width, half_width + du/2, du)
    K = kernel_fun(t)
    # normalise to area 1 for K (or zero mean for derivatives if desired)
    K = K / (K.sum() * du)
    m = len(K)

    # FFT-based zero-padded convolution
    L = 1 << (n + m - 1).bit_length()
    F = np.fft.rfft(f, L)
    G = np.fft.rfft(K, L)
    h = np.fft.irfft(F * G, L)[: n + m - 1] * du
    pad = (m - 1) // 2
    return h[pad : pad + n]

def convolve_uniform_with_kernel(u: np.ndarray, f: np.ndarray, t: np.ndarray, K: np.ndarray) -> np.ndarray:
    """Convolve f with a *given* kernel (already sampled on t with step = du).
    K is assumed to have unit area (sum(K)*du = 1).
    """
    du = u[1] - u[0]
    n, m = len(u), len(K)
    L = 1 << (n + m - 1).bit_length()
    F = np.fft.rfft(f, L)
    G = np.fft.rfft(K, L)
    h = np.fft.irfft(F * G, L)[: n + m - 1] * du
    pad = (m - 1) // 2
    return h[pad : pad + n]

# --------------------------
# li-series surrogate S_li(u) = 1 + 1/u + 2/u^2 + 6/u^3 + ...
# (truncated for plotting)
# --------------------------
def S_li_series(u: np.ndarray, order: int = 6) -> np.ndarray:
    terms = [math.factorial(k) / (u ** k) for k in range(order)]
    return np.sum(terms, axis=0)

# --------------------------
# Main script
# --------------------------
def main():
    # Parameters
    N_max = 10_000_000               # sieve limit for primes
    u_min, u_max, n_u = 2.0, math.log(N_max), 2001
    half_width = 16.0                # effective kernel half-width in u

    # Grid and data
    u = np.linspace(u_min, u_max, n_u)          # uniform in u = log x
    x = np.exp(u)
    primes = primes_upto(N_max)
    pi_vals = np.searchsorted(primes, x, side='right')
    R = pi_vals * u / x                           # PNT-normalised ratio

    # Torsion-smoothed curve and derivatives by convolution with K', K''
    R_tor = convolve_uniform(u, R, K_logistic, half_width)

    # Build K, K', K'' on the same step for derivative checks
    du = u[1] - u[0]
    t = np.arange(-half_width, half_width + du/2, du)
    K = K_logistic(t); K = K / (K.sum() * du)    # unit area
    K1 = K1_logistic(t); K1 -= (K1.sum() * du) / (t[-1]-t[0])  # enforce tiny zero-mean
    K2 = K2_logistic(t); K2 -= (K2.sum() * du) / (t[-1]-t[0])  # enforce tiny zero-mean

    R1_conv = convolve_uniform_with_kernel(u, R, t, K1)
    R2_conv = convolve_uniform_with_kernel(u, R, t, K2)

    # Finite-difference derivatives of the smoothed curve
    R1_fd = np.gradient(R_tor, du)
    R2_fd = np.gradient(R1_fd, du)

    # Interior window to avoid edge effects from finite x-range
    mask = (u >= 8.0) & (u <= 10.5)

    # Figure 1: comparison with li-series surrogate
    S = S_li_series(u, order=6)
    plt.figure(figsize=(10,6))
    plt.plot(u[mask], R_tor[mask], lw=2, label='torsion-smoothed $R(u)$')
    plt.plot(u[mask], S[mask], ls='--', label='$S_{\\mathrm{li}}(u)$ (to $5!/u^5$)')
    plt.axhline(1, color='k', lw=0.8, ls=':')
    plt.xlabel('$u=\\log x$'); plt.ylabel('$R(u)$')
    plt.title('Torsion–smoothed prime ratio vs. li-series surrogate (interior window)')
    plt.legend(); plt.tight_layout()
    plt.savefig('torsion_li_comparison.png', dpi=200); plt.close()

    # Figure 2: derivative check (quantified smoothness)
    plt.figure(figsize=(10,5))
    plt.plot(u[mask], R1_conv[mask], label='conv $K\' * R$')
    plt.plot(u[mask], R1_fd[mask], ls='--', label='finite diff of $K*R$')
    plt.xlabel('$u$'); plt.ylabel('first derivative')
    plt.title('Derivative of torsion–smoothed prime ratio')
    plt.legend(); plt.tight_layout()
    plt.savefig('torsion_derivative_check.png', dpi=200); plt.close()

    # Print tiny diagnostic on agreement of derivatives
    sup1 = np.max(np.abs(R1_conv[mask]-R1_fd[mask]))
    sup2 = np.max(np.abs(R2_conv[mask]-R2_fd[mask]))
    print({'sup_norm_diff_first_derivative': float(sup1),
           'sup_norm_diff_second_derivative': float(sup2)})

if __name__ == "__main__":
    main()
