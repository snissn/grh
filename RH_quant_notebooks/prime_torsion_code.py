import os, math
import numpy as np
import matplotlib.pyplot as plt

# ---------- primes and pi(x)
def primes_upto(n: int):
    sieve = np.ones(n+1, dtype=bool); sieve[:2] = False
    for i in range(2, int(math.isqrt(n)) + 1):
        if sieve[i]:
            sieve[i*i:n+1:i] = False
    return np.nonzero(sieve)[0]

# ---------- torsion kernel and simple conv on uniform u-grid
def logistic_K(t):
    return 0.25 / (np.cosh(t/2.0)**2)

def conv_discrete_same(u, f, K):
    pad = len(K)//2
    return np.convolve(np.pad(f, pad, mode='edge'), K, 'same')[pad:-pad]

# ---------- data and grids
N = 10_000_000
u = np.linspace(2, math.log(N), 2001)          # u = log x (uniform)
x = np.exp(u)
pr = primes_upto(N)
pi = np.searchsorted(pr, x, side='right')
R = pi * u / x                                  # PNT-normalized ratio

# ---------- build kernel on u-grid step and normalize to unit mass
du = u[1] - u[0]
t = np.arange(-16.0, 16.0 + du/2, du)           # wide enough support
K = logistic_K(t)
K /= K.sum()                                    # preserves constant mode

# ---------- (A1) Overview: straightening plot
R_tor = conv_discrete_same(u, R, K)
plt.figure(figsize=(10,6))
plt.plot(u, R, alpha=.5, label='raw $R(u)$')
plt.plot(u, R_tor, lw=2, label='torsion-smoothed $R(u)$')
plt.axhline(1, ls='--', c='k')
plt.xlabel('$u=\\log x$'); plt.ylabel('$R(u)$')
plt.legend()
os.makedirs('figures', exist_ok=True)
plt.savefig('figures/fig_torsion_straighten.png', dpi=240, bbox_inches='tight')
plt.close()

# ---------- (A2-left) li-series comparison on interior window
def S_li(u, order=6):
    # 1 + 1/u + 2/u^2 + 6/u^3 + ...
    return sum(math.factorial(k) / (u**k) for k in range(order))

S = S_li(u, order=6)
mask = (u >= 8.0) & (u <= 10.5)
plt.figure(figsize=(10,6))
plt.plot(u[mask], R_tor[mask], lw=2, label='torsion-smoothed $R(u)$')
plt.plot(u[mask], S[mask], ls='--', label='$S_{\\mathrm{li}}(u)$ (to $5!/u^5$)')
plt.axhline(1, ls='--', c='k')
plt.xlabel('$u=\\log x$'); plt.ylabel('$R(u)$')
plt.legend(); plt.tight_layout()
plt.savefig('figures/torsion_li_comparison.png', dpi=240, bbox_inches='tight')
plt.close()

# ---------- (A2-right) derivative check on interior window
# Derivative via convolution with K' vs finite differences of (K*R)
# Use numerical derivative of K for robustness and enforce zero mean.
K1 = np.gradient(K, du)
K1 -= (K1.sum() * du) / (t[-1]-t[0])            # near-zero mean
R1_conv = conv_discrete_same(u, R, K1)
R1_fd = np.gradient(R_tor, du)

plt.figure(figsize=(10,5))
plt.plot(u[mask], R1_conv[mask], label='conv $K\' * R$')
plt.plot(u[mask], R1_fd[mask], ls='--', label='finite diff of $K*R$')
plt.xlabel('$u$'); plt.ylabel('first derivative')
plt.title('Derivative of torsion–smoothed ratio (interior window)')
plt.legend(); plt.tight_layout()
plt.savefig('figures/torsion_derivative_check.png', dpi=240, bbox_inches='tight')
plt.close()

