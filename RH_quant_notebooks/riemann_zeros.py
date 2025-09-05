# Generate a Riemann zeros file (first ~120 zeros on the critical line)
# We compute zeros of zeta(1/2 + i t) by finding real zeros of the
# Hardy Z-function Z(t) = e^{i theta(t)} zeta(1/2 + i t), which is real for real t.

from mpmath import mp

mp.dps = 60  # working precision (decimal digits)

def theta(t):
    # Riemann–Siegel theta function
    return mp.im(mp.loggamma(mp.mpf('0.25') + 0.5j*t)) - (t/2)*mp.log(mp.pi)

def Z(t):
    # Hardy's Z function; real-valued on R
    return mp.re(mp.e**(1j*theta(t)) * mp.zeta(mp.mpf('0.5') + 1j*t))

def find_zeros(n_target=120, t_start=0.1, t_max=320.0, step=0.1):
    zeros = []
    t = mp.mpf(t_start)
    f_prev = Z(t)
    t_prev = t
    while t <= t_max and len(zeros) < n_target:
        t_next = t + step
        f_next = Z(t_next)
        # Detect sign change or near-zero
        if f_prev == 0:
            root = t_prev
        elif f_prev * f_next <= 0:
            # bracketed root; refine with secant using the bracket
            try:
                root = mp.findroot(Z, (t_prev, t_next), tol=mp.mpf('1e-30'), maxsteps=50)
            except:  # fallback: mid-point + Newton
                mid = (t_prev + t_next) / 2
                try:
                    root = mp.findroot(Z, mid, tol=mp.mpf('1e-30'), maxsteps=50)
                except:
                    root = None
        else:
            root = None

        if root is not None:
            # de-duplicate if we saw this zero already (within tolerance)
            if not zeros or abs(root - zeros[-1]) > mp.mpf('1e-10'):
                zeros.append(root)
        # advance
        t_prev, f_prev = t_next, f_next
        t = t_next
    return zeros

zeros = find_zeros(n_target=12000, t_start=0.1, t_max=320.0, step=0.1)

# Prepare three common formats:
# 1) gamma-only (one value per line)
# 2) CSV with header gamma
# 3) CSV with header beta,gamma (beta fixed at 0.5)
# Also a small README about the format

# Choose a consistent number formatting
def fmt(x, digits=40):
    # mp.nstr returns a string with the requested significant digits
    return mp.nstr(x, n=digits)

gamma_only_path = "riemann_zeros_first_120_gamma.txt"
gamma_csv_path = "riemann_zeros_first_120_gamma.csv"
beta_gamma_csv_path = "riemann_zeros_first_120_beta_gamma.csv"
readme_path = "README_riemann_zeros_formats.txt"

with open(gamma_only_path, "w") as f:
    for g in zeros:
        f.write(fmt(g, 40) + "\n")

with open(gamma_csv_path, "w") as f:
    f.write("gamma\n")
    for g in zeros:
        f.write(fmt(g, 40) + "\n")

with open(beta_gamma_csv_path, "w") as f:
    f.write("beta,gamma\n")
    for g in zeros:
        f.write("0.5," + fmt(g, 40) + "\n")

with open(readme_path, "w") as f:
    f.write(
        "Riemann zeta nontrivial zeros (first ~120) on the critical line\n"
        "s = 1/2 + i*gamma, ordered by increasing gamma.\n\n"
        "Files provided:\n"
        "  1) riemann_zeros_first_120_gamma.txt\n"
        "     - Plain text, one gamma per line (no header).\n"
        "  2) riemann_zeros_first_120_gamma.csv\n"
        "     - CSV with header 'gamma'; one gamma per line.\n"
        "  3) riemann_zeros_first_120_beta_gamma.csv\n"
        "     - CSV with header 'beta,gamma'; beta fixed at 0.5 for each row.\n\n"
        "Precision: ~40 significant digits per entry. Computed with mpmath at 60 dps\n"
        "using Hardy's Z-function and bracketed root refinement.\n"
        "If your script expects a different format, the gamma-only file is the most\n"
        "widely supported. Rename as needed.\n"
    )

# Return the first few zeros so the user sees them inline as a quick check
zeros_preview = [fmt(z, 25) for z in zeros[:10]]
zeros[-1], len(zeros), zeros_preview

