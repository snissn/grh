#!/usr/bin/env python3
from decimal import Decimal
import argparse

from grh_verifier import verify_certificate, TestFunction
from verifier_trace import compute_trace
from reporting import plot_and_report
from zeta_gl1_cert import (
    build_zeta_gl1_certificate_rs,
    build_zeta_gl1_certificate_bl_evaluate,
    build_zeta_gl1_certificate_heat_evaluate,
    run_proof_of_grh_rs,
)
from lfunc_examples import (
    build_gl2_newform_11a1_bl_evaluate,
    build_gl2_11a1_heat_evaluate,
    build_gl3_sym2_11a1_bl_evaluate,
    build_dirichlet_bl_evaluate,
    build_dirichlet_heat_evaluate,
)


def main():
    parser = argparse.ArgumentParser(description="GRH Certificate Verifier CLI (Zeta GL(1) examples)")
    parser.add_argument("--mode", choices=["rs", "bl", "heat", "gl2", "gl2-heat", "sym2", "dir", "dir-heat", "all", "max"], default="max",
                        help="Which mode to run: rs (RS-positivity), bl (BL evaluate), heat (heat evaluate), all (rs+bl+heat), max (most comprehensive)")
    parser.add_argument("--ap-file", type=str, default=None, help="Path to a file with 'p a_p' lines for GL2 11a1 (optional)")
    parser.add_argument("--q", type=int, default=3, help="Dirichlet modulus (prime) for 'dir' modes")
    parser.add_argument("--r", type=int, default=1, help="Dirichlet character exponent r (chi(g)=exp(2πi r/(q-1)))")
    parser.add_argument("--X", type=float, default=6.0, help="Band-limit X (for BL) or cutoff u for heat tail")
    parser.add_argument("--a", type=float, default=0.8, help="Heat test parameter a (>1/2)")
    parser.add_argument("--tau", type=float, default=2.0, help="Heat Gaussian tau")
    parser.add_argument("--proof", action="store_true", help="Run proof_of_GRH test net (RS baseline)")
    parser.add_argument("--report", dest="report", action="store_true", help="Generate a detailed report with plots (default)")
    parser.add_argument("--no-report", dest="report", action="store_false", help="Disable report generation")
    parser.set_defaults(report=True)

    args = parser.parse_args()

    X = Decimal(str(args.X))
    a = Decimal(str(args.a))
    tau = Decimal(str(args.tau))

    def run_and_print(tag, cert):
        rep = verify_certificate(cert)
        print(f"\n[{tag}] -> {rep.result} (phase {rep.phase_passed})")
        print("details:", rep.details)
        if args.report:
            trace = compute_trace(cert)
            out_dir = plot_and_report(tag, trace)
            print(f"report written to: {out_dir}")

    if args.mode in ("rs", "all", "max"):
        cert = build_zeta_gl1_certificate_rs()
        # Reflect CLI cutoff into RS cert for context naming
        cert.band_limit_X = X
        if args.proof or args.mode == "max":
            cert.proof_mode = "proof_of_GRH"
            cert.test_net = [
                TestFunction(family="heat", a=Decimal("0.7")),
                TestFunction(family="heat", a=Decimal("0.8")),
                TestFunction(family="heat", a=Decimal("0.9")),
            ]
            tag = f"RS-proof_GL{cert.m}_{cert.K}_X{X}"
            run_and_print(tag, cert)
        else:
            tag = f"RS_GL{cert.m}_{cert.K}_X{X}"
            run_and_print(tag, cert)

    if args.mode in ("bl", "all", "max"):
        cert = build_zeta_gl1_certificate_bl_evaluate(X=X, a=a, tau=tau)
        tag = f"BL-evaluate_GL{cert.m}_{cert.K}_X{X}_a{a}_tau{tau}"
        run_and_print(tag, cert)

    if args.mode in ("heat", "all", "max"):
        cert = build_zeta_gl1_certificate_heat_evaluate(X=X, a=a, tau=tau)
        tag = f"Heat-evaluate_GL{cert.m}_{cert.K}_X{X}_a{a}_tau{tau}"
        run_and_print(tag, cert)

    # GL2 (holomorphic newform 11a1)
    if args.mode in ("gl2", "all", "max"):
        cert = build_gl2_newform_11a1_bl_evaluate(X=X, a=a, tau=tau, ap_file=args.ap_file)
        tag = f"GL2_11a1_BL_GL{cert.m}_{cert.K}_X{X}_a{a}_tau{tau}"
        run_and_print(tag, cert)
    if args.mode in ("gl2-heat", "all", "max"):
        cert = build_gl2_11a1_heat_evaluate(X=X, a=a, tau=tau, ap_file=args.ap_file)
        tag = f"GL2_11a1_Heat_GL{cert.m}_{cert.K}_X{X}_a{a}_tau{tau}"
        run_and_print(tag, cert)

    # GL3 Sym^2 of 11a1
    if args.mode in ("sym2", "all", "max"):
        use_amp = args.sym2_amplifier is not None and args.sym2_amplifier > 0
        amp_L = int(args.sym2_amplifier) if use_amp else 0
        cert = build_gl3_sym2_11a1_bl_evaluate(X=X, a=a, tau=tau, use_amplifier=use_amp, amp_L=amp_L)
        tag = f"GL3_Sym2_11a1_BL_GL{cert.m}_{cert.K}_X{X}_a{a}_tau{tau}"
        run_and_print(tag, cert)

    # Dirichlet L(s,chi) mod q
    if args.mode in ("dir", "all", "max"):
        cert = build_dirichlet_bl_evaluate(q=args.q, r=args.r, X=X, a=a, tau=tau)
        tag = f"Dirichlet_BL_mod{args.q}_r{args.r}_X{X}_a{a}_tau{tau}"
        run_and_print(tag, cert)
    if args.mode in ("dir-heat", "all", "max"):
        cert = build_dirichlet_heat_evaluate(q=args.q, r=args.r, X=X, a=a, tau=tau)
        tag = f"Dirichlet_Heat_mod{args.q}_r{args.r}_X{X}_a{a}_tau{tau}"
        run_and_print(tag, cert)


    parser.add_argument("--sym2-amplifier", type=int, default=None,
                        help="Enable Sym^2 amplifier with Fejér length L (prototype). If omitted or <=0, amplifier is off.")
if __name__ == "__main__":
    main()
