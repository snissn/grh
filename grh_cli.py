#!/usr/bin/env python3
from decimal import Decimal
import argparse

from grh_verifier import verify_certificate
from zeta_gl1_cert import (
    build_zeta_gl1_certificate_rs,
    build_zeta_gl1_certificate_bl_evaluate,
    build_zeta_gl1_certificate_heat_evaluate,
    run_proof_of_grh_rs,
)


def main():
    parser = argparse.ArgumentParser(description="GRH Certificate Verifier CLI (Zeta GL(1) examples)")
    parser.add_argument("--mode", choices=["rs", "bl", "heat", "all", "max"], default="max",
                        help="Which mode to run: rs (RS-positivity), bl (BL evaluate), heat (heat evaluate), all (rs+bl+heat), max (most comprehensive)")
    parser.add_argument("--X", type=float, default=6.0, help="Band-limit X (for BL) or cutoff u for heat tail")
    parser.add_argument("--a", type=float, default=0.8, help="Heat test parameter a (>1/2)")
    parser.add_argument("--tau", type=float, default=2.0, help="Heat Gaussian tau")
    parser.add_argument("--proof", action="store_true", help="Run proof_of_GRH test net (RS baseline)")

    args = parser.parse_args()

    X = Decimal(str(args.X))
    a = Decimal(str(args.a))
    tau = Decimal(str(args.tau))

    def run_and_print(tag, cert):
        rep = verify_certificate(cert)
        print(f"\n[{tag}] -> {rep.result} (phase {rep.phase_passed})")
        print("details:", rep.details)

    if args.mode in ("rs", "all", "max"):
        cert = build_zeta_gl1_certificate_rs()
        if args.proof or args.mode == "max":
            run_proof_of_grh_rs()
        else:
            run_and_print("RS", cert)

    if args.mode in ("bl", "all", "max"):
        cert = build_zeta_gl1_certificate_bl_evaluate(X=X, a=a, tau=tau)
        run_and_print("BL-evaluate", cert)

    if args.mode in ("heat", "all", "max"):
        cert = build_zeta_gl1_certificate_heat_evaluate(X=X, a=a, tau=tau)
        run_and_print("Heat-evaluate", cert)


if __name__ == "__main__":
    main()

