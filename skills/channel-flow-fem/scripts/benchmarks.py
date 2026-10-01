"""Closed-form benchmarks for the channel-flow solvers. Run these before trusting any result.

    python benchmarks.py steady                 # parallel-plate (Poiseuille) wall shear
    python benchmarks.py womersley              # oscillatory channel flow, time-step convergence
    python benchmarks.py all --out-dir results/

steady     Solves steady Navier-Stokes in a ridge-free 500 um channel at 6 mL/min and compares
           the wall shear stress with the exact value tau = 6 mu Q / (w h^2). Passes when the
           relative error is below 1e-8 (the discretisation contains the exact solution, so
           anything larger means a bug, not a resolution limit).

womersley  Drives the same channel with the closed-form Womersley profile at 1 Hz (Wo ~ 0.75)
           and 100 Hz (Wo ~ 7.5) and halves the time step twice. Passes when the 1 Hz error
           falls by a factor of ~4 per halving (BDF2 is second order). The 100 Hz case is
           reported, not graded: on the default mesh its error stalls near 3e-3 because the
           47 um Stokes layer is under-resolved -- a spatial floor, not a time-stepping one.
           Refining the wall spacing to ~3 um (about 15 elements per Stokes-layer thickness)
           brings it to ~5e-4.

Exit status is non-zero if a graded benchmark fails.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import argparse, json, tempfile
import numpy as np

H_BENCH = 500e-6


def control_mesh(workdir, lc_wall=8e-6):
    """Ridge-free 500 um channel, 2 mm long, used by every benchmark."""
    from make_mesh import build
    path = os.path.join(workdir, f"control_{lc_wall * 1e6:g}um.msh")
    if not os.path.exists(path):
        build(path, n_ridge=0, H=H_BENCH, L_IN=1.0e-3, L_OUT=1.0e-3, lc_bulk=30e-6, lc_wall=lc_wall)
    return path


def steady(workdir, q_ml_min=6.0, tol=1e-8):
    from solve_ns import solve_case, wall_shear, MU, W_CHAN
    res = solve_case(control_mesh(workdir), q_ml_min, H=H_BENCH, verbose=False)
    x, _, tau = wall_shear(res, "bottom")
    mid = (x > 1.0e-3) & (x < 1.6e-3)                  # away from inlet and outlet
    exact = 6 * MU * (q_ml_min * 1e-6 / 60.0) / (W_CHAN * H_BENCH ** 2)
    num = float(tau[mid].mean())
    rel = abs(num - exact) / exact
    out = dict(tau_numeric_Pa=num, tau_exact_Pa=float(exact), rel_error=float(rel),
               passed=bool(rel < tol))
    print(f"steady    tau = {num:.6f} Pa  exact {exact:.6f} Pa  rel err {rel:.1e}  "
          f"{'PASS' if out['passed'] else 'FAIL'}")
    return out


def womersley_case(workdir, f_hz, amp, steps=(25, 50, 100)):
    from solve_unsteady import womersley, womersley_number, stokes_layer, solve_transient
    om = 2 * np.pi * f_hz
    Wo, dl = womersley_number(H_BENCH, om), stokes_layer(om)
    prof = lambda z, t: womersley(z, t, H_BENCH, amp, om)
    mesh, rows = control_mesh(workdir), []
    for nst in steps:
        dt = (1.0 / f_hz) / nst
        r = solve_transient(mesh, prof, dt, nst, H=H_BENCH, u0_fn=prof, linear=True)
        ub, uh = r["ub"], r["u"]
        i1 = np.arange(0, ub.N, 2)                     # x-component dofs
        xs, zs = ub.doflocs[0, i1], ub.doflocs[1, i1]
        sel = (xs > 0.9e-3) & (xs < 1.1e-3)
        ex = prof(zs[sel], r["t"])
        rows.append((nst, dt, float(np.abs(uh[i1][sel] - ex).max() / np.abs(ex).max())))
    ratios = [rows[i][2] / rows[i + 1][2] for i in range(len(rows) - 1)]
    print(f"womersley {f_hz:5g} Hz  Wo={Wo:.2f}  Stokes layer {dl * 1e6:.1f} um  errors "
          + "  ".join(f"{e:.2e}" for *_, e in rows)
          + "  ratios " + " ".join(f"{x:.2f}" for x in ratios))
    return dict(Wo=float(Wo), stokes_layer_um=float(dl * 1e6),
                errors=[{"n_steps": n, "dt": d, "rel_err": e} for n, d, e in rows],
                order_ratios=[float(x) for x in ratios])


def womersley(workdir):
    out = {"1Hz": womersley_case(workdir, 1.0, 1.0e4),
           "100Hz": womersley_case(workdir, 100.0, 1.0e6)}
    lo = out["1Hz"]
    ok = all(3.5 < x < 5.5 for x in lo["order_ratios"]) and lo["errors"][-1]["rel_err"] < 1e-5
    print(f"womersley second-order convergence at Wo={lo['Wo']:.2f}: {'PASS' if ok else 'FAIL'}"
          f"  (100 Hz reported only -- see module docstring)")
    return out, ok


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("which", choices=("steady", "womersley", "all"))
    ap.add_argument("--out-dir", default=None, help="write steady_benchmark.json / "
                    "womersley_verification.json here")
    ap.add_argument("--workdir", default=None, help="where to keep meshes (default: temp dir)")
    a = ap.parse_args()
    wd = a.workdir or tempfile.mkdtemp(prefix="fiducial-bench-")
    os.makedirs(wd, exist_ok=True)
    failed = False
    if a.which in ("steady", "all"):
        s = steady(wd); failed |= not s["passed"]
        if a.out_dir:
            os.makedirs(a.out_dir, exist_ok=True)
            json.dump(s, open(os.path.join(a.out_dir, "steady_benchmark.json"), "w"), indent=1)
    if a.which in ("womersley", "all"):
        w, ok = womersley(wd); failed |= not ok
        if a.out_dir:
            os.makedirs(a.out_dir, exist_ok=True)
            json.dump(w, open(os.path.join(a.out_dir, "womersley_verification.json"), "w"), indent=1)
    sys.exit(1 if failed else 0)
