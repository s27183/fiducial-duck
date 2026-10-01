"""Unsteady incompressible Navier-Stokes in the same geometry, Taylor-Hood P2/P1.

Extension of solve_ns.py to time-dependent forcing: oscillatory, pulsatile or arbitrary
waveform-driven flow instead of a steady flow rate.

What is reused unchanged from the steady solver
-----------------------------------------------
* make_mesh.py geometry and its boundary tags ("inlet", "outlet", "top", "bottom")
* the P2/P1 bases, the viscous and divergence forms, and the Picard treatment of convection
* wall_shear() traction extraction -- applied per time step instead of once

What is new
-----------
* a velocity mass matrix and BDF2 time stepping (BDF1 startup when no second history state)
* a time-dependent inlet condition
* cycle-resolved wall-shear metrics: TAWSS, OSI

Verification
------------
The closed-form Womersley solution for oscillatory flow between parallel plates is an exact
solution of the full Navier-Stokes equations (the convective term vanishes identically for
unidirectional flow), so it verifies the time discretisation directly, exactly as the steady
parallel-plate solution verifies the steady solver. See benchmarks.py.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from skfem import (MeshTri, Basis, ElementVector, ElementTriP2, ElementTriP1,
                   BilinearForm, asm)
from skfem.helpers import grad, ddot, div, dot

from solve_ns import RHO, MU, W_CHAN, H_CHAN, vlap, divu, conv, wall_shear

NU = MU / RHO


@BilinearForm
def mass(u, v, w):
    return dot(u, v)


# ---------------------------------------------------------------- analytical solution
def womersley(z, t, H, dpdx_amp, omega, rho=RHO, mu=MU):
    """Oscillatory plane-channel flow driven by -dp/dx = Re{dpdx_amp * exp(i*omega*t)}.

    Exact solution of the unsteady Navier-Stokes equations for unidirectional flow between
    plates at z = 0 and z = H, with no-slip at both walls.
    """
    lam = np.sqrt(1j * omega / (mu / rho))
    zc = np.asarray(z) - 0.5 * H
    uhat = (dpdx_amp / (1j * omega * rho)) * (1.0 - np.cosh(lam * zc) / np.cosh(lam * 0.5 * H))
    return np.real(uhat * np.exp(1j * omega * t))


def womersley_number(H, omega, nu=NU):
    """Wo = (H/2) * sqrt(omega / nu); the unsteadiness parameter for a channel."""
    return 0.5 * H * np.sqrt(omega / nu)


def stokes_layer(omega, nu=NU):
    return np.sqrt(2.0 * nu / omega)


# ---------------------------------------------------------------- time stepping
def solve_transient(mesh_path, inlet_fn, dt, n_steps, H=H_CHAN, u0_fn=None,
                    linear=False, picard=2, on_step=None, verbose=False):
    """March the unsteady problem.

    inlet_fn(z, t) -> streamwise velocity on the inlet facet at time t.
    u0_fn(z)       -> initial streamwise velocity profile (uniform in x); zero if None.
    linear=True    -> drop convection and reuse one LU factorisation for every step. Valid
                      for the unidirectional verification case, where convection vanishes
                      identically; never for the device geometry.
    on_step(k, t, uh, ub, m) -> optional per-step callback for metric accumulation.
    """
    m = MeshTri.load(mesh_path)
    ub = Basis(m, ElementVector(ElementTriP2()), intorder=4)
    pb = Basis(m, ElementTriP1(), intorder=4)
    Nu, Np = ub.N, pb.N

    A = asm(vlap, ub)
    M = asm(mass, ub)
    B = asm(divu, ub, pb)
    Z = sp.csr_matrix((Np, Np))

    din = ub.get_dofs("inlet")
    d1, d2 = din.all(["u^1"]), din.all(["u^2"])
    D = np.unique(np.concatenate([din.all(),
                                  ub.get_dofs("top").all(),
                                  ub.get_dofs("bottom").all()]))
    z_in = ub.doflocs[1, d1]
    I = np.setdiff1d(np.arange(Nu + Np), D)      # free dofs: all pressures, non-Dirichlet velocities

    # initial state, and a second history state one step earlier so BDF2 can start immediately
    uh = np.zeros(Nu)
    if u0_fn is not None:
        i1 = ub.get_dofs().all(["u^1"]) if False else np.arange(0, Nu, 2)
        uh[i1] = u0_fn(ub.doflocs[1, i1], 0.0)
    u_prev = uh.copy()
    if u0_fn is not None:
        i1 = np.arange(0, Nu, 2)
        u_prev[i1] = u0_fn(ub.doflocs[1, i1], -dt)
        bdf2_from_start = True
    else:
        bdf2_from_start = False

    lu_cache = {}

    def build(coef, wind_vec):
        """System matrix for mass coefficient `coef` (= c/dt) and a given convecting field."""
        if linear:
            key = round(coef, 12)
            if key not in lu_cache:
                S = sp.bmat([[coef * RHO * M + MU * A, -B.T], [B, Z]], format="csr")
                lu_cache[key] = (splu(S[I][:, I].tocsc()), S)
            return lu_cache[key]
        C = asm(conv, ub, wind=ub.interpolate(wind_vec))
        S = sp.bmat([[coef * RHO * M + MU * A + RHO * C, -B.T], [B, Z]], format="csr")
        return splu(S[I][:, I].tocsc()), S

    for k in range(1, n_steps + 1):
        t = k * dt
        xd = np.zeros(Nu + Np)                       # Dirichlet data at this time
        xd[d1] = inlet_fn(z_in, t)
        xd[d2] = 0.0

        if k == 1 and not bdf2_from_start:           # implicit Euler startup
            coef, hist = 1.0 / dt, (M @ uh) / dt
        else:                                        # BDF2
            coef = 1.5 / dt
            hist = (2.0 * (M @ uh) - 0.5 * (M @ u_prev)) / dt
        rhs = np.zeros(Nu + Np)
        rhs[:Nu] = RHO * hist

        u_new = uh
        for _ in range(1 if linear else picard):
            lu, S = build(coef, u_new)
            r = rhs - S @ xd
            sol = xd.copy()
            sol[I] = lu.solve(r[I])
            u_new = sol[:Nu]
        u_prev, uh = uh, u_new
        if verbose and k % max(1, n_steps // 10) == 0:
            print(f"    step {k:4d}/{n_steps}  t={t:.5g}  |u|={np.linalg.norm(uh):.4e}", flush=True)
        if on_step is not None:
            on_step(k, t, uh, ub, m)

    return dict(mesh=m, ub=ub, pb=pb, u=uh, t=n_steps * dt)


# ---------------------------------------------------------------- cycle metrics
class ShearAccumulator:
    """Accumulate TAWSS and OSI on a wall over a cycle.

    TAWSS = mean of |tau| over the cycle.
    OSI   = (1 - |mean(tau_vector)| / mean(|tau|)) / 2, the standard oscillatory shear index;
            0 for unidirectional shear, 0.5 for shear that fully reverses.
    """

    def __init__(self, boundary="bottom"):
        self.boundary = boundary
        self.x = None
        self.sum_abs = None
        self.sum_signed = None
        self.n = 0

    def __call__(self, k, t, uh, ub, m):
        x, z, tau = wall_shear(dict(mesh=m, ub=ub, u=uh), self.boundary)
        if self.x is None:
            self.x, self.sum_abs, self.sum_signed = x, np.zeros_like(tau), np.zeros_like(tau)
        self.sum_abs += np.abs(tau)
        self.sum_signed += tau
        self.n += 1

    def result(self):
        tawss = self.sum_abs / self.n
        with np.errstate(divide="ignore", invalid="ignore"):
            osi = 0.5 * (1.0 - np.abs(self.sum_signed / self.n) / tawss)
        osi = np.where(np.isfinite(osi), osi, 0.0)
        return self.x, tawss, osi
